"""Implementation of BP7 criteria."""

from typing import List, Optional, Tuple

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.annonars_variant import VariantResult
from src.defs.auto_acmg import BP7
from src.defs.exceptions import AutoAcmgBaseException, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.mehari import Exon
from src.defs.seqvar import SeqVar
from src.utils import SeqVarTranscriptsHelper, SplicingPrediction


class AutoBP7:
    """Class for automatic BP7 prediction."""

    def __init__(
        self,
        seqvar: SeqVar,
        variant_info: VariantResult,
        *,
        config: Optional[Config] = None,
    ):
        #: Configuration to use.
        self.config: Config = config or Config()
        #: Sequence variant to predict.
        self.seqvar: SeqVar = seqvar
        #: Variant information.
        self.variant_info: VariantResult = variant_info
        #: Annonars client.
        self.annonars_client: AnnonarsClient = AnnonarsClient(
            api_base_url=self.config.api_base_url_annonars
        )
        #: Prediction result.
        self.prediction: Optional[BP7] = None
        #: Comment to store the prediction explanation.
        self.comment: str = ""

    def _get_var_info(self, seqvar: SeqVar) -> Optional[VariantResult]:
        """
        Get variant information from Annonars.

        Args:
            seqvar: The variant to get information for.

        Returns:
            Optional[VariantResult]: The variant information. None if the information is not
            available.
        """
        try:
            logger.debug("Getting variant information for {}.", seqvar)
            return self.annonars_client.get_variant_info(seqvar).result
        except AutoAcmgBaseException as e:
            logger.error("Failed to get variant information. Error: {}", e)
            return None

    def _check_proximity_to_pathogenic_vars(self, seqvar: SeqVar) -> bool:
        """
        Check for pathogenic variants +/- 2bp of the position in ClinVar.

        Check if there are pathogenic variants within 2bp of the position in ClinVar.

        Args:
            seqvar: The variant to check.

        Returns:
            bool: True if pathogenic variants are found, False otherwise.
        """
        response = self.annonars_client.get_variant_from_range(
            seqvar, seqvar.pos - 2, seqvar.pos + 2
        )
        if response and response.result.clinvar:
            pathogenic_variants = [
                v
                for v in response.result.clinvar
                if v.records
                and v.records[0].classifications
                and v.records[0].classifications.germlineClassification
                and v.records[0].classifications.germlineClassification.description
                in [
                    "Pathogenic",
                    "Likely pathogenic",
                ]
            ]
            return len(pathogenic_variants) > 0
        return False

    def _check_proximity_to_ss(self, seqvar: SeqVar) -> bool:
        """
        Check if the variant is closer than 2bp to a splice site.

        Check if the variant is located within 2bp of a splice site.

        Args:
            seqvar: The variant to check.

        Returns:
            bool: True if the variant is within 2bp of a splice site, False otherwise.
        """
        # Fetch transcript data
        seqvar_transcript_helper = SeqVarTranscriptsHelper(seqvar, config=self.config)
        seqvar_transcript_helper.initialize()
        (
            _,
            gene_transcript,
            _,
            _,
            _,
        ) = seqvar_transcript_helper.get_ts_info()
        if (
            not gene_transcript
            or not gene_transcript.genomeAlignments
            or not gene_transcript.genomeAlignments[0].exons
        ):
            logger.warning("No exons found for the transcript.")
            return False

        # Check if the variant is within 2bp of a start or end of an exon
        for exon in gene_transcript.genomeAlignments[0].exons:
            if abs(seqvar.pos - exon.altCdsStartI) <= 2 or abs(exon.altCdsEndI - seqvar.pos) <= 2:
                return True
        return False

    def _pred_spliceai(self, seqvar: SeqVar) -> bool:
        """
        Predict splice site alterations using SpliceAI.

        Fetch the SpliceAI data from CADD and predict if the variant is a splice site alteration.
        If any of SpliceAI scores are greater than 0, the variant is considered a splice site
        alteration.

        Args:
            seqvar: The variant to check.

        Returns:
            bool: True if the variant is a splice site alteration, False otherwise.
        """
        variant_info = self._get_var_info(seqvar)
        if not variant_info or not variant_info.cadd:
            raise MissingDataError("Missing CADD data for variant.")
        acc_gain = variant_info.cadd.SpliceAI_acc_gain
        acc_loss = variant_info.cadd.SpliceAI_acc_loss
        don_gain = variant_info.cadd.SpliceAI_don_gain
        don_loss = variant_info.cadd.SpliceAI_don_loss
        if (
            (acc_gain and acc_gain > 0)
            or (acc_loss and acc_loss > 0.5)
            or (don_gain and don_gain > 0)
            or (don_loss and don_loss > 0.5)
        ):
            return True
        return False

    def predict(self) -> Tuple[Optional[BP7], str]:
        """Predict BP7 criterion."""
        self.prediction = BP7()
        try:
            if self.seqvar.chrom == "MT":
                # skipped according to McCormick et al. (2020).
                self.comment = "Variant is in the mitochondrial genome. BP7 is not met."
                logger.debug("Variant is in the mitochondrial genome. BP7 is not met.")
                self.prediction.BP7 = False
                return self.prediction, self.comment

            self.comment = "Checking for pathogenic variants in the range of 2bp. => \n"
            logger.debug("Checking for pathogenic variants in the range of 2bp.")
            if self._check_proximity_to_pathogenic_vars(self.seqvar):
                self.comment += "Found pathogenic variants in the range of 2bp. BP7 is not met."
                logger.debug("Found pathogenic variants in the range of 2bp. BP7 is not met.")
                self.prediction.BP7 = False
                return self.prediction, self.comment
            else:
                self.comment += "No pathogenic variants found in the range of 2bp. => \n"
                logger.debug("No pathogenic variants found in the range of 2bp. => \n")

            self.comment += "Checking for proximity to splice site. => \n"
            logger.debug("Checking for proximity to splice site.")
            if self._check_proximity_to_ss(self.seqvar):
                self.comment += "Variant is within 2bp of a splice site. BP7 is not met."
                logger.debug("Variant is within 2bp of a splice site. BP7 is not met.")
                self.prediction.BP7 = False
                return self.prediction, self.comment
            else:
                self.comment += "Variant is not within 2bp of a splice site. => \n"
                logger.debug("Variant is not within 2bp of a splice site. => \n")

            self.comment += "Predicting splice site alterations using SpliceAI. => \n"
            logger.debug("Predicting splice site alterations using SpliceAI.")
            if self._pred_spliceai(self.seqvar):
                self.comment += "Variant is a splice site alteration. BP7 is not met."
                logger.debug("Variant is a splice site alteration. BP7 is not met.")
                self.prediction.BP7 = False
            else:
                self.comment += "Variant is not a splice site alteration. BP7 is met."
                logger.debug("Variant is not a splice site alteration. BP7 is met.")
                self.prediction.BP7 = True
        except AutoAcmgBaseException as e:
            self.comment += f"Failed to predict BP7 criterion. Error: {e}"
            logger.error("Failed to predict BP7 criterion. Error: {}", e)
            self.prediction = None

        return self.prediction, self.comment
