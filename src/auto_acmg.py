"""Implementations of the PVS1 algorithm."""

from typing import List, Optional, Union

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.criteria.auto_bp7 import AutoBP7
from src.criteria.auto_criteria import AutoACMGCriteria
from src.criteria.auto_pm1 import AutoPM1
from src.criteria.auto_pm2_ba1_bs1_bs2 import AutoPM2BA1BS1BS2
from src.criteria.auto_pm4_bp3 import AutoPM4BP3
from src.criteria.auto_pp2_bp1 import AutoPP2BP1
from src.criteria.auto_pp3_bp4 import AutoPP3BP4
from src.criteria.auto_ps1_pm5 import AutoPS1PM5
from src.defs.annonars_variant import VariantResult
from src.defs.auto_acmg import ACMGPrediction, AutoACMGPrediction, AutoACMGResult
from src.defs.auto_pvs1 import PVS1Prediction, PVS1PredictionPathMapping, PVS1PredictionStrucVarPath
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, ParseError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar, SeqVarResolver
from src.defs.strucvar import StrucVar, StrucVarResolver
from src.pvs1.auto_pvs1 import AutoPVS1
from src.utils import SeqVarTranscriptsHelper

#: Pathogenic PVS1 predictions of sequence variants.
PVS1_POSITIVE_SEQVAR_PREDICTIONS = [
    PVS1Prediction.PVS1,
    PVS1Prediction.PVS1_Strong,
    PVS1Prediction.PVS1_Moderate,
    PVS1Prediction.PVS1_Supporting,
]

#: Criteria, which have no automatic prediction yet.
NOT_IMPLEMENTED_CRITERIA: List[str] = [
    "PS2",  # De novo (both parents not tested)
    "PS3",  # Well-established in vitro or in vivo functional studies
    "PS4",  # Increased prevalence in affected individuals vs. controls
    "PM3",  # Recessive disorder (zygosity unknown)
    "PM6",  # Assumed de novo (both parents not tested)
    "PP1",  # Cosegregation with disease in multiple affected family members
    "PP4",  # Patient's phenotype or family history specificity
    "BS3",  # Well-established in vitro or in vivo functional studies
    "BS4",  # Lack of segregation in affected family members
    "BP2",  # Observed in healthy individuals (zygosity unknown)
    "BP5",  # Case with an alternate molecular basis for disease
]


class AutoACMG(AutoPS1PM5, AutoPM1, AutoPM2BA1BS1BS2, AutoPM4BP3, AutoPP2BP1, AutoPP3BP4, AutoBP7):
    """Class for predicting ACMG criteria.

    This class handles both sequence variants and structural variants to determine their potential
    impact under the various criteria of the ACMG guidelines for variant classification. Currently
    it only implements the PVS1 criterion (not finished yet).

    Attributes:
        variant_name (str): The name or identifier of the variant being analyzed.
        genome_release (GenomeRelease): The genome release version, defaults to GRCh38.
    """

    def __init__(
        self,
        variant_name: str,
        genome_release: GenomeRelease = GenomeRelease.GRCh38,
        *,
        config: Optional[Config] = None,
    ):
        """Initializes the AutoACMG with the specified variant and genome release.

        Args:
            variant_name: The name or identifier of the variant.
            genome_release (Optional): The genome release version, such as GRCh38 or GRCh37.
        """
        #: Configuration to use.
        self.config = config or Config()
        #: Annonars client.
        self.annonars_client: AnnonarsClient = AnnonarsClient(
            api_base_url=self.config.api_base_url_annonars
        )
        #: The name or identifier of the variant.
        self.variant_name = variant_name
        #: The genome release version.
        self.genome_release = genome_release
        #: The resolved sequence variant.
        self.seqvar: Optional[SeqVar] = None
        #: Prediction result.
        self.result: AutoACMGResult = AutoACMGResult()

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

    def _convert_score_val(self, score_value: Optional[Union[str, float, int]]) -> Optional[float]:
        """
        Convert score value to float.

        Since the score values can be represented as strings (with ";" as separator), we pick the
        maximum value that is not empty ("."). If the value is already numeric, we return it as is.

        Args:
            score_value (Optional[Union[str, float, int]]): Score value to convert.

        Returns:
            Optional[float]: Converted score value.

        Raises:
            AlgorithmError: If the score value cannot be converted to float.
        """
        if score_value is None:
            return None
        if isinstance(score_value, (float, int)):
            return float(score_value)
        try:
            return max(float(score) for score in score_value.split(";") if score != ".")
        except ValueError as e:
            logger.error("Failed to convert score value to float. Error: {}", e)
            raise AlgorithmError("Failed to convert score value to float.") from e

    def resolve_variant(self) -> Optional[SeqVar]:
        """Attempts to resolve the specified variant as either a sequence or structural variant.

        This method first tries to resolve the variant as a sequence variant. If it fails, it then
        attempts to resolve it as a structural variant.

        Returns:
            SeqVar or None: The resolved variant object or None if resolution fails.

        Raises:
            Exception: Specific exceptions are caught and logged, but generic exceptions may be
            raised if both resolutions fail.
        """
        logger.debug("Resolving variant: {}", self.variant_name)
        try:
            seqvar_resolver = SeqVarResolver(config=self.config)
            seqvar: SeqVar = seqvar_resolver.resolve_seqvar(self.variant_name, self.genome_release)
            logger.debug("Resolved sequence variant: {}", seqvar)
            return seqvar
        except ParseError:
            logger.exception("Failed to resolve sequence variant.")
            return None
        except AutoAcmgBaseException as e:
            logger.error("An unexpected error occurred: {}", e)
            return None

    def parse_data(self, seqvar: SeqVar) -> None:
        """Parses the data for the prediction."""
        # Mehari data
        ts_helper = SeqVarTranscriptsHelper(seqvar)
        ts_helper.initialize()
        seqvar_transcript, gene_transcript, _, _, _ = ts_helper.get_ts_info()
        if not seqvar_transcript or not gene_transcript:
            raise AutoAcmgBaseException("Transcript information is missing.")

        self.result.data.consequence.mehari = seqvar_transcript.consequences
        self.result.data.gene_symbol = seqvar_transcript.gene_symbol
        self.result.data.hgnc_id = seqvar_transcript.gene_id
        self.result.data.transcript_id = seqvar_transcript.feature_id
        self.result.data.transcript_tags = seqvar_transcript.feature_tag
        self.result.data.cds_start = gene_transcript.genomeAlignments[0].cdsStart
        self.result.data.cds_end = gene_transcript.genomeAlignments[0].cdsEnd
        self.result.data.strand = gene_transcript.genomeAlignments[0].strand
        self.result.data.exons = gene_transcript.genomeAlignments[0].exons

        # Annonars data
        variant_info = self._get_var_info(seqvar)
        if not variant_info:
            logger.error("Failed to get variant information.")
            raise AutoAcmgBaseException("Failed to get variant information.")

        if cadd := variant_info.cadd:
            self.result.data.consequence.cadd = cadd.ConsDetail or ""
            self.result.data.consequence.cadd_consequence = cadd.Consequence or ""
            self.result.data.scores.cadd.phyloP100 = cadd.verPhyloP
            self.result.data.scores.cadd.spliceAI_acceptor_gain = cadd.SpliceAI_acc_gain
            self.result.data.scores.cadd.spliceAI_acceptor_loss = cadd.SpliceAI_acc_loss
            self.result.data.scores.cadd.spliceAI_donor_gain = cadd.SpliceAI_don_gain
            self.result.data.scores.cadd.spliceAI_donor_loss = cadd.SpliceAI_don_loss
            self.result.data.scores.cadd.ada = cadd.dbscSNV_ada
            self.result.data.scores.cadd.rf = cadd.dbscSNV_rf
        if dbsnfp := variant_info.dbnsfp:
            self.result.data.pHGVS = dbsnfp.HGVSp_VEP or ""
            self.result.data.scores.dbnsfp.alpha_missense = self._convert_score_val(
                dbsnfp.AlphaMissense_score
            )
            self.result.data.scores.dbnsfp.metaRNN = self._convert_score_val(dbsnfp.MetaRNN_score)
            self.result.data.scores.dbnsfp.bayesDel_noAF = self._convert_score_val(
                dbsnfp.BayesDel_noAF_score
            )
            self.result.data.scores.dbnsfp.revel = self._convert_score_val(dbsnfp.REVEL_score)
            self.result.data.scores.dbnsfp.phyloP100 = self._convert_score_val(
                dbsnfp.phyloP100way_vertebrate
            )
        if dbscsnv := variant_info.dbscsnv:
            self.result.data.scores.dbscsnv.ada = dbscsnv.ada_score
            self.result.data.scores.dbscsnv.rf = dbscsnv.rf_score

        # Thresholds
        pass

    def predict(self) -> Optional[AutoACMGResult]:
        """
        Predict ACMG criteria for the specified variant.

        This method first resolves the variant, then predicts the PVS1 criterion for sequence
        variants and other ACMG criteria.

        Note:
            The method can resolve both sequence and structural variants, but currently only
            sequence variants are supported for ACMG criteria prediction.

        Raises:
            Exception: Specific exceptions are caught and logged, but generic exceptions may be
            raised if the prediction fails.
        """
        logger.info("Predicting ACMG criteria for variant: {}", self.variant_name)
        self.seqvar = self.resolve_variant()
        if not self.seqvar:
            logger.error("Unable to make a prediction for the variant: {}", self.variant_name)
            return None

        if isinstance(self.seqvar, SeqVar):
            # ====== Setup the data ======
            self.parse_data(self.seqvar)

            # ====== Predict ======
            # PP5 and BP6 criteria are depricated
            logger.warning("Note, that PP5 and BP6 criteria are depricated and not predicted.")
            # Not implemented criteria
            logger.warning(
                "Some criteria are not implemented yet: {}",
                NOT_IMPLEMENTED_CRITERIA,
            )

            # PVS1
            pass

            # PS1 and PM5
            self.result.criteria.ps1, self.result.criteria.pm5 = self.predict_ps1pm5(
                self.seqvar, self.result.data
            )

            # PM1
            self.result.criteria.pm1 = self.predict_pm1(self.seqvar, self.result.data)

            # PM2, BA1, BS1 and BS2
            (
                self.result.criteria.pm2,
                self.result.criteria.ba1,
                self.result.criteria.bs1,
                self.result.criteria.bs2,
            ) = self.predict_pm2ba1bs1bs2(self.seqvar, self.result.data)

            # PM4 and BP3
            self.result.criteria.pm4, self.result.criteria.bp3 = self.predict_pm4bp3(
                self.seqvar, self.result.data
            )

            # PP2 and BP1
            self.result.criteria.pp2, self.result.criteria.bp1 = self.predict_pp2bp1(
                self.seqvar, self.result.data
            )

            # PP3 and BP4
            self.result.criteria.pp3, self.result.criteria.bp4 = self.predict_pp3bp4(
                self.seqvar, self.result.data
            )

            # BP7
            self.result.criteria.bp7 = self.predict_bp7(self.seqvar, self.result.data)

            logger.info("ACMG criteria prediction completed.")
            return self.result

        else:
            logger.info("Structural variants are not supported for ACMG criteria prediction yet.")
            return None
