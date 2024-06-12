"""Implementation of PP2 and BP1 criteria."""

from typing import Optional, Tuple

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.annonars_range import ClinvarItem
from src.defs.annonars_variant import VariantResult
from src.defs.auto_acmg import PP2BP1
from src.defs.auto_pvs1 import SeqVarConsequence, SeqvarConsequenceMapping
from src.defs.exceptions import (
    AlgorithmError,
    AutoAcmgBaseException,
    InvalidAPIResposeError,
    MissingDataError,
)
from src.defs.genome_builds import GenomeRelease
from src.defs.mehari import TranscriptGene, TranscriptSeqvar
from src.defs.seqvar import SeqVar
from src.utils import SeqVarTranscriptsHelper


class AutoPP2BP1:
    """Class for automatic PP2 and BP1 prediction."""

    def __init__(
        self,
        seqvar: SeqVar,
        genome_release: GenomeRelease,
        variant_info: VariantResult,
        *,
        config: Config,
    ):
        #: Configuration to use.
        self.config = config or Config()
        #: Sequence variant to predict.
        self.seqvar = seqvar
        #: Genome release.
        self.genome_release = genome_release
        #: Variant information.
        self.variant_info = variant_info
        #: Annonars client.
        self.annonars_client = AnnonarsClient(api_base_url=config.api_base_url_annonars)
        #: Prediction result.
        self.prediction: PP2BP1 | None = None
        #: Comment to store the prediction explanation.
        self.comment: str = ""

    def _calculate_range(self, gene_transcript: TranscriptGene) -> Tuple[int, int]:
        """
        Calculate the range for the missense variants. The range is calculated based on the CDS
        start and end positions of the gene transcript and includes the entire transcript region.

        Args:
            gene_transcript: The gene transcript.

        Returns:
            Tuple[int, int]: The start and end positions of the range.

        Raises:
            MissingDataError: If the CDS start or end is not set.
        """
        if (
            not gene_transcript.genomeAlignments
            or not gene_transcript.genomeAlignments[0].cdsStart
            or not gene_transcript.genomeAlignments[0].cdsEnd
        ):
            raise MissingDataError("CDS start or end is not set.")
        return (
            gene_transcript.genomeAlignments[0].cdsStart,
            gene_transcript.genomeAlignments[0].cdsEnd,
        )

    def _is_missense(self, variant: ClinvarItem) -> bool:
        """
        Check if the variant is missense. The method retrieves information for the variant and
        checks if the consequence is missense.

        Note:
            If the consequence cannot be determined, the variant is not counted (False is
            returned).

        Args:
            variant: ClinVar information of the variant to check.

        Returns:
            bool: True if the variant is missense, False otherwise.
        """
        if variant.records and variant.records[0].sequenceLocation:
            # Get the seqvar information
            info = variant.records[0].sequenceLocation
            genome_release = (
                GenomeRelease.from_string(info.assembly)
                if info.assembly
                else self.seqvar.genome_release
            )
            chromosome = (
                info.chr.lower().replace("chromosome_", "") if info.chr else self.seqvar.chrom
            )
            start = info.start
            reference = info.referenceAlleleVcf
            alternate = info.alternateAlleleVcf

            if not start or not reference or not alternate:
                raise MissingDataError("Start, reference, or alternate allele is not set.")
            try:
                seqvar = SeqVar(
                    chrom=chromosome,
                    pos=start,
                    delete=reference,
                    insert=alternate,
                    genome_release=genome_release,
                )
                # Fetch transcript data and check the consequence
                seqvar_transcript_helper = SeqVarTranscriptsHelper(seqvar, config=self.config)
                seqvar_transcript_helper.initialize()
                (
                    _,
                    _,
                    _,
                    _,
                    consequence,
                ) = seqvar_transcript_helper.get_ts_info()
                return consequence == SeqVarConsequence.Missense
            except AutoAcmgBaseException:
                # Don't count the variant if the consequence cannot be determined
                return False
        else:
            return False

    def _get_missense_variants(
        self, seqvar: SeqVar, start_pos: int, end_pos: int
    ) -> Tuple[int, int, int]:
        """
        Counts pathogenic, benign, and total missense variants in the specified range.

        The method retrieves variants from the specified range and iterates through the ClinVar data
        of each variant to count the number of pathogenic variants, benign variants, and the total
        number of missense variants.

        Args:
            seqvar: The sequence variant being analyzed.
            start_pos: The start position of the range.
            end_pos: The end position of the range.

        Returns:
            Tuple[int, int, int]: The number of pathogenic variants, benign variants, and the total
            number of variants.

        Raises:
            InvalidAPIResposeError: If the API response is invalid or cannot be processed.
        """
        logger.debug("Counting missense variants in the range {} - {}.", start_pos, end_pos)
        if end_pos < start_pos:
            logger.error("End position is less than the start position.")
            logger.debug("Positions given: {} - {}", start_pos, end_pos)
            raise AlgorithmError("End position is less than the start position.")

        response = self.annonars_client.get_variant_from_range(seqvar, start_pos, end_pos)
        if response and response.result.clinvar:
            pathogenic_variants = []
            benign_variants = []
            missense_variants = []

            for v in response.result.clinvar:
                if (
                    v.records
                    and v.records[0].classifications
                    and v.records[0].classifications.germlineClassification
                ):
                    significance = v.records[0].classifications.germlineClassification.description
                    if not self._is_missense(v):
                        continue
                    missense_variants.append(v)
                    if significance in [
                        "Pathogenic",
                        "Likely pathogenic",
                    ]:
                        pathogenic_variants.append(v)
                    elif significance in [
                        "Benign",
                        "Likely benign",
                    ]:
                        benign_variants.append(v)

            return len(pathogenic_variants), len(benign_variants), len(missense_variants)
        else:
            logger.error("Failed to get variant from range. No ClinVar data.")
            raise InvalidAPIResposeError("Failed to get variant from range. No ClinVar data.")

    def predict(self) -> Tuple[Optional[PP2BP1], str]:
        """Predict PP2 and BP1 criteria."""
        self.prediction = PP2BP1()
        try:
            if self.seqvar.chrom == "MT":
                self.comment = "Variant is in mitochondrial DNA. PP2 and BP1 criteria are not met."
                logger.info("Variant is in mitochondrial DNA. PP2 and BP1 criteria are not met.")
                # skipped according to McCormick et al. (2020).
                self.prediction.PP2 = False
                self.prediction.BP1 = False
                return self.prediction, self.comment

            # Fetch transcript data
            self.comment = "Fetching transcript data. => \n"
            logger.debug("Fetching transcript data.")
            seqvar_transcript_helper = SeqVarTranscriptsHelper(self.seqvar, config=self.config)
            seqvar_transcript_helper.initialize()
            (
                _,
                gene_transcript,
                _,
                _,
                consequence,
            ) = seqvar_transcript_helper.get_ts_info()

            if not gene_transcript or consequence == SeqVarConsequence.NotSet:
                logger.error("Transcript data is not set. Cannot initialize the PVS1 class.")
                raise MissingDataError(
                    "Transcript data is not fully set. Cannot initialize the PVS1 class."
                )
            if consequence != SeqVarConsequence.Missense:
                self.comment += (
                    f"Consequence is not missense. Consequence: {consequence}. "
                    "PP2 and BP1 criteria are not met."
                )
                logger.info(
                    "Consequence is not missense: {}. PP2 and BP1 criteria are not met.",
                    consequence,
                )
                self.prediction.PP2 = False
                self.prediction.BP1 = False
                return self.prediction, self.comment

            start_pos, end_pos = self._calculate_range(gene_transcript)
            self.comment += f"Count missense variants on range: {start_pos} - {end_pos}. => \n"
            logger.debug("Count missense variants on range: {} - {}.", start_pos, end_pos)
            pathogenic_count, benign_count, total_count = self._get_missense_variants(
                self.seqvar, start_pos, end_pos
            )
            pathogenic_ratio = pathogenic_count / total_count
            benign_ratio = benign_count / total_count
            self.comment += (
                f"Pathogenic missense variants: {pathogenic_count}, "
                f"Benign missense variants: {benign_count}, "
                f"Total missense variants: {total_count}. => \n"
                f"Pathogenic ratio: {pathogenic_ratio}, Benign ratio: {benign_ratio}. => \n"
            )
            logger.debug(
                "Pathogenic missense variants: {}, Benign missense variants: {}, "
                "Total missense variants: {}",
                pathogenic_count,
                benign_count,
                total_count,
            )

            if pathogenic_ratio > 0.808:
                self.comment += "Pathogenic ratio is greater than 0.808. PP2 is met."
                logger.info("Pathogenic ratio is greater than 0.808. PP2 is met.")
                self.prediction.PP2 = True
            else:
                self.comment += "Pathogenic ratio is less than 0.808. PP2 is not met."
                logger.info("Pathogenic ratio is less than 0.808. PP2 is not met.")
            if benign_ratio > 0.569:
                self.comment += "Benign ratio is greater than 0.569. BP1 is met."
                logger.info("Benign ratio is greater than 0.569. BP1 is met.")
                self.prediction.BP1 = True
            else:
                self.comment += "Benign ratio is less than 0.569. BP1 is not met."
                logger.info("Benign ratio is less than 0.569. BP1 is not met.")

        except AutoAcmgBaseException as e:
            self.comment += f"Error occurred during PP2 and BP1 prediction. Error: {e}"
            logger.error("Failed to predict PP2 and BP1 criteria. Error: {}", e)
            self.prediction = None

        return self.prediction, self.comment
