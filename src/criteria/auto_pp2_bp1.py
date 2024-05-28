"""Implementation of PP2 and BP1 criteria."""

from typing import Optional, Tuple

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.annonars_variant import VariantResult
from src.defs.auto_acmg import PP2BP1
from src.defs.auto_pvs1 import SeqVarConsequence, SeqvarConsequenceMapping
from src.defs.exceptions import AlgorithmError, InvalidAPIResposeError, MissingDataError
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
        """Initialize the class."""
        self.seqvar = seqvar
        self.genome_release = genome_release
        self.variant_info = variant_info
        self.config = config
        self.annonars_client = AnnonarsClient(api_base_url=config.api_base_url_annonars)
        self.prediction: PP2BP1 | None = None

    @staticmethod
    def _get_consequence(seqvar_transcript: TranscriptSeqvar | None) -> SeqVarConsequence:
        """Get the consequence of the sequence variant.

        Args:
            seqvar_transcript: The sequence variant transcript.

        Returns:
            SeqVarConsequence: The consequence of the sequence variant.
        """
        if not seqvar_transcript:
            return SeqVarConsequence.NotSet
        else:
            for consequence in seqvar_transcript.consequences:
                if consequence in SeqvarConsequenceMapping:
                    return SeqvarConsequenceMapping[consequence]
            return SeqVarConsequence.NotSet

    @staticmethod
    def _calculate_range(gene_transcript: TranscriptGene) -> Tuple[int, int]:
        """
        Calculate the range for the missense variants.

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

    def _get_missense_variants(
        self, seqvar: SeqVar, start_pos: int, end_pos: int
    ) -> Tuple[int, int, int]:
        """
        Counts pathogenic, benign, and total missense variants in the specified range.

        The method retrieves variants from the specified range and iterates through the ClinVar data
        of each variant to count the number of pathogenic variants, benign variants, and the total number
        of missense variants.

        Args:
            seqvar: The sequence variant being analyzed.
            start_pos: The start position of the range.
            end_pos: The end position of the range.

        Returns:
            Tuple[int, int, int]: The number of pathogenic variants, benign variants, and the total number of variants.

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

            # TODO: Check if the variant is missense
            for variant in response.result.clinvar:
                if (
                    variant.referenceAssertions
                    and variant.referenceAssertions[0].clinicalSignificance
                ):
                    significance = variant.referenceAssertions[0].clinicalSignificance
                    missense_variants.append(variant)
                    if significance in [
                        "CLINICAL_SIGNIFICANCE_LIKELY_PATHOGENIC",
                        "CLINICAL_SIGNIFICANCE_PATHOGENIC",
                    ]:
                        pathogenic_variants.append(variant)
                    elif significance in [
                        "CLINICAL_SIGNIFICANCE_BENIGN",
                        "CLINICAL_SIGNIFICANCE_LIKELY_BENIGN",
                    ]:
                        benign_variants.append(variant)

            logger.debug(
                "Pathogenic missense variants: {}, Benign missense variants: {}, Total missense variants: {}",
                len(pathogenic_variants),
                len(benign_variants),
                len(missense_variants),
            )
            return len(pathogenic_variants), len(benign_variants), len(missense_variants)
        else:
            logger.error("Failed to get variant from range. No ClinVar data.")
            raise InvalidAPIResposeError("Failed to get variant from range. No ClinVar data.")

    def predict(self) -> Optional[PP2BP1]:
        """Predict PP2 and BP1 criteria."""
        self.prediction = PP2BP1()
        try:
            if self.seqvar.chrom == "MT":
                # skipped according to McCormick et al. (2020).
                self.prediction.PP2 = False
                self.prediction.BP1 = False
                return self.prediction

            # Fetch transcript data
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
                self.prediction.PP2 = False
                self.prediction.BP1 = False
                return self.prediction

            start_pos, end_pos = self._calculate_range(gene_transcript)
            pathogenic_count, benign_count, total_count = self._get_missense_variants(
                self.seqvar, start_pos, end_pos
            )
            pathogenic_ratio = pathogenic_count / total_count
            benign_ratio = benign_count / total_count

            if pathogenic_ratio > 0.808:
                self.prediction.PP2 = True
            if benign_ratio > 0.569:
                self.prediction.BP1 = True

        except:
            logger.error("Failed to predict PP2 and BP1 criteria.")
            self.prediction = None

        return self.prediction
