"""Implementation of PM1 criteria."""

from typing import Optional, Tuple

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.annonars_range import VariantInfo
from src.defs.auto_acmg import PM1
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, InvalidAPIResposeError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


class AutoPM1:
    """Class for automatic PM1 prediction."""

    def __init__(
        self,
        seqvar: SeqVar,
        genome_release: GenomeRelease,
        variant_info: VariantInfo,
        *,
        config: Config,
    ):
        """Initialize the class."""
        self.seqvar = seqvar
        self.genome_release = genome_release
        self.variant_info = variant_info
        self.annonars_client = AnnonarsClient(api_base_url=config.api_base_url_annonars)
        self.config = config
        self.prediction: PM1 = PM1()

    def _count_pathogenic_variants(
        self, seqvar: SeqVar, start_pos: int, end_pos: int
    ) -> Tuple[int, int]:
        """
        Counts pathogenic variants in the specified range.

        The method retrieves variants from the specified range and iterates through the ClinVar data
        of each variant to count the number of pathogenic variants and the total number of variants.

        Args:
            seqvar: The sequence variant being analyzed.
            start_pos: The start position of the range.
            end_pos: The end position of the range.

        Returns:
            Tuple[int, int]: The number of pathogenic variants and the total number of variants.

        Raises:
            InvalidAPIResposeError: If the API response is invalid or cannot be processed.
        """
        logger.debug("Counting pathogenic variants in the range {} - {}.", start_pos, end_pos)
        if end_pos < start_pos:
            logger.error("End position is less than the start position.")
            logger.debug("Positions given: {} - {}", start_pos, end_pos)
            raise AlgorithmError("End position is less than the start position.")

        response = self.annonars_client.get_variant_from_range(seqvar, start_pos, end_pos)
        if response and response.result.clinvar:
            pathogenic_variants = [
                v
                for v in response.result.clinvar
                if v.referenceAssertions
                and v.referenceAssertions[0].clinicalSignificance
                in [
                    "CLINICAL_SIGNIFICANCE_LIKELY_PATHOGENIC",
                    "CLINICAL_SIGNIFICANCE_PATHOGENIC",
                ]
            ]
            logger.debug(
                "Pathogenic variants: {}, Total variants: {}",
                len(pathogenic_variants),
                len(response.result.clinvar),
            )
            return len(pathogenic_variants), len(response.result.clinvar)
        else:
            logger.error("Failed to get variant from range. No ClinVar data.")
            raise InvalidAPIResposeError("Failed to get variant from range. No ClinVar data.")

    @staticmethod
    def _get_uniprot_domain(variant_info: VariantInfo) -> Optional[Tuple[int, int]]:
        """Check if the variant is in a UniProt domain."""
        # TODO: Implement this method
        return None

    def predict(self) -> PM1:
        """Predict PM1 criteria."""
        try:
            if self.seqvar.chrom == "MT":
                # skipped according to McCormick et al. (2020).
                self.prediction.PM1 = False
                return self.prediction

            pathogenic_count, _ = self._count_pathogenic_variants(
                self.seqvar, self.seqvar.pos - 25, self.seqvar.pos + 25
            )

            if pathogenic_count >= 4:
                self.prediction.PM1 = True
                return self.prediction

            uniprot_domain = self._get_uniprot_domain(self.variant_info)
            if not uniprot_domain:
                self.prediction.PM1 = False
                return self.prediction
            start_pos, end_pos = uniprot_domain
            pathogenic_count, _ = self._count_pathogenic_variants(self.seqvar, start_pos, end_pos)
            if pathogenic_count >= 2:
                self.prediction.PM1 = True
                return self.prediction

            self.prediction.PM1 = False
        except AutoAcmgBaseException as e:
            logger.error("Error occurred during PM1 prediction. Error: {}", e)
            self.prediction.PM1 = False

        return self.prediction
