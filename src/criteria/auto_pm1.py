"""Implementation of PM1 criteria."""

import os
from typing import Optional, Tuple

import tabix  # type: ignore
from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config, settings
from src.defs.annonars_variant import VariantResult
from src.defs.auto_acmg import PM1
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, InvalidAPIResposeError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


class AutoPM1:
    """Class for automatic PM1 prediction."""

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
        self.prediction: Optional[PM1] = None
        #: Comment to store the prediction explanation.
        self.comment: str = ""

    def _count_pathogenic_vars(
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
        if end_pos < start_pos:
            logger.error("End position is less than the start position.")
            logger.debug("Positions given: {} - {}", start_pos, end_pos)
            raise AlgorithmError("End position is less than the start position.")

        response = self.annonars_client.get_variant_from_range(seqvar, start_pos, end_pos)
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
            return len(pathogenic_variants), len(response.result.clinvar)
        else:
            logger.error("Failed to get variant from range. No ClinVar data.")
            raise InvalidAPIResposeError("Failed to get variant from range. No ClinVar data.")

    def _get_uniprot_domain(self, seqvar: SeqVar) -> Optional[Tuple[int, int]]:
        """Check if the variant is in a UniProt domain."""
        self.comment += "Check if the variant is in a UniProt domain.\n"
        try:
            # Find path to the lib file
            if seqvar.genome_release == GenomeRelease.GRCh37:
                path = os.path.join(
                    settings.PATH_TO_ROOT, "lib", "uniprot", "grch37", "uniprot.bed.gz"
                )
            else:
                path = os.path.join(
                    settings.PATH_TO_ROOT, "lib", "uniprot", "grch38", "uniprot.bed.gz"
                )
            tb = tabix.open(path)
            records = tb.query(f"chr{seqvar.chrom}", seqvar.pos - 1, seqvar.pos)
            # Return the first record
            for record in records:
                return int(record[1]), int(record[2])
            return None

        except tabix.TabixError as e:
            logger.error("Failed to check if the variant is in a UniProt domain. Error: {}", e)
            raise AlgorithmError("Failed to check if the variant is in a UniProt domain.") from e

    def predict(self) -> Tuple[Optional[PM1], str]:
        """Predict PM1 criteria."""
        self.prediction = PM1()
        try:
            if self.seqvar.chrom == "MT":
                # skipped according to McCormick et al. (2020).
                self.comment = "The variant is in the mitochondrial genome. PM1 is not met."
                self.prediction.PM1 = False
                return self.prediction, self.comment

            self.comment = (
                "Counting pathogenic variants in the range of 50bp."
                f"The range is {self.seqvar.pos - 25} - {self.seqvar.pos + 25}. => \n"
            )
            logger.debug(
                "Counting pathogenic variants in the range of 50bp."
                f"The range is {self.seqvar.pos - 25} - {self.seqvar.pos + 25}."
            )
            pathogenic_count, _ = self._count_pathogenic_vars(
                self.seqvar, self.seqvar.pos - 25, self.seqvar.pos + 25
            )

            self.comment += f"Found {pathogenic_count} Pathogenic variants. => \n"
            logger.debug("Found {} Pathogenic variants.", pathogenic_count)
            if pathogenic_count >= 4:
                self.comment += "Found 4 or more pathogenic variants. PM1 is met."
                logger.debug("Found 4 or more pathogenic variants. PM1 is met.")
                self.prediction.PM1 = True
                return self.prediction, self.comment
            else:
                self.comment += "Found less than 4 pathogenic variants."
                logger.debug("Found less than 4 pathogenic variants.")

            self.comment += "Checking if the variant is in a UniProt domain. => \n"
            logger.debug("Checking if the variant is in a UniProt domain.")
            uniprot_domain = self._get_uniprot_domain(self.seqvar)
            if not uniprot_domain:
                self.comment += "The variant is not in a UniProt domain."
                logger.debug("The variant is not in a UniProt domain.")
                self.prediction.PM1 = False
                return self.prediction, self.comment

            start_pos, end_pos = uniprot_domain
            self.comment += (
                "Counting pathogenic variants in the UniProt domain. "
                f"The range is {start_pos} - {end_pos}. => \n"
            )
            logger.debug(
                "Counting pathogenic variants in the UniProt domain. "
                f"The range is {start_pos} - {end_pos}."
            )
            pathogenic_count, _ = self._count_pathogenic_vars(self.seqvar, start_pos, end_pos)
            if pathogenic_count >= 2:
                self.comment += (
                    "Found 2 or more pathogenic variants in the UniProt domain. PM1 is met."
                )
                logger.debug(
                    "Found 2 or more pathogenic variants in the UniProt domain. PM1 is met."
                )
                self.prediction.PM1 = True
                return self.prediction, self.comment
            else:
                self.comment += (
                    "Found less than 2 pathogenic variants in the UniProt domain. "
                    "PM1 is not met."
                )
                logger.debug(
                    "Found less than 2 pathogenic variants in the UniProt domain. PM1 is not met."
                )
                self.prediction.PM1 = False
        except AutoAcmgBaseException as e:
            self.comment += f"Error occurred during PM1 prediction. Error: {e}"
            logger.error("Error occurred during PM1 prediction. Error: {}", e)
            self.prediction = None

        return self.prediction, self.comment
