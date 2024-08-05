"""Implementation of PM1 criteria."""

import os
from typing import Optional, Tuple

import tabix  # type: ignore
from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config, settings
from src.criteria.auto_criteria import AutoACMGCriteria
from src.defs.annonars_variant import VariantResult
from src.defs.auto_acmg import PM1, AutoACMGData, AutoACMGPrediction
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, InvalidAPIResposeError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.utils import AutoACMGHelper


class AutoPM1(AutoACMGHelper):
    """Class for automatic PM1 prediction."""

    def __init__(self):
        super().__init__()
        #: Prediction result.
        self.prediction_pm1: Optional[PM1] = None
        #: comment_pm1 to store the prediction explanation.
        self.comment_pm1: str = ""

    def _count_vars(self, seqvar: SeqVar, start_pos: int, end_pos: int) -> Tuple[int, int]:
        """
        Counts pathogenic and benign variants in the specified range.

        The method retrieves variants from the specified range and iterates through the ClinVar data
        of each variant to count the number of pathogenic and benign variants.

        Args:
            seqvar: The sequence variant being analyzed.
            start_pos: The start position of the range.
            end_pos: The end position of the range.

        Returns:
            Tuple[int, int]: The number of pathogenic and benign variants.

        Raises:
            InvalidAPIResposeError: If the API response is invalid or cannot be processed.
        """
        if end_pos < start_pos:
            logger.error("End position is less than the start position.")
            logger.debug("Positions given: {} - {}", start_pos, end_pos)
            raise AlgorithmError("End position is less than the start position.")

        response = self.annonars_client.get_variant_from_range(seqvar, start_pos, end_pos)
        if response and response.clinvar:
            pathogenic_variants = [
                v
                for v in response.clinvar
                if v.records
                and v.records[0].classifications
                and v.records[0].classifications.germlineClassification
                and v.records[0].classifications.germlineClassification.description
                in [
                    "Pathogenic",
                    "Likely pathogenic",
                ]
                and v.records[0].variationType == "VARIATION_TYPE_SNV"
            ]
            benign_variants = [
                v
                for v in response.clinvar
                if v.records
                and v.records[0].classifications
                and v.records[0].classifications.germlineClassification
                and v.records[0].classifications.germlineClassification.description
                in [
                    "Benign",
                    "Likely benign",
                ]
                and v.records[0].variationType == "VARIATION_TYPE_SNV"
            ]
            return len(pathogenic_variants), len(benign_variants)
        else:
            self.comment_pm1 += "Missing ClinVar data."
            logger.error("Failed to get variant from range. No ClinVar data.")
            raise InvalidAPIResposeError("Failed to get variant from range. No ClinVar data.")

    def _get_uniprot_domain(self, seqvar: SeqVar) -> Optional[Tuple[int, int]]:
        """Check if the variant is in a UniProt domain."""
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

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Predict PM1 criteria."""
        self.prediction = PM1()
        if seqvar.chrom == "MT":
            # skipped according to McCormick et al. (2020).
            self.comment_pm1 = "The variant is in the mitochondrial genome. PM1 is not met."
            self.prediction.PM1 = False
        else:
            try:
                pathogenic_count, benign_count = self._count_vars(
                    seqvar, seqvar.pos - 25, seqvar.pos + 25
                )
                if pathogenic_count >= var_data.thresholds.pm1_pathogenic:
                    self.comment_pm1 += "Found 8 or more pathogenic variants. PM1 is met."
                    self.prediction.PM1 = True
                    return self.prediction, self.comment_pm1
                else:
                    self.comment_pm1 += "Found less than 8 pathogenic variants."

                uniprot_domain = self._get_uniprot_domain(seqvar)
                if not uniprot_domain:
                    self.prediction.PM1 = False
                    return self.prediction, self.comment_pm1

                start_pos, end_pos = uniprot_domain
                pathogenic_count, benign_count = self._count_vars(seqvar, start_pos, end_pos)
                self.comment_pm1 += f"Found {pathogenic_count} Pathogenic variants and {benign_count} Benign variants."
                if pathogenic_count >= (end_pos - start_pos) / 4:
                    self.prediction.PM1 = True
                    return self.prediction, self.comment_pm1
                else:
                    self.prediction.PM1 = False
            except AutoAcmgBaseException as e:
                self.comment_pm1 += f"Error occurred during PM1 prediction. Error: {e}"
                logger.error("Error occurred during PM1 prediction. Error: {}", e)
                self.prediction = None

        return AutoACMGCriteria(
            name="PM1",
            summary=self.comment_pm1,
            prediction=(
                AutoACMGPrediction.Met
                if self.prediction.PM1
                else (
                    AutoACMGPrediction.NotMet
                    if self.prediction.PM1 is False
                    else AutoACMGPrediction.Failed
                )
            ),
        )
