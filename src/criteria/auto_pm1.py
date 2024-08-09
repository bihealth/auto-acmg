"""Implementation of PM1 criteria."""

import os
from typing import Optional, Tuple

import tabix
from loguru import logger

from src.core.config import settings
from src.defs.auto_acmg import (
    PM1,
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
)
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
            AlgorithmError: If end position is less than the start position.
            InvalidAPIResposeError: If the API response is invalid or cannot be processed.
        """
        if end_pos < start_pos:
            logger.error("End position is less than the start position.")
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
            logger.error("Failed to get variant from range. No ClinVar data.")
            raise InvalidAPIResposeError("Failed to get variant from range. No ClinVar data.")

    @staticmethod
    def _get_uniprot_domain(seqvar: SeqVar) -> Optional[Tuple[int, int]]:
        """
        Retrieve the UniProt domain for the variant and return the start and end positions if found
        or None otherwise.

        Args:
            seqvar: The sequence variant being analyzed.

        Returns:
            Optional[Tuple[int, int]]: The start and end positions of the UniProt domain if found,
            None otherwise.
        """
        # Find path to the lib file
        if seqvar.genome_release == GenomeRelease.GRCh37:
            path = os.path.join(settings.PATH_TO_ROOT, "lib", "uniprot", "grch37", "uniprot.bed.gz")
        else:
            path = os.path.join(settings.PATH_TO_ROOT, "lib", "uniprot", "grch38", "uniprot.bed.gz")
        tb = tabix.open(path)
        records = tb.query(f"chr{seqvar.chrom}", seqvar.pos - 1, seqvar.pos)
        # Return the first record
        for record in records:
            return int(record[1]), int(record[2])
        return None

    def verify_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> Tuple[Optional[PM1], str]:
        """Predict PM1 criteria."""
        self.prediction_pm1 = PM1()
        self.comment_pm1 = ""
        if seqvar.chrom == "MT":
            # skipped according to McCormick et al. (2020).
            self.comment_pm1 = "The variant is in the mitochondrial genome. PM1 is not met."
            self.prediction_pm1.PM1 = False
        else:
            try:
                pathogenic_count, benign_count = self._count_vars(
                    seqvar, seqvar.pos - 25, seqvar.pos + 25
                )
                if pathogenic_count >= var_data.thresholds.pm1_pathogenic:
                    self.comment_pm1 += (
                        f"Found {var_data.thresholds.pm1_pathogenic} or more pathogenic variants. "
                        "PM1 is met."
                    )
                    self.prediction_pm1.PM1 = True
                    return self.prediction_pm1, self.comment_pm1
                else:
                    self.comment_pm1 += (
                        f"Found less than {var_data.thresholds.pm1_pathogenic} pathogenic "
                        "variants. "
                    )

                uniprot_domain = self._get_uniprot_domain(seqvar)
                if not uniprot_domain:
                    self.comment_pm1 += "Variant is not in a UniProt domain. PM1 is not met."
                    self.prediction_pm1.PM1 = False
                    return self.prediction_pm1, self.comment_pm1

                start_pos, end_pos = uniprot_domain
                pathogenic_count, benign_count = self._count_vars(seqvar, start_pos, end_pos)
                self.comment_pm1 += (
                    f"Found {pathogenic_count} Pathogenic variants "
                    f"and {benign_count} Benign variants "
                    f"in Uniprot Domain: {start_pos}-{end_pos}. "
                )
                if pathogenic_count >= (end_pos - start_pos) / 4:
                    self.comment_pm1 += "PM1 is met."
                    self.prediction_pm1.PM1 = True
                    return self.prediction_pm1, self.comment_pm1
                else:
                    self.comment_pm1 += "PM1 is not met."
                    self.prediction_pm1.PM1 = False
            except AutoAcmgBaseException as e:
                logger.error("Error occurred during PM1 prediction. Error: {}", e)
                self.comment_pm1 = f"Error occurred during PM1 prediction. Error: {e}"
                self.prediction_pm1 = None
        return self.prediction_pm1, self.comment_pm1

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Predict PM1 criteria."""
        pred, comment = self.verify_pm1(seqvar, var_data)
        if pred:
            pm1_pred = (
                AutoACMGPrediction.Met
                if pred.PM1
                else (AutoACMGPrediction.NotMet if pred.PM1 is False else AutoACMGPrediction.Failed)
            )
            pm1_strength = pred.PM1_strength
        else:
            pm1_pred = AutoACMGPrediction.Failed
            pm1_strength = AutoACMGStrength.PathogenicModerate
        return AutoACMGCriteria(
            name="PM1",
            prediction=pm1_pred,
            strength=pm1_strength,
            summary=comment,
        )
