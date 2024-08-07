"""Implementation of PP2 and BP1 criteria."""

from typing import Optional, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    PP2BP1,
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
)
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, InvalidAPIResposeError
from src.defs.seqvar import SeqVar
from src.utils import AutoACMGHelper


class AutoPP2BP1(AutoACMGHelper):
    """Class for automatic PP2 and BP1 prediction."""

    def __init__(self):
        super().__init__()
        #: Prediction result.
        self.prediction_pp2bp1: Optional[PP2BP1] = None
        #: Comment to store the prediction explanation.
        self.comment_pp2bp1: str = ""

    def _get_missense_vars(
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
            number of missense variants.

        Raises:
            AlgorithmError: If end position is less than the start position.
            InvalidAPIResposeError: If the API response is invalid or cannot be processed.
        """
        if end_pos < start_pos:
            logger.error(
                "End position is less than the start position. Positions given: {} - {}",
                start_pos,
                end_pos,
            )
            raise AlgorithmError("End position is less than the start position.")

        response = self.annonars_client.get_variant_from_range(seqvar, start_pos, end_pos)
        if response and response.clinvar:

            pathogenic_variants = [
                v
                for v in response.clinvar
                if (
                    v.records
                    and v.records[0].classifications
                    and v.records[0].classifications.germlineClassification
                    and v.records[0].classifications.germlineClassification.description
                    in ["Pathogenic"]
                    and v.records[0].variationType == "VARIATION_TYPE_SNV"
                )
            ]
            benign_variants = [
                v
                for v in response.clinvar
                if (
                    v.records
                    and v.records[0].classifications
                    and v.records[0].classifications.germlineClassification
                    and v.records[0].classifications.germlineClassification.description
                    in ["Benign"]
                    and v.records[0].variationType == "VARIATION_TYPE_SNV"
                )
            ]

            return (
                len(pathogenic_variants),
                len(benign_variants),
                len(pathogenic_variants) + len(benign_variants),
            )
        else:
            logger.error("Failed to get variant from range. No ClinVar data.")
            raise InvalidAPIResposeError("Failed to get variant from range. No ClinVar data.")

    @staticmethod
    def _is_missense(var_data: AutoACMGData) -> bool:
        """
        Check if the variant is a missense variant.

        Args:
            var_data (AutoACMGData): The variant information.

        Returns:
            bool: True if the variant is a missense variant, False otherwise.
        """
        if "missense" in var_data.consequence.cadd:
            return True
        if any("missense" in cons for cons in var_data.consequence.mehari):
            return True
        return False

    def verify_pp2bp1(self, seqvar: SeqVar, var_data: AutoACMGData) -> Tuple[Optional[PP2BP1], str]:
        """Predict PP2 and BP1 criteria."""
        self.prediction_pp2bp1 = PP2BP1()
        self.comment_pp2bp1 = ""
        if seqvar.chrom == "MT":
            self.comment_pp2bp1 = (
                "Variant is in mitochondrial DNA. PP2 and BP1 criteria are not met."
            )
            self.prediction_pp2bp1.PP2, self.prediction_pp2bp1.BP1 = False, False
        else:
            try:
                if not self._is_missense(var_data):
                    self.comment_pp2bp1 = (
                        "Variant is not a missense variant. PP2 and BP1 criteria are not met."
                    )
                    self.prediction_pp2bp1.PP2 = False
                    self.prediction_pp2bp1.BP1 = False
                    return self.prediction_pp2bp1, self.comment_pp2bp1

                start_pos, end_pos = min(var_data.cds_start, var_data.cds_end), max(
                    var_data.cds_start, var_data.cds_end
                )
                pathogenic_count, benign_count, total_count = self._get_missense_vars(
                    seqvar, start_pos, end_pos
                )
                pathogenic_ratio = pathogenic_count / total_count
                benign_ratio = benign_count / total_count
                self.comment_pp2bp1 = (
                    f"Found pathogenic missense variants: {pathogenic_count}, "
                    f"benign missense variants: {benign_count}, "
                    f"total missense variants: {total_count} "
                    f"in the range {start_pos}-{end_pos}. "
                    f"Pathogenic ratio: {pathogenic_ratio}, Benign ratio: {benign_ratio}. "
                )

                if pathogenic_ratio > var_data.thresholds.pp2bp1_pathogenic:
                    self.comment_pp2bp1 += (
                        f"Pathogenic ratio is greater than {var_data.thresholds.pp2bp1_pathogenic}. "
                        "PP2 is met. "
                    )
                    self.prediction_pp2bp1.PP2 = True
                else:
                    self.comment_pp2bp1 += (
                        f"Pathogenic ratio is less than {var_data.thresholds.pp2bp1_pathogenic}. "
                        "PP2 is not met. "
                    )
                if benign_ratio > var_data.thresholds.pp2bp1_benign:
                    self.comment_pp2bp1 += (
                        f"Benign ratio is greater than {var_data.thresholds.pp2bp1_benign}. "
                        "BP1 is met."
                    )
                    self.prediction_pp2bp1.BP1 = True
                else:
                    self.comment_pp2bp1 += (
                        f"Benign ratio is less than {var_data.thresholds.pp2bp1_benign}. "
                        "BP1 is not met."
                    )

            except AutoAcmgBaseException as e:
                logger.error("Failed to predict PP2 and BP1 criteria. Error: {}", e)
                self.comment_pp2bp1 = f"Error occurred during PP2 and BP1 prediction. Error: {e}"
                self.prediction_pp2bp1 = None
        return self.prediction_pp2bp1, self.comment_pp2bp1

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Predict PP2 and BP1 criteria for the provided sequence variant.

        Returns:
            AutoACMGCriteria: PP2 and BP1 prediction.
        """
        pred, comment = self.verify_pp2bp1(seqvar, var_data)
        if pred:
            pp2_pred = (
                AutoACMGPrediction.Met
                if pred.PP2
                else (AutoACMGPrediction.NotMet if pred.PP2 is False else AutoACMGPrediction.Failed)
            )
            bp1_pred = (
                AutoACMGPrediction.Met
                if pred.BP1
                else (AutoACMGPrediction.NotMet if pred.BP1 is False else AutoACMGPrediction.Failed)
            )
            pp2_strength = pred.PP2_strength
            bp1_strength = pred.BP1_strength
        else:
            pp2_pred = AutoACMGPrediction.Failed
            bp1_pred = AutoACMGPrediction.Failed
            pp2_strength = AutoACMGStrength.PathogenicSupporting
            bp1_strength = AutoACMGStrength.BenignSupporting
        return (
            AutoACMGCriteria(
                name="PP2",
                prediction=pp2_pred,
                strength=pp2_strength,
                summary=comment,
            ),
            AutoACMGCriteria(
                name="BP1",
                prediction=bp1_pred,
                strength=bp1_strength,
                summary=comment,
            ),
        )
