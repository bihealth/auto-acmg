"""Implementation of the PP3 and BP4 criteria."""

from typing import Optional, Tuple

from loguru import logger

from src.criteria.auto_criteria import AutoACMGCriteria
from src.defs.auto_acmg import PP3BP4, AutoACMGData, AutoACMGPrediction
from src.defs.exceptions import AutoAcmgBaseException, MissingDataError
from src.defs.seqvar import SeqVar
from src.utils import AutoACMGHelper


class AutoPP3BP4(AutoACMGHelper):
    """Class for automatic PP3 and BP4 prediction."""

    def __init__(self):
        super().__init__()
        #: Prediction result.
        self.prediction_pp3bp4: Optional[PP3BP4] = None
        #: Comment to store the prediction explanation.
        self.comment_pp3bp4: str = ""

    def _splice_variant(self, var_data: AutoACMGData) -> bool:
        """
        Check if the variant is a splice variant.

        Args:
            var_data (AutoACMGData): The variant information.

        Returns:
            bool: True if the variant is a splice variant, False otherwise.
        """
        if "splice" in var_data.consequence.cadd:
            return True
        if any("splice" in cons for cons in var_data.consequence.mehari):
            return True
        return False

    def _is_pathogenic_score(self, var_data: AutoACMGData) -> bool:
        """
        Check if any of the pathogenic scores meet the threshold.

        Check if any of the pathogenic scores meet the threshold. If the variant is pathogenic
        based on the scores, return True.

        Args:
            variant_info (VariantResult): Variant information.

        Returns:
            bool: True if the variant is pathogenic, False otherwise.

        Raises:
            MissingDataError: If the variant information is missing.
        """
        if (
            var_data.scores.dbnsfp.metaRNN
            and var_data.scores.dbnsfp.metaRNN >= var_data.thresholds.metaRNN_pathogenic
        ):
            return True
        if (
            var_data.scores.dbnsfp.bayesDel_noAF
            and var_data.scores.dbnsfp.bayesDel_noAF >= var_data.thresholds.bayesDel_noAF_pathogenic
        ):
            return True
        return False

    def _is_benign_score(self, var_data: AutoACMGData) -> bool:
        """
        Check if any of the benign scores meet the threshold.

        Check if any of the benign scores meet the threshold. If the variant is benign
        based on the scores, return True.

        Args:
            variant_info (VariantResult): Variant information.

        Returns:
            bool: True if the variant is benign, False otherwise.

        Raises:
            MissingDataError: If the variant information is missing.
        """
        if (
            var_data.scores.dbnsfp.metaRNN
            and var_data.scores.dbnsfp.metaRNN <= var_data.thresholds.metaRNN_benign
        ):
            return True
        if (
            var_data.scores.dbnsfp.bayesDel_noAF
            and var_data.scores.dbnsfp.bayesDel_noAF <= var_data.thresholds.bayesDel_noAF_benign
        ):
            return True
        return False

    def _is_pathogenic_splicing(self, var_data: AutoACMGData) -> bool:
        """
        Check if the variant is pathogenic based on splicing scores.

        Args:
            variant_info (VariantResult): Variant information.

        Returns:
            bool: True if the variant is pathogenic, False otherwise.

        Raises:
            MissingDataError: If the variant information is missing.
        """
        ada = var_data.scores.dbscsnv.ada or var_data.scores.cadd.ada
        rf = var_data.scores.dbscsnv.rf or var_data.scores.cadd.rf
        if not ada and not rf:
            raise MissingDataError("Missing Ada and RF scores.")
        if ada >= var_data.thresholds.ada or rf >= var_data.thresholds.rf:
            return True
        return False

    def _is_benign_splicing(self, var_data: AutoACMGData) -> bool:
        """
        Check if the variant is benign based on splicing scores.

        Args:
            variant_info (VariantResult): Variant information.

        Returns:
            bool: True if the variant is benign, False otherwise.

        Raises:
            MissingDataError: If the variant information is missing.
        """
        ada = var_data.scores.dbscsnv.ada or var_data.scores.cadd.ada
        rf = var_data.scores.dbscsnv.rf or var_data.scores.cadd.rf
        if not ada and not rf:
            raise MissingDataError("Missing Ada and RF scores.")
        if ada <= var_data.thresholds.ada and rf <= var_data.thresholds.rf:
            return True
        return False

    def predict_pp3bp4(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Predict PP3 and BP4 criteria."""
        self.prediction_pp3bp4 = PP3BP4()

        if seqvar.chrom == "MT":
            self.comment_pp3bp4 = (
                "Variant is in mitochondrial DNA. PP3 and BP4 criteria are not met."
            )
            self.prediction_pp3bp4.PP3, self.prediction_pp3bp4.BP4 = False, False
        else:
            try:
                if self._splice_variant(var_data):
                    self.comment_pp3bp4 += "Variant is a splice variant."
                    self.prediction_pp3bp4.PP3 = self._is_pathogenic_splicing(var_data)
                    self.prediction_pp3bp4.BP4 = self._is_benign_splicing(var_data)
                else:
                    self.comment_pp3bp4 += "Variant is not a splice variant."
                    self.prediction_pp3bp4.PP3 = self._is_pathogenic_score(var_data)
                    self.prediction_pp3bp4.BP4 = self._is_benign_score(var_data)

            except AutoAcmgBaseException as e:
                self.comment_pp3bp4 += f"An error occurred during prediction. Error: {e}"
                logger.error("Failed to predict PP3 and BP4 criteria. Error: {}", e)
                self.prediction_pp3bp4 = None

        return (
            AutoACMGCriteria(
                name="PP3",
                summary=self.comment_pp3bp4,
                prediction=(
                    AutoACMGPrediction.Met
                    if self.prediction_pp3bp4.PP3
                    else (
                        AutoACMGPrediction.NotMet
                        if self.prediction_pp3bp4.PP3 is False
                        else AutoACMGPrediction.Failed
                    )
                ),
            ),
            AutoACMGCriteria(
                name="BP4",
                summary=self.comment_pp3bp4,
                prediction=(
                    AutoACMGPrediction.Met
                    if self.prediction_pp3bp4.BP4
                    else (
                        AutoACMGPrediction.NotMet
                        if self.prediction_pp3bp4.BP4 is False
                        else AutoACMGPrediction.Failed
                    )
                ),
            ),
        )
