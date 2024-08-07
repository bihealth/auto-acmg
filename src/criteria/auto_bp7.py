"""Implementation of BP7 criteria."""

from typing import Optional, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    BP7,
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
)
from src.defs.exceptions import AutoAcmgBaseException, MissingDataError
from src.defs.seqvar import SeqVar
from src.utils import AutoACMGHelper


class AutoBP7(AutoACMGHelper):
    """Class for BP7 prediction."""

    def __init__(self):
        super().__init__()
        #: Prediction result.
        self.prediction_bp7: Optional[BP7] = None
        #: Comment to store the prediction explanation.
        self.comment_bp7: str = ""

    @staticmethod
    def _spliceai_impact(var_data: AutoACMGData) -> bool:
        """
        Predict splice site alterations using SpliceAI.

        If any of SpliceAI scores are greater than specific thresholds, the variant is considered a
        splice site alteration.

        Args:
            var_data: The data containing variant scores and thresholds.

        Returns:
            bool: True if the variant is a splice site alteration, False otherwise.
        """
        score_checks = {
            "spliceAI_acceptor_gain": var_data.thresholds.spliceAI_acceptor_gain,
            "spliceAI_acceptor_loss": var_data.thresholds.spliceAI_acceptor_loss,
            "spliceAI_donor_gain": var_data.thresholds.spliceAI_donor_gain,
            "spliceAI_donor_loss": var_data.thresholds.spliceAI_donor_loss,
        }
        return any(
            getattr(var_data.scores.cadd, score_type, 0) > threshold
            for score_type, threshold in score_checks.items()
        )

    @staticmethod
    def _is_conserved(var_data: AutoACMGData) -> bool:
        """
        Predict if the variant is conserved.

        Check if the variant is conserved using the phyloP100 score.

        Args:
            variant_info: The variant information.

        Returns:
            bool: True if the variant is not conserved, False otherwise.
        """
        phylop = var_data.scores.cadd.phyloP100 or var_data.scores.dbnsfp.phyloP100
        if not phylop:
            raise MissingDataError("Missing phyloP100 score.")
        if phylop >= var_data.thresholds.phyloP100:
            return True
        return False

    def verify_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> Tuple[Optional[BP7], str]:
        """Predict BP7 criterion."""
        self.prediction_bp7 = BP7()
        self.comment_bp7 = ""
        try:
            if seqvar.chrom == "MT":
                self.comment_bp7 = "Variant is in the mitochondrial genome. BP7 is not met."
                self.prediction_bp7.BP7 = False
                return self.prediction_bp7, self.comment_bp7

            if not self._is_conserved(var_data) and not self._spliceai_impact(var_data):
                self.comment_bp7 += (
                    "Variant is not conserved and not predicted to affect splicing. "
                    f"PhyloP100 score: {var_data.scores.cadd.phyloP100}. "
                    f"SpliceAI scores: {var_data.scores.cadd.spliceAI_acceptor_gain},"
                    f"{var_data.scores.cadd.spliceAI_acceptor_loss}, "
                    f"{var_data.scores.cadd.spliceAI_donor_gain}, "
                    f"{var_data.scores.cadd.spliceAI_donor_loss}. "
                    "BP7 is met."
                )
                self.prediction_bp7.BP7 = True
            else:
                self.comment_bp7 += (
                    "Variant is conserved or predicted to affect splicing. "
                    f"PhyloP100 score: {var_data.scores.cadd.phyloP100}. "
                    f"SpliceAI scores: {var_data.scores.cadd.spliceAI_acceptor_gain},"
                    f"{var_data.scores.cadd.spliceAI_acceptor_loss}, "
                    f"{var_data.scores.cadd.spliceAI_donor_gain}, "
                    f"{var_data.scores.cadd.spliceAI_donor_loss}. "
                    "BP7 is not met."
                )
                self.prediction_bp7.BP7 = False

        except AutoAcmgBaseException as e:
            logger.error("Failed to predict BP7 criterion. Error: {}", e)
            self.comment_bp7 = f"Failed to predict BP7 criterion. Error: {e}"
            self.prediction_bp7 = None

        return self.prediction_bp7, self.comment_bp7

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Predict BP7 criterion."""
        pred, comment = self.verify_bp7(seqvar, var_data)
        if pred:
            pred_bp7 = (
                AutoACMGPrediction.Met
                if pred.BP7
                else (AutoACMGPrediction.NotMet if pred.BP7 is False else AutoACMGPrediction.Failed)
            )
        else:
            pred_bp7 = AutoACMGPrediction.Failed
        return AutoACMGCriteria(
            name="BP7",
            prediction=pred_bp7,
            strength=pred.BP7_strength if pred else AutoACMGStrength.BenignSupporting,
            summary=comment,
        )
