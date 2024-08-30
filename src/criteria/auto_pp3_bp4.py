"""Implementation of the PP3 and BP4 criteria."""

from typing import Optional, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    PP3BP4,
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
)
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

    def _is_splice_variant(self, var_data: AutoACMGData) -> bool:
        """
        Check if the variant's consequence is a splice related.

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

    def _is_inframe_indel(self, var_data: AutoACMGData) -> bool:
        """
        Check if the variant's consequence is an inframe indel.

        Args:
            var_data (AutoACMGData): The variant information.

        Returns:
            bool: True if the variant is an inframe indel, False otherwise.
        """
        if "inframe" in var_data.consequence.cadd:
            return True
        if any("inframe" in cons for cons in var_data.consequence.mehari):
            return True
        return False

    def _is_missense_variant(self, var_data: AutoACMGData) -> bool:
        """
        Check if the variant's consequence is a missense variant.

        Args:
            var_data (AutoACMGData): The variant information.

        Returns:
            bool: True if the variant is a missense variant, False otherwise.
        """
        if "missense" in var_data.consequence.cadd:
            return True
        if "missense_variant" in var_data.consequence.mehari:
            return True
        return False

    def _is_synonymous_variant(self, var_data: AutoACMGData) -> bool:
        """
        Check if the variant's consequence is a synonymous variant.

        Args:
            var_data (AutoACMGData): The variant information.

        Returns:
            bool: True if the variant is a synonymous variant, False otherwise.
        """
        if "synonymous" in var_data.consequence.cadd:
            return True
        if "synonymous_variant" in var_data.consequence.mehari:
            return True
        return False

    def _is_intron_variant(self, var_data: AutoACMGData) -> bool:
        """
        Check if the variant's consequence is an intron variant.

        Args:
            var_data (AutoACMGData): The variant information.

        Returns:
            bool: True if the variant is an intron variant, False otherwise.
        """
        if "intron" in var_data.consequence.cadd:
            return True
        if any("intron" in cons for cons in var_data.consequence.mehari):
            return True
        return False

    def _is_utr_variant(self, var_data: AutoACMGData) -> bool:
        """
        Check if the variant's consequence is an UTR variant.

        Args:
            var_data (AutoACMGData): The variant information.

        Returns:
            bool: True if the variant is an UTR variant, False otherwise.
        """
        if (x in var_data.consequence.cadd for x in ["UTR", "utr"]):
            return True
        if any("utr" in cons for cons in var_data.consequence.mehari) or any(
            "UTR" in cons for cons in var_data.consequence.mehari
        ):
            return True
        return False

    def _is_pathogenic_score(
        self, var_data: AutoACMGData, *score_threshold_pairs: Tuple[str, float]
    ) -> bool:
        """
        Check if any of the specified scores meet their corresponding threshold.

        Args:
            var_data (AutoACMGData): Variant data containing scores and thresholds.
            score_threshold_pairs (Tuple[str, float]): Pairs of score attributes and their corresponding pathogenic thresholds.

        Returns:
            bool: True if any of the specified scores meet their corresponding threshold, False otherwise.
        """
        for score_attr, threshold in score_threshold_pairs:
            score_value = getattr(var_data.scores.dbnsfp, score_attr, None)
            if score_value is not None and score_value >= threshold:
                return True
        return False

    def _is_benign_score(
        self, var_data: AutoACMGData, *score_threshold_pairs: Tuple[str, float]
    ) -> bool:
        """
        Check if any of the specified scores meet their corresponding threshold.

        Args:
            var_data (AutoACMGData): Variant data containing scores and thresholds.
            score_threshold_pairs (Tuple[str, float]): Pairs of score attributes and their corresponding benign thresholds.

        Returns:
            bool: True if any of the specified scores meet their corresponding threshold, False otherwise.
        """
        for score_attr, threshold in score_threshold_pairs:
            score_value = getattr(var_data.scores.dbnsfp, score_attr, None)
            if score_value is not None and score_value <= threshold:
                return True
        return False

    def _affect_spliceAI(self, var_data: AutoACMGData) -> bool:
        """
        Predict splice site alterations using SpliceAI.

        If any of SpliceAI scores are greater than specific thresholds, the variant is considered a
        splice site alteration. The thresholds are defined in the variant data thresholds.

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
            (getattr(var_data.scores.cadd, score_name) or 0) > threshold
            for score_name, threshold in score_checks.items()
        )

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
        if ada:
            if ada > var_data.thresholds.ada:
                return True
        if rf:
            if rf > var_data.thresholds.rf:
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
        if ada:
            if ada < var_data.thresholds.ada:
                return True
        if rf:
            if rf < var_data.thresholds.rf:
                return True
        return False

    def verify_pp3bp4(self, seqvar: SeqVar, var_data: AutoACMGData) -> Tuple[Optional[PP3BP4], str]:
        """Predict PP3 and BP4 criteria."""
        self.prediction_pp3bp4 = PP3BP4()
        self.comment_pp3bp4 = ""
        if seqvar.chrom == "MT":
            self.comment_pp3bp4 = (
                "Variant is in mitochondrial DNA. PP3 and BP4 criteria are not met."
            )
            self.prediction_pp3bp4.PP3, self.prediction_pp3bp4.BP4 = False, False
        else:
            try:
                # if self._is_splice_variant(var_data):
                #     self.comment_pp3bp4 = "Variant is a splice variant."
                #     self.prediction_pp3bp4.PP3 = self._is_pathogenic_splicing(var_data)
                #     self.prediction_pp3bp4.BP4 = self._is_benign_splicing(var_data)
                #     self.comment_pp3bp4 += (
                #         f"Ada score: {var_data.scores.dbscsnv.ada}, "
                #         f"Ada threshold: {var_data.thresholds.ada}. "
                #         f"RF score: {var_data.scores.dbscsnv.rf}, "
                #         f"RF threshold: {var_data.thresholds.rf}. "
                #     )
                # else:
                #     self.comment_pp3bp4 = "Variant is not a splice variant."
                if (score := var_data.thresholds.pp3bp4_strategy) == "default":
                    self.prediction_pp3bp4.PP3 = self._is_pathogenic_score(
                        var_data,
                        ("metaRNN", var_data.thresholds.metaRNN_pathogenic),
                        ("bayesDel_noAF", var_data.thresholds.bayesDel_noAF_pathogenic),
                    )
                    self.prediction_pp3bp4.BP4 = self._is_benign_score(
                        var_data,
                        ("metaRNN", var_data.thresholds.metaRNN_benign),
                        ("bayesDel_noAF", var_data.thresholds.bayesDel_noAF_benign),
                    )
                    self.comment_pp3bp4 += (
                        f"MetaRNN score: {var_data.scores.dbnsfp.metaRNN}, "
                        f"MetaRNN threshold: {var_data.thresholds.metaRNN_pathogenic}. "
                        f"BayesDel_noAF score: {var_data.scores.dbnsfp.bayesDel_noAF}, "
                        f"BayesDel_noAF threshold: {var_data.thresholds.bayesDel_noAF_pathogenic}. "
                    )
                else:
                    self.prediction_pp3bp4.PP3 = self._is_pathogenic_score(
                        var_data,
                        (score, getattr(var_data.thresholds, f"{score}_pathogenic")),
                    )
                    self.prediction_pp3bp4.BP4 = self._is_benign_score(
                        var_data,
                        (score, getattr(var_data.thresholds, f"{score}_benign")),
                    )

                self.prediction_pp3bp4.PP3 = self.prediction_pp3bp4.PP3 or self._is_pathogenic_splicing(var_data)
                self.prediction_pp3bp4.BP4 = self.prediction_pp3bp4.BP4 or self._is_benign_splicing(var_data)
                self.comment_pp3bp4 += (
                    f"Ada score: {var_data.scores.dbscsnv.ada}, "
                    f"Ada threshold: {var_data.thresholds.ada}. "
                    f"RF score: {var_data.scores.dbscsnv.rf}, "
                    f"RF threshold: {var_data.thresholds.rf}. "
                )

            except AutoAcmgBaseException as e:
                logger.error("Failed to predict PP3 and BP4 criteria. Error: {}", e)
                self.comment_pp3bp4 = f"An error occurred during prediction. Error: {e}"
                self.prediction_pp3bp4 = None
        return self.prediction_pp3bp4, self.comment_pp3bp4

    def predict_pp3bp4(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Predict PP3 and BP4 criteria."""
        logger.info("Predict PP3 and BP4")
        pred, comment = self.verify_pp3bp4(seqvar, var_data)
        if pred:
            pp3_pred = (
                AutoACMGPrediction.Met
                if pred.PP3
                else (AutoACMGPrediction.NotMet if pred.PP3 is False else AutoACMGPrediction.Failed)
            )
            bp4_pred = (
                AutoACMGPrediction.Met
                if pred.BP4
                else (AutoACMGPrediction.NotMet if pred.BP4 is False else AutoACMGPrediction.Failed)
            )
            pp3_strength = pred.PP3_strength
            bp4_strength = pred.BP4_strength
        else:
            pp3_pred = AutoACMGPrediction.Failed
            bp4_pred = AutoACMGPrediction.Failed
            pp3_strength = AutoACMGStrength.PathogenicSupporting
            bp4_strength = AutoACMGStrength.BenignSupporting
        return (
            AutoACMGCriteria(
                name="PP3",
                prediction=pp3_pred,
                strength=pp3_strength,
                summary=comment,
            ),
            AutoACMGCriteria(
                name="BP4",
                prediction=bp4_pred,
                strength=bp4_strength,
                summary=comment,
            ),
        )
