"""Implementation of PM4 and BP3 rules for sequence variants."""

import os
from typing import Optional, Tuple

import tabix  # type: ignore
from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config, settings
from src.defs.annonars_variant import AnnonarsVariantResponse, VariantResult
from src.defs.auto_acmg import (
    PM4BP3,
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
)
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.utils import AutoACMGHelper


class AutoPM4BP3(AutoACMGHelper):
    """Predicts PM4 and BP3 criteria for sequence variants."""

    def __init__(self):
        super().__init__()
        #: Prediction result.
        self.prediction_pm4bp3: Optional[PM4BP3] = None
        #: Comment to store the prediction explanation.
        self.comment_pm4bp3: str = ""

    @staticmethod
    def _in_repeat_region(seqvar: SeqVar) -> bool:
        """Check if the variant is in a repeat region.

        Check if the variant is in a repeat region using the RepeatMasker track.

        Args:
            seqvar (SeqVar): Sequence variant.

        Returns:
            bool: True if the variant is in a repeat region, False otherwise.
        """
        try:
            # Find path to the lib file
            if seqvar.genome_release == GenomeRelease.GRCh37:
                path = os.path.join(settings.PATH_TO_ROOT, "lib", "rmsk", "grch37", "rmsk.bed.gz")
            else:
                path = os.path.join(settings.PATH_TO_ROOT, "lib", "rmsk", "grch38", "rmsk.bed.gz")
            tb = tabix.open(path)
            records = tb.query(f"chr{seqvar.chrom}", seqvar.pos - 1, seqvar.pos)
            # Check if iterator is not empty
            if any(True for _ in records):
                return True
            else:
                return False
        except tabix.TabixError as e:
            logger.error("Failed to check if the variant is in a repeat region. Error: {}", e)
            raise AlgorithmError("Failed to check if the variant is in a repeat region.") from e

    @staticmethod
    def _is_stop_loss(var_data: AutoACMGData) -> bool:
        """Check if the variant is a stop-loss variant.

        Args:
            var_data (AutoACMGData): The variant information.

        Returns:
            bool: True if the variant is a stop-loss variant, False otherwise.
        """
        if "stop_loss" in var_data.consequence.cadd:
            return True
        if any("stop_loss" in cons for cons in var_data.consequence.mehari):
            return True
        return False

    @staticmethod
    def is_inframe_delins(var_data: AutoACMGData) -> bool:
        """Check if the variant is an in-frame deletion/insertion.

        Args:
            var_data (AutoACMGData): The variant information.

        Returns:
            bool: True if the variant is an in-frame deletion/insertion, False otherwise.
        """
        if "inframe_deletion" in var_data.consequence.cadd:
            return True
        if "inframe_insertion" in var_data.consequence.cadd:
            return True
        if any("inframe_deletion" in cons for cons in var_data.consequence.mehari):
            return True
        if any("inframe_insertion" in cons for cons in var_data.consequence.mehari):
            return True
        return False

    def verify_pm4bp3(self, seqvar: SeqVar, var_data: AutoACMGData) -> Tuple[Optional[PM4BP3], str]:
        """Predicts PM4 and BP3 criteria for the provided sequence variant.

        Implementation of the rule:
        - If the variant is a stop-loss variant, PM4 is True and BP3 is False.
        - If the variant is an in-frame deletion/insertion:
        - If the variant is not in a repeat region, PM4 is True and BP3 is False.
        - If the variant is in a repeat region, PM4 is False and BP3 is True.
        - Otherwise, PM4 and BP3 are False.

        Note:
            Rules:
            PM4: Protein length changes due to in-frame deletions/insertions in a non-repeat region
            or stop-loss variants.
            BP3: In-frame deletions/insertions in a repetitive region without a known function.

        Returns:
            PM4BP3: PM4 and BP3 prediction.
        """
        self.prediction_pm4bp3 = PM4BP3()
        self.comment_pm4bp3 = ""
        try:
            # Stop-loss variants are considered as PM4
            if self._is_stop_loss(var_data):
                self.comment_pm4bp3 = "Variant consequence is stop-loss. PM4 is met."
                self.prediction_pm4bp3.PM4 = True
            # In-frame deletions/insertions
            elif self.is_inframe_delins(var_data):
                self.comment_pm4bp3 += f"Variant consequence is in-frame deletion/insertion.\n"
                if not self._in_repeat_region(seqvar):
                    self.comment_pm4bp3 += (
                        "Variant is not in a repeat region or a conserved domain. PM4 is met."
                    )
                    self.prediction_pm4bp3.PM4 = True
                else:
                    self.comment_pm4bp3 += (
                        "Variant is in a repeat region or not in a conserved domain. BP3 is met."
                    )
                    self.prediction_pm4bp3.BP3 = True
            else:
                self.comment_pm4bp3 = (
                    "Variant consequence is not indel or stop-loss. PM4 and BP3 are not met."
                )
                self._in_repeat_region(seqvar)
                self.prediction_pm4bp3.PM4 = False
                self.prediction_pm4bp3.BP3 = False

        except AutoAcmgBaseException as e:
            logger.error("Failed to predict PM4 and BP3 criteria. Error: {}", e)
            self.comment_pm4bp3 = f"An error occured while predicting PM4 and BP3 criteria: {e}"
            self.prediction_pm4bp3 = None

        return self.prediction_pm4bp3, self.comment_pm4bp3

    def predict_pm4bp3(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Predict PM4 and BP3 criteria for the provided sequence variant.

        Returns:
            AutoACMGCriteria: PM4 and BP3 prediction.
        """
        pred, comment = self.verify_pm4bp3(seqvar, var_data)
        if pred:
            pm4_pred = (
                AutoACMGPrediction.Met
                if pred.PM4
                else (AutoACMGPrediction.NotMet if pred.PM4 is False else AutoACMGPrediction.Failed)
            )
            bp3_pred = (
                AutoACMGPrediction.Met
                if pred.BP3
                else (AutoACMGPrediction.NotMet if pred.BP3 is False else AutoACMGPrediction.Failed)
            )
        else:
            pm4_pred = AutoACMGPrediction.Failed
            bp3_pred = AutoACMGPrediction.Failed
        return (
            AutoACMGCriteria(
                name="PM4",
                prediction=pm4_pred,
                strength=pred.PM4_strength if pred else AutoACMGStrength.PathogenicModerate,
                summary=comment,
            ),
            AutoACMGCriteria(
                name="BP3",
                prediction=bp3_pred,
                strength=pred.BP3_strength if pred else AutoACMGStrength.BenignSupporting,
                summary=comment,
            ),
        )
