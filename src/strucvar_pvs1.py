"""PVS1 criteria for Structural Variants (StrucVar)."""

from src.defs.autopvs1 import PVS1Prediction, PVS1PredictionStrucVarPath
from src.defs.strucvar import StrucVar, StrucVarType


class StrucVarPVS1:
    """PVS1 prediction for structural variants."""

    def __init__(self, variant: StrucVar):
        self.variant = variant
        self.prediction: PVS1Prediction = PVS1Prediction.NotPVS1
        self.prediction_path: PVS1PredictionStrucVarPath = PVS1PredictionStrucVarPath.NotSet

    def initialize(self):
        """Initialize the PVS1 prediction."""
        if self.variant.sv_type == StrucVarType.DEL:
            self.prediction_path = PVS1PredictionStrucVarPath.DEL1
        elif self.variant.sv_type == StrucVarType.DUP:
            self.prediction_path = PVS1PredictionStrucVarPath.DUP1

    def verify_PVS1(self):
        """Verify PVS1 prediction."""
        if self.variant.sv_type == StrucVarType.DEL:
            self.prediction = PVS1Prediction.PVS1
        elif self.variant.sv_type == StrucVarType.DUP:
            self.prediction = PVS1Prediction.PVS1

    def get_prediction(self) -> tuple[PVS1Prediction, PVS1PredictionStrucVarPath]:
        """Get the PVS1 prediction."""
        return self.prediction, self.prediction_path

    def __repr__(self):
        return f"{self.variant.user_repr} - {self.prediction.name}"
