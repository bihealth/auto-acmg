"""PVS1 criteria for Structural Variants (StrucVar)."""

from loguru import logger

from src.defs.auto_pvs1 import PVS1Prediction, PVS1PredictionStrucVarPath
from src.defs.strucvar import StrucVar, StrucVarType


class StrucVarPVS1Helper:
    """Helper methods for PVS1 criteria for structural variants."""

    @staticmethod
    def _full_gene_deletion() -> bool:
        """Check if the variant is a full gene deletion."""
        return False

    @staticmethod
    def _deletion_disrupts_rf() -> bool:
        """Check if the single or multiple exon deletion disrupts the reading frame."""
        return False

    @staticmethod
    def _undergo_nmd() -> bool:
        """Check if the deletion is expected to undergo NMD."""
        return False

    @staticmethod
    def _in_biologically_relevant_transcript() -> bool:
        """Check if the deletion is in a biologically relevant transcript."""
        return False

    @staticmethod
    def _critical4protein_function() -> bool:
        """Check if the deletion is critical for protein function."""
        return False

    @staticmethod
    def _lof_is_frequent_in_population() -> bool:
        """Check if loss-of-function is frequent in the population."""
        return False

    @staticmethod
    def _lof_removes_more_then_10_percent_of_protein() -> bool:
        """Check if the loss-of-function removes more than 10% of the protein."""
        return False

    @staticmethod
    def _proven_in_tandem() -> bool:
        """Check if the duplication is proven in tandem."""
        return False

    @staticmethod
    def _presumed_in_tandem() -> bool:
        """Check if the duplication is presumed in tandem."""
        return False

    @staticmethod
    def _duplication_disrupts_rf() -> bool:
        """Check if the duplication disrupts the reading frame."""
        return False


class StrucVarPVS1(StrucVarPVS1Helper):
    """PVS1 prediction for structural variants."""

    def __init__(self, variant: StrucVar):
        self.variant = variant
        self.prediction: PVS1Prediction = PVS1Prediction.NotSet
        self.prediction_path: PVS1PredictionStrucVarPath = PVS1PredictionStrucVarPath.NotSet

    def initialize(self):
        """Initialize the PVS1 prediction."""
        pass

    def verify_PVS1(self):
        """Verify PVS1 prediction."""
        if self.variant.sv_type == StrucVarType.DEL:
            if self._full_gene_deletion():
                self.prediction = PVS1Prediction.PVS1
                self.prediction_path = PVS1PredictionStrucVarPath.DEL1
            elif self._deletion_disrupts_rf() and self._undergo_nmd():
                if self._in_biologically_relevant_transcript():
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL2
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL3
            elif self._deletion_disrupts_rf() and not self._undergo_nmd():
                if self._critical4protein_function():
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL4
                else:
                    if (
                        self._lof_is_frequent_in_population()
                        or not self._in_biologically_relevant_transcript()
                    ):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionStrucVarPath.DEL5_1
                    else:
                        if self._lof_removes_more_then_10_percent_of_protein():
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL6_1
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL7_1
            else:
                if self._critical4protein_function():
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL8
                else:
                    if (
                        self._lof_is_frequent_in_population()
                        or not self._in_biologically_relevant_transcript()
                    ):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionStrucVarPath.DEL5_2
                    else:
                        if self._lof_removes_more_then_10_percent_of_protein():
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL6_2
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL7_2

        elif self.variant.sv_type == StrucVarType.DUP:
            if self._proven_in_tandem():
                if self._duplication_disrupts_rf() and self._undergo_nmd():
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DUP1
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DUP2_1
            elif self._presumed_in_tandem():
                if self._duplication_disrupts_rf() and self._undergo_nmd():
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionStrucVarPath.DUP3
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DUP2_2

            else:
                self.prediction = PVS1Prediction.NotPVS1
                self.prediction_path = PVS1PredictionStrucVarPath.DUP4

        else:
            self.prediction = PVS1Prediction.NotSet
            self.prediction_path = PVS1PredictionStrucVarPath.NotSet
            logger.error("Unsupported structural variant type: {}", self.variant.sv_type)

    def get_prediction(self) -> tuple[PVS1Prediction, PVS1PredictionStrucVarPath]:
        """Get the PVS1 prediction."""
        return self.prediction, self.prediction_path
