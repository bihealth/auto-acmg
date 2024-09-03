"""PVS1 criteria for Structural Variants (StrucVar)."""

from typing import Dict, Tuple

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGPrediction, AutoACMGStrength
from src.defs.auto_pvs1 import PVS1Prediction, PVS1PredictionPathMapping, PVS1PredictionStrucVarPath
from src.utils import AutoACMGHelper


class StrucVarHelper(AutoACMGHelper):
    """Helper methods for PVS1 criteria for Structural Variants (StrucVar)."""

    def __init__(self):
        super().__init__()
        #: Comment to store the prediction explanation.
        self.comment_pvs1: str = ""

    @staticmethod
    def full_gene_del() -> bool:
        """Check if the variant is a full gene deletion."""
        return False

    @staticmethod
    def del_disrupt_rf() -> bool:
        """Check if the single or multiple exon deletion disrupts the reading frame."""
        return False

    @staticmethod
    def dup_disrupt_rf() -> bool:
        """Check if the duplication disrupts the reading frame."""
        return False

    def undergo_nmd(self) -> bool:
        """Check if the variant undergoes NMD."""
        return True

    @staticmethod
    def in_bio_relevant_tsx() -> bool:
        """Check if the deletion is in a biologically relevant transcript."""
        return False

    @staticmethod
    def crit4prot_func() -> bool:
        """Check if the deletion is critical for protein function."""
        return False

    @staticmethod
    def lof_freq_in_pop() -> bool:
        """Check if loss-of-function is frequent in the population."""
        return False

    @staticmethod
    def lof_rm_gt_10pct_of_prot() -> bool:
        """Check if the loss-of-function removes more than 10% of the protein."""
        return False

    @staticmethod
    def proven_in_tandem() -> bool:
        """Check if the duplication is proven in tandem."""
        return False

    @staticmethod
    def presumed_in_tandem() -> bool:
        """Check if the duplication is presumed in tandem."""
        return False


class AutoPVS1(StrucVarHelper):
    """Handles the PVS1 criteria assesment for structural variants."""

    def __init__(self):
        super().__init__()
        self.prediction: PVS1Prediction = PVS1Prediction.NotPVS1
        self.prediction_path: PVS1PredictionStrucVarPath = PVS1PredictionStrucVarPath.NotSet

    def verify_pvs1(self) -> Tuple[PVS1Prediction, PVS1PredictionStrucVarPath, str]:
        """Make the PVS1 prediction.

        The prediction is based on the PVS1 criteria for structural variants.
        """
        return self.prediction, self.prediction_path, self.comment_pvs1

    def predict_pvs1(self) -> AutoACMGCriteria:
        """Predict the PVS1 criteria for structural variants."""
        pred, path, comment = self.verify_pvs1()
        evidence_strength_mapping: Dict[PVS1Prediction, AutoACMGStrength] = {
            PVS1Prediction.PVS1: AutoACMGStrength.PathogenicVeryStrong,
            PVS1Prediction.PVS1_Strong: AutoACMGStrength.PathogenicStrong,
            PVS1Prediction.PVS1_Moderate: AutoACMGStrength.PathogenicModerate,
            PVS1Prediction.PVS1_Supporting: AutoACMGStrength.PathogenicSupporting,
            PVS1Prediction.NotPVS1: AutoACMGStrength.PathogenicVeryStrong,
            PVS1Prediction.UnsupportedConsequence: AutoACMGStrength.PathogenicVeryStrong,
            PVS1Prediction.NotSet: AutoACMGStrength.PathogenicVeryStrong,
        }
        return AutoACMGCriteria(
            name="PVS1",
            prediction=(
                AutoACMGPrediction.Met
                if pred
                in [
                    PVS1Prediction.PVS1,
                    PVS1Prediction.PVS1_Strong,
                    PVS1Prediction.PVS1_Moderate,
                    PVS1Prediction.PVS1_Supporting,
                ]
                else AutoACMGPrediction.NotMet
            ),
            strength=evidence_strength_mapping[pred],
            summary=comment,
            description=f"Prediction path: {PVS1PredictionPathMapping[path]}",
        )
