"""PVS1 criteria for Structural Variants (StrucVar)."""

from typing import Optional

from loguru import logger

from src.core.config import Config
from src.defs.auto_pvs1 import PVS1Prediction, PVS1PredictionStrucVarPath
from src.defs.strucvar import StrucVar, StrucVarType


class StrucVarPVS1Helper:
    """Helper methods for PVS1 criteria for structural variants."""

    def __init__(self, *, config: Optional[Config] = None):
        #: Configuration to use.
        self.config: Config = config or Config()
        #: Comment to store the prediction explanation.
        self.comment: str = ""

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

    @staticmethod
    def undergo_nmd() -> bool:
        """Check if the deletion is expected to undergo NMD."""
        return False

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


class StrucVarPVS1(StrucVarPVS1Helper):
    """PVS1 prediction for structural variants."""

    def __init__(self, variant: StrucVar, *, config: Optional[Config] = None):
        super().__init__()

        # === Attributes ===
        self.variant = variant

        # === Prediction attributes ===
        self.prediction: PVS1Prediction = PVS1Prediction.NotSet
        self.prediction_path: PVS1PredictionStrucVarPath = PVS1PredictionStrucVarPath.NotSet

    def initialize(self):
        """Initialize the PVS1 prediction."""
        pass

    def verify_PVS1(self):
        """Verify PVS1 prediction."""
        if self.variant.sv_type == StrucVarType.DEL:
            if self.full_gene_del():
                self.prediction = PVS1Prediction.PVS1
                self.prediction_path = PVS1PredictionStrucVarPath.DEL1
            elif self.del_disrupt_rf() and self.undergo_nmd():
                if self.in_bio_relevant_tsx():
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL2
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL3
            elif self.del_disrupt_rf() and not self.undergo_nmd():
                if self.crit4prot_func():
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL4
                else:
                    if self.lof_freq_in_pop() or not self.in_bio_relevant_tsx():
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionStrucVarPath.DEL5_1
                    else:
                        if self.lof_rm_gt_10pct_of_prot():
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL6_1
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL7_1
            else:
                if self.crit4prot_func():
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL8
                else:
                    if self.lof_freq_in_pop() or not self.in_bio_relevant_tsx():
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionStrucVarPath.DEL5_2
                    else:
                        if self.lof_rm_gt_10pct_of_prot():
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL6_2
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL7_2

        elif self.variant.sv_type == StrucVarType.DUP:
            if self.proven_in_tandem():
                if self.dup_disrupt_rf() and self.undergo_nmd():
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DUP1
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DUP2_1
            elif self.presumed_in_tandem():
                if self.dup_disrupt_rf() and self.undergo_nmd():
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

        return self.prediction, self.prediction_path, self.comment
