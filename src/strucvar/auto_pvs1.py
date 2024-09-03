"""PVS1 criteria for Structural Variants (StrucVar)."""

from typing import Dict, List, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGStrength,
    AutoACMGStrucVarData,
)
from src.defs.auto_pvs1 import PVS1Prediction, PVS1PredictionPathMapping, PVS1PredictionStrucVarPath
from src.defs.exceptions import MissingDataError
from src.defs.mehari import Exon
from src.defs.strucvar import StrucVar, StrucVarType
from src.utils import AutoACMGHelper


class StrucVarHelper(AutoACMGHelper):
    """Helper methods for PVS1 criteria for Structural Variants (StrucVar)."""

    def __init__(self):
        super().__init__()
        #: Comment to store the prediction explanation.
        self.comment_pvs1: str = ""

    def full_gene_del(self, strucvar: StrucVar, exons: List[Exon]) -> bool:
        """
        Check if the variant is a full gene deletion.

        Args:
            strucvar: The structural variant.
            exons: The exons of the gene.

        Returns:
            True if the variant is a full gene deletion, False otherwise.

        Raises:
            MissingDataError: If exons are not available.
        """
        if not exons:
            raise MissingDataError(
                "Exons are not available. Cannot determine if the variant is a full gene deletion."
            )

        gene_start = min(
            exons[0].altStartI, exons[0].altEndI, exons[-1].altStartI, exons[-1].altEndI
        )
        gene_end = max(exons[0].altStartI, exons[0].altEndI, exons[-1].altStartI, exons[-1].altEndI)
        self.comment_pvs1 += f"Gene start: {gene_start}, gene end: {gene_end}."
        return strucvar.start <= gene_start and strucvar.stop >= gene_end

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

    def verify_pvs1(
        self, strucvar: StrucVar, var_data: AutoACMGStrucVarData
    ) -> Tuple[PVS1Prediction, PVS1PredictionStrucVarPath, str]:
        """Make the PVS1 prediction.

        The prediction is based on the PVS1 criteria for structural variants.
        """
        if strucvar.sv_type == StrucVarType.DEL:
            self.comment_pvs1 = "Analysing the deletion variant. => "
            if self.full_gene_del(strucvar, var_data.exons):
                self.prediction = PVS1Prediction.PVS1
                self.prediction_path = PVS1PredictionStrucVarPath.DEL1
            elif self.del_disrupt_rf() and self.undergo_nmd():
                self.comment_pvs1 += " =>"
                if self.in_bio_relevant_tsx():
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL2
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL3
            elif self.del_disrupt_rf() and not self.undergo_nmd():
                self.comment_pvs1 += " =>"
                if self.crit4prot_func():
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL4
                else:
                    self.comment_pvs1 += " =>"
                    if self.lof_freq_in_pop() or not self.in_bio_relevant_tsx():
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionStrucVarPath.DEL5_1
                    else:
                        self.comment_pvs1 += " =>"
                        if self.lof_rm_gt_10pct_of_prot():
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL6_1
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL7_1
            else:
                self.comment_pvs1 += " =>"
                if self.crit4prot_func():
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL8
                else:
                    self.comment_pvs1 += " =>"
                    if self.lof_freq_in_pop() or not self.in_bio_relevant_tsx():
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionStrucVarPath.DEL5_2
                    else:
                        self.comment_pvs1 += " =>"
                        if self.lof_rm_gt_10pct_of_prot():
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL6_2
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL7_2

        elif strucvar.sv_type == StrucVarType.DUP:
            self.comment_pvs1 = "Analysing the duplication variant. => "
            if self.proven_in_tandem():
                self.comment_pvs1 += " =>"
                if self.dup_disrupt_rf() and self.undergo_nmd():
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DUP1
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DUP2_1
            elif self.presumed_in_tandem():
                self.comment_pvs1 += " =>"
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
            logger.error("Unsupported structural variant type: {}", strucvar.sv_type)
            self.comment_pvs1 = "Unsupported structural variant type."

        return self.prediction, self.prediction_path, self.comment_pvs1

    def predict_pvs1(self, strucvar: StrucVar, var_data: AutoACMGStrucVarData) -> AutoACMGCriteria:
        """Predict the PVS1 criteria for structural variants."""
        pred, path, comment = self.verify_pvs1(strucvar, var_data)
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
