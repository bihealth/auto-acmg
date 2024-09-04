"""
Predictor for Glaucoma VCEP.
Included gene: MYOC (HGNC:7610).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN019
"""

from typing import Optional, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    PS1PM5,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.exceptions import MissingDataError
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

#: VCEP specification for MYOC.
SPEC = VcepSpec(
    identifier="GN019",
    version="1.1.0",
)


class GlaucomaPredictor(DefaultSeqVarPredictor):

    def predict_pvs1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """PVS1 is not applicable."""
        logger.info("Predict PVS1")
        return AutoACMGCriteria(
            name="PVS1",
            prediction=AutoACMGPrediction.NotApplicable,
            strength=AutoACMGStrength.PathogenicVeryStrong,
            summary="PVS1 is not applicable for the gene.",
        )

    def verify_ps1pm5(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PS1PM5], str]:
        """Override PS1/PM5 for MYOC."""
        self.prediction_ps1pm5, self.comment_ps1pm5 = super().verify_ps1pm5(seqvar, var_data)
        var_data.thresholds.spliceAI_acceptor_gain = 0.2
        var_data.thresholds.spliceAI_acceptor_loss = 0.2
        var_data.thresholds.spliceAI_donor_gain = 0.2
        var_data.thresholds.spliceAI_donor_loss = 0.2
        if self.prediction_ps1pm5 and self._affect_splicing(var_data):
            self.prediction_ps1pm5.PS1 = False
            self.comment_ps1pm5 = "Variant affects splicing. PS1 is not applicable."
            self.prediction_ps1pm5.PM5 = False
            self.comment_ps1pm5 = "Variant affects splicing. PM5 is not applicable."
        return self.prediction_ps1pm5, self.comment_ps1pm5

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override predict_pm1 to return a not applicable status for PM1."""
        logger.info("Predict PM1")

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotApplicable,
            strength=AutoACMGStrength.PathogenicSupporting,
            summary="PM1 is not applicable for MYOC.",
        )

    def _bs2_not_applicable(self, var_data: AutoACMGSeqVarData) -> bool:
        """BS2 is not applicable for Glaucoma VCEP."""
        return True

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.0001
        var_data.thresholds.ba1_benign = 0.01
        var_data.thresholds.bs1_benign = 0.001
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 is not applicable for Glaucoma VCEP."""
        return True

    def _is_conserved(self, var_data: AutoACMGSeqVarData) -> bool:
        """
        Predict if the variant is conserved.

        Check if the variant is conserved using the GERP score.

        Args:
            variant_info: The variant information.

        Returns:
            bool: True if the variant is conserved, False otherwise.

        Raises:
            MissingDataError: If the GERP score is missing
        """
        if not var_data.scores.cadd.gerp:
            raise MissingDataError("GERP score is missing.")
        return var_data.scores.cadd.gerp >= var_data.thresholds.gerp

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to return a not applicable status for PP2 and BP1."""
        return (
            AutoACMGCriteria(
                name="PP2",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary="PP2 is not applicable for the gene.",
            ),
            AutoACMGCriteria(
                name="BP1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary="BP1 is not applicable for the gene.",
            ),
        )

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Change BP7 thresholds for Glaucoma VCEP."""
        var_data.thresholds.spliceAI_acceptor_gain = 0.2
        var_data.thresholds.spliceAI_acceptor_loss = 0.2
        var_data.thresholds.spliceAI_donor_gain = 0.2
        var_data.thresholds.spliceAI_donor_loss = 0.2
        return super().predict_bp7(seqvar, var_data)
