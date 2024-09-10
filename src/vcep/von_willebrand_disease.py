"""
Predictor for von Willebrand Disease VCEP.
Included gene: VWF (HGNC:12726).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN081
https://cspec.genome.network/cspec/ui/svi/doc/GN090
"""

from typing import Optional, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    PP3BP4,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.exceptions import AutoAcmgBaseException
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

#: VCEP specification for von Willebrand Disease.
SPEC: VcepSpec = VcepSpec(
    identifier="GN081",
    version="1.0.0",
)


class VonWillebrandDiseasePredictor(DefaultSeqVarPredictor):
    def predict_pvs1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """PVS1 is not applicable."""
        logger.info("Predict PVS1")
        return AutoACMGCriteria(
            name="PVS1",
            prediction=AutoACMGPrediction.NotApplicable,
            strength=AutoACMGStrength.PathogenicVeryStrong,
            summary="PVS1 is not applicable for the gene.",
        )

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to return a not applicable status for PM1.
        """
        logger.info("Predict PM1")

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotApplicable,
            strength=AutoACMGStrength.PathogenicSupporting,
            summary=f"PM1 is not applicable for {var_data.hgnc_id}.",
        )

    def _bs2_not_applicable(self, var_data: AutoACMGSeqVarData) -> bool:
        """BS2 is not applicable for von Willebrand Disease VCEP."""
        return True

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.0001
        var_data.thresholds.ba1_benign = 0.1
        var_data.thresholds.bs1_benign = 0.01
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 it not applicable."""
        return True

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

    def verify_pp3bp4(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PP3BP4], str]:
        """Predict PP3 and BP4 criteria."""
        self.prediction_pp3bp4 = PP3BP4()
        self.comment_pp3bp4 = ""
        try:
            score = "revel"
            var_data.thresholds.revel_pathogenic = 0.644
            var_data.thresholds.revel_benign = 0.29
            self.prediction_pp3bp4.PP3 = self._is_pathogenic_score(
                var_data,
                (score, getattr(var_data.thresholds, f"{score}_pathogenic")),
            )
            self.prediction_pp3bp4.BP4 = self._is_benign_score(
                var_data,
                (score, getattr(var_data.thresholds, f"{score}_benign")),
            )

            var_data.thresholds.spliceAI_acceptor_gain = 0.5
            var_data.thresholds.spliceAI_acceptor_loss = 0.5
            var_data.thresholds.spliceAI_donor_gain = 0.5
            var_data.thresholds.spliceAI_donor_loss = 0.5
            self.prediction_pp3bp4.PP3 = self.prediction_pp3bp4.PP3 or self._affect_spliceAI(
                var_data
            )
            var_data.thresholds.spliceAI_acceptor_gain = 0.1
            var_data.thresholds.spliceAI_acceptor_loss = 0.1
            var_data.thresholds.spliceAI_donor_gain = 0.1
            var_data.thresholds.spliceAI_donor_loss = 0.1
            self.prediction_pp3bp4.BP4 = self.prediction_pp3bp4.BP4 and not self._affect_spliceAI(
                var_data
            )

        except AutoAcmgBaseException as e:
            self.comment_pp3bp4 = f"An error occurred during prediction. Error: {e}"
            self.prediction_pp3bp4 = None
        return self.prediction_pp3bp4, self.comment_pp3bp4
