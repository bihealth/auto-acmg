"""
Predictor for Platelet Disorders VCEP.
Included genes:
ITGA2B (HGNC:6138),
ITGB3 (HGNC:6156).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN011
"""

from typing import Tuple

from loguru import logger

from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

#: VCEP specification for Platelet Disorders.
SPEC: VcepSpec = VcepSpec(
    identifier="GN011",
    version="2.1.0",
)


class PlateletDisordersPredictor(DefaultSeqVarPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """
        Override PM1 prediction for Platelet Disorders to return a not met status for ITGA2B and
        ITGB3.
        """
        logger.info("Predict PM1")

        if var_data.hgnc_id in ["HGNC:6138", "HGNC:6156"]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotMet,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary="PM1 is not met for ITGA2B and ITGB3 due to genes being highly polymorphic.",
            )

        return super().predict_pm1(seqvar, var_data)

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.0001
        var_data.thresholds.ba1_benign = 0.0024
        var_data.thresholds.bs1_benign = 0.00158
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

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

    def predict_pp3bp4(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Use REVEL for PP3 and BP4 for Epilepsy Sodium Channel."""
        var_data.thresholds.pp3bp4_strategy = "revel"
        var_data.thresholds.revel_pathogenic = 0.7
        var_data.thresholds.revel_benign = 0.25
        return super().predict_pp3bp4(seqvar, var_data)
