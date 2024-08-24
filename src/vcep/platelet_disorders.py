"""
Predictor for Platelet Disorders VCEP.
Included genes:
ITGA2B (HGNC:6138),
ITGB3 (HGNC:6156).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN011
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar


class PlateletDisordersPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override PM1 prediction for Platelet Disorders."""
        logger.info("Predict PM1")

        if var_data.hgnc_id in ["HGNC:6138", "HGNC:6156"]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotMet,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary="PM1 is not met for ITGA2B and ITGB3 due to genes being highly polymorphic.",
            )

        return super().predict_pm1(seqvar, var_data)
