"""
Predictor for CDH1 VCEP.
Included gene: CDH1 (HGNC:1748).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN007
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar


class CDH1Predictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override predict_pm1 to include VCEP-specific logic for CDH1."""
        logger.info("Predict PM1")

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotApplicable,
            strength=AutoACMGStrength.PathogenicSupporting,
            summary="PM1 is not applicable for CDH1.",
        )

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override donor and acceptor positions for CDH1 VCEP."""
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
