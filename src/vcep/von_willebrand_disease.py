"""
Predictor for von Willebrand Disease VCEP.
Included gene: VWF (HGNC:12726).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN081
https://cspec.genome.network/cspec/ui/svi/doc/GN090
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar


class VonWillebrandDiseasePredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for von Willebrand Disease.
        """
        logger.info("Predict PM1")

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotApplicable,
            strength=AutoACMGStrength.PathogenicSupporting,
            summary=f"PM1 is not applicable for {var_data.hgnc_id}.",
        )
