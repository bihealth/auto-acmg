"""
Predictor for ENIGMA BRCA1 and BRCA2 VCEP.
Included genes:
BRCA1 (HGNC:1100),
BRCA2 (HGNC:1101).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN092
https://cspec.genome.network/cspec/ui/svi/doc/GN097
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar


class ENIGMAPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override PM1 prediction for ENIGMA BRCA1 and BRCA2."""
        logger.info("Predict PM1")

        if var_data.hgnc_id in ["HGNC:1100", "HGNC:1101"]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicModerate,
                summary="PM1 is not applicable for BRCA1 and BRCA2.",
            )

        return super().predict_pm1(seqvar, var_data)
