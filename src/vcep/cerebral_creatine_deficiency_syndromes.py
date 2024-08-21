"""
Predictor for Cerebral Creatine Deficiency Syndromes VCEP.
Included genes:
GATM (HGNC:4175),
GAMT (HGNC:4136),
SLC6A8 (HGNC:11055).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN025
https://cspec.genome.network/cspec/ui/svi/doc/GN026
https://cspec.genome.network/cspec/ui/svi/doc/GN027
"""

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar


class CerebralCreatineDeficiencySyndromesPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Cerebral Creatine Deficiency
        Syndromes.
        """
        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotApplicable,
            strength=AutoACMGStrength.PathogenicModerate,
            summary="PM1 is not applicable for Cerebral Creatine Deficiency Syndromes.",
        )