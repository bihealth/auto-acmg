"""
Predictor for Heriditary Breast, Ovarian and Pancreatic Cancer VCEP.
Included genes:
ATM (HGNC:795),
PALB2 (HGNC:26144).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN020
https://cspec.genome.network/cspec/ui/svi/doc/GN077
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar


class HBOPCPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Heriditary Breast, Ovarian and
        Pancreatic Cancer.
        """
        logger.info("Predict PM1")

        if var_data.hgnc_id in [
            "HGNC:795",  # ATM
            "HGNC:26144",  # PALB2
        ]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=f"PM1 is not applicable for {var_data.hgnc_id}.",
            )

        return super().predict_pm1(seqvar, var_data)

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override donor and acceptor positions for ATM and PALB2."""
        if var_data.hgnc_id == "HGNC:26144":
            var_data.thresholds.bp7_donor = 7
            var_data.thresholds.bp7_acceptor = 21
        elif var_data.hgnc_id == "HGNC:795":
            var_data.thresholds.bp7_donor = 7
            var_data.thresholds.bp7_acceptor = 40
        return super().predict_bp7(seqvar, var_data)
