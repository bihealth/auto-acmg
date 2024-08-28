"""
Predictor for InSIGHT Hereditary Colorectal Cancer/Polyposis VCEP.
Included genes:
APC (HGNC:583),
MLH1 (HGNC:7127),
MSH2 (HGNC:7325),
MSH6 (HGNC:7329),
PMS2 (HGNC:9122).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN089
https://cspec.genome.network/cspec/ui/svi/doc/GN115
https://cspec.genome.network/cspec/ui/svi/doc/GN137
https://cspec.genome.network/cspec/ui/svi/doc/GN138
https://cspec.genome.network/cspec/ui/svi/doc/GN139
"""

from typing import Tuple

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar


class InsightColorectalCancerPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override PM1 prediction for InSIGHT Hereditary Colorectal Cancer/Polyposis."""
        logger.info("Predict PM1")

        if var_data.hgnc_id in [
            "HGNC:583",  # APC
            "HGNC:7127",  # MLH1
            "HGNC:7325",  # MSH2
            "HGNC:7329",  # MSH6
            "HGNC:9122",  # PMS2
        ]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=f"PM1 is not applicable for {var_data.hgnc_id}.",
            )

        return super().predict_pm1(seqvar, var_data)

    def predict_pm4bp3(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        return (
            AutoACMGCriteria(
                name="PM4",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=(
                    "PM4 is not applicable for the InSIGHT Hereditary Colorectal Cancer/Polyposis "
                    "VCEP."
                ),
            ),
            AutoACMGCriteria(
                name="BP3",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.BenignSupporting,
                summary=(
                    "BP3 is not applicable for the InSIGHT Hereditary Colorectal Cancer/Polyposis "
                    "VCEP."
                ),
            ),
        )

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override donor and acceptor positions for InSIGHT Hereditary Colorectal Cancer/Polyposis
        VCEP.
        """
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
