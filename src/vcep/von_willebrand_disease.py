"""
Predictor for von Willebrand Disease VCEP.
Included gene: VWF (HGNC:12726).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN081
https://cspec.genome.network/cspec/ui/svi/doc/GN090
"""

from typing import Tuple

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.seqvar import SeqVar

#: VCEP specification for von Willebrand Disease.
SPEC: VcepSpec = VcepSpec(
    identifier="GN081",
    version="1.0.0",
)


class VonWillebrandDiseasePredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
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

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.0001
        var_data.thresholds.ba1_benign = 0.1
        var_data.thresholds.bs1_benign = 0.0
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGData) -> bool:
        """BP3 it not applicable."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
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
