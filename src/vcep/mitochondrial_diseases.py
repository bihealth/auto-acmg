"""
Predictor for Mitochondrial Diseases VCEP.
Included genes:
ETHE1 (HGNC:23287),
PDHA1 (HGNC:8806),
POLG (HGNC:9179),
SLC19A3 (HGNC:16266),
other Mitochondrial Disease genes.  (Mostly not implemented in this snippet)
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN014
https://cspec.genome.network/cspec/ui/svi/doc/GN015
"""

from typing import List, Tuple

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

#: VCEP specifications for Mitochondrial Diseases.
SPECs: List[VcepSpec] = [
    VcepSpec(
        identifier="GN014",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN015",
        version="1.0.0",
    ),
]


class MitochondrialDiseasesPredictor(DefaultSeqVarPredictor):
    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to return not applicable status for ETHE1, PDHA1, POLG, and SLC19A3.
        """
        logger.info("Predict PM1")

        if var_data.hgnc_id in [
            "HGNC:23287",  # ETHE1
            "HGNC:8806",  # PDHA1
            "HGNC:9179",  # POLG
            "HGNC:16266",  # SLC19A3
        ]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=f"PM1 is not applicable for {var_data.hgnc_id}.",
            )

        return super().predict_pm1(seqvar, var_data)

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.00005
        var_data.thresholds.ba1_benign = 0.001
        var_data.thresholds.bs1_benign = 0.0005
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to return not applicable status for PP2 and BP1."""
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
        var_data.thresholds.revel_pathogenic = 0.75
        var_data.thresholds.revel_benign = 0.15
        return super().predict_pp3bp4(seqvar, var_data)
