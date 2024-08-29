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

from typing import List, Tuple

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.seqvar import SeqVar

#: VCEP specifications for Cerebral Creatine Deficiency Syndromes.
SPECs: List[VcepSpec] = [
    VcepSpec(
        identifier="GN025",
        version="1.1.0",
    ),
    VcepSpec(
        identifier="GN026",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN027",
        version="1.1.0",
    ),
]


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

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGData) -> bool:
        """Override _bp3_not_applicable for Cerebral Creatine Deficiency Syndromes."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """
        Override predict_pp2bp1 to include VCEP-specific logic for Cerebral Creatine Deficiency
        Syndromes.
        """
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

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override donor and acceptor positions for Cerebral Creatine Deficiency Syndromes VCEP."""
        if var_data.hgnc_id == "HGNC:4136":
            var_data.thresholds.bp7_donor = 7
            var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
