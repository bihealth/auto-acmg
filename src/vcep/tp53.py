"""
Predictor for TP53 VCEP.
Included gene: TP53 (HGNC:11998).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN009
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

#: VCEP specification for TP53.
SPEC: VcepSpec = VcepSpec(
    identifier="GN009",
    version="2.0.0",
)

# Define the critical regions for PM1 consideration in TP53
PM1_CLUSTER = {
    # "NM_000546.4": {
    "HGNC:11998": {
        "residues": [175, 245, 248, 249, 273, 282],
    }
}


class TP53Predictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override PM1 prediction to specify critical residues for TP53."""
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check if variant is in the moderate cluster
        if var_data.prot_pos in gene_cluster["residues"]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=f"Variant affects a critical residue in TP53 at position {var_data.prot_pos}. "
                f"PM1 is met at the Moderate level.",
            )

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary="Variant does not meet the PM1 criteria for TP53.",
        )

    def predict_pm4bp3(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override PM4 and BP3 to return not applicable for TP53 VCEP."""
        return (
            AutoACMGCriteria(
                name="PM4",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicModerate,
                summary="PM4 is not applicable for TP53 VCEP.",
            ),
            AutoACMGCriteria(
                name="BP3",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.BenignSupporting,
                summary="BP3 is not applicable for TP53 VCEP.",
            ),
        )

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to return not applicable status for TP53."""
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
        """Override donor and acceptor positions for TP53 VCEP."""
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
