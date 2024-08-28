"""
Predictor for Thrombosis VCEP.
Included gene: SERPINC1 (HGNC:775).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN084
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

# fmt: off
PM1_CLUSTER = {
    "HGNC:775": {  # SERPINC1
        "residues": [
            # Cysteine residues involved in disulfide bridges
            40, 53, 127, 160, 279, 462,
            # Heparin binding site residues
            39, 56, 73, 79,
            # Reactive site residues
            414, 416,
        ],
    }
}


class ThrombosisPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Thrombosis.
        """
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check if variant is in the cluster
        if var_data.prot_pos in gene_cluster["residues"]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicModerate,
                summary="Variant meets the PM1 criteria for SERPINC1.",
            )

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary="Variant does not meet the PM1 criteria for SERPINC1.",
        )

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGData) -> bool:
        """Override BP3 for thrombosis VCEP."""
        return True
