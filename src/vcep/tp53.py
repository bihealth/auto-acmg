"""
Predictor for TP53 VCEP.
Included gene: TP53 (HGNC:11998).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN009
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

# Define the critical regions for PM1 consideration in TP53
PM1_CLUSTER = {
    # "NM_000546.4": {
    "HGNC:11998": {
        "residues": [175, 245, 248, 249, 273, 282],
    }
}


class TP53Predictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for TP53.
        """
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
