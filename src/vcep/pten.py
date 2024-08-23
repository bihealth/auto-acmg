"""
Predictor for PTEN VCEP.
Included gene: PTEN (HGNC:9588).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN003
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

PM1_CLUSTER = {
    "HGNC:9588": {  # PTEN
        "moderate":
        # catalytic motifs
        list(range(90, 95))
        + list(range(123, 131))
        + list(range(166, 169)),
    }
}


class PTENPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for PTEN.
        """
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check if the variant falls within the specified catalytic motifs
        if var_data.prot_pos in gene_cluster.get("moderate", []):
            comment = (
                f"Variant affects a critical residue within the catalytic motifs of PTEN "
                f"at position {var_data.prot_pos}. PM1 is met at the Moderate level."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=comment,
            )

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary="Variant does not meet the PM1 criteria for PTEN.",
        )
