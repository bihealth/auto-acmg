"""
Predictor for Malignant Hyperthermia Susceptibility VCEP.
Included gene: RYR1 (HGNC:10483).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN012
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

PM1_CLUSTER = {
    "HGNC:10483": {
        "moderate": {
            "regions": [
                (1, 552),  # N-terminal region
                (2101, 2458),  # Central region
            ],
        },
        "supporting": {
            "regions": [
                (1, 552),  # N-terminal region (if PS1/PM5 applicable)
                (2101, 2458),  # Central region (if PS1/PM5 applicable)
                (4631, 4991),  # C-terminal region
            ],
        },
    }
}


class MalignantHyperthermiaPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Malignant Hyperthermia Susceptibility.
        """
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check moderate level regions
        for start, end in gene_cluster.get("moderate", {}).get("regions", []):
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant falls within a critical region in RYR1 between positions {start}-{end}. "
                    f"PM1 is met at the Moderate level."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicModerate,
                    summary=comment,
                )

        # Check supporting level regions
        for start, end in gene_cluster.get("supporting", {}).get("regions", []):
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant falls within a critical region in RYR1 between positions {start}-{end}. "
                    f"PM1 is met at the Supporting level."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicSupporting,
                    summary=comment,
                )

        # If no criteria match
        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary="Variant does not meet the PM1 criteria for Malignant Hyperthermia Susceptibility.",
        )
