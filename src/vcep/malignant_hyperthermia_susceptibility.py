"""
Predictor for Malignant Hyperthermia Susceptibility VCEP.
Included gene: RYR1 (HGNC:10483).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN012
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

#: VCEP specification for Malignant Hyperthermia Susceptibility.
SPEC: VcepSpec = VcepSpec(
    identifier="GN012",
    version="2.0.0",
)

PM1_CLUSTER = {
    "HGNC:10483": {
        "moderate": {
            "domains": [
                (1, 552),  # N-terminal region
                (2101, 2458),  # Central region
            ],
        },
        "supporting": {
            "domains": [
                (1, 552),  # N-terminal region (if PS1/PM5 applicable)
                (2101, 2458),  # Central region (if PS1/PM5 applicable)
                (4631, 4991),  # C-terminal region
            ],
        },
    }
}


class MalignantHyperthermiaPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override PM1 to specify critical domains for Malignant Hyperthermia Susceptibility."""
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check moderate level domains
        for start, end in gene_cluster.get("moderate", {}).get("domains", []):
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

        # Check supporting level domains
        for start, end in gene_cluster.get("supporting", {}).get("domains", []):
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

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """
        Override PP2 and BP1 for Malignant Hyperthermia Susceptibility to return not applicable
        status.
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
