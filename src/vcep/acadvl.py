"""
Predictor for ACADVL VCEP.
Included gene: ACADVL (HGNC:92).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN021
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

PM1_CLUSTER_REGIONS = [
    (214, 223),  # Nucleotide and substrate binding
    (249, 251),  # Nucleotide and substrate binding
    (460, 466),  # Nucleotide and substrate binding
    (481, 516),  # Membrane binding
    (1, 40),  # Mitochondrial signal peptide
]


class ACADVLPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override predict_pm1 to include VCEP-specific logic for ACADVL."""
        logger.info("Predict PM1")
        # Check if variant falls within critical regions
        for start, end in PM1_CLUSTER_REGIONS:
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant falls within a critical region for ACADVL between positions "
                    f"{start}-{end}. PM1 is met."
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
            summary="Variant does not fall within any critical region for ACADVL. PM1 is not met.",
        )