"""
Predictor for ACADVL VCEP.
Included gene: ACADVL (HGNC:92).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN021
"""

from typing import Tuple

from loguru import logger

from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGSeqVarResult,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

#: VCEP specification for ACADVL.
SPEC: VcepSpec = VcepSpec(
    identifier="GN021",
    version="1.0.0",
)

PM1_CLUSTER = [
    (214, 223),  # Nucleotide and substrate binding
    (249, 251),  # Nucleotide and substrate binding
    (460, 466),  # Nucleotide and substrate binding
    (481, 516),  # Membrane binding
    (1, 40),  # Mitochondrial signal peptide
]


class ACADVLPredictor(DefaultSeqVarPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override predict_pm1 to return a not applicable status for PM1."""
        logger.info("Predict PM1")
        # Check if variant falls within critical regions
        for start, end in PM1_CLUSTER:
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

    def _bs2_not_applicable(self, var_data: AutoACMGSeqVarData) -> bool:
        """BS2 is not applicable for ACADVL."""
        return True

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.001
        var_data.thresholds.ba1_benign = 0.007
        var_data.thresholds.bs1_benign = 0.0035
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 is not applicable for ACADVL."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
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
