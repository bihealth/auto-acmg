"""
Predictor for Thrombosis VCEP.
Included gene: SERPINC1 (HGNC:775).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN084
"""

from typing import Tuple

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

#: VCEP specification for Thrombosis.
SPEC: VcepSpec = VcepSpec(
    identifier="GN084",
    version="1.0.0",
)

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


class ThrombosisPredictor(DefaultSeqVarPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override PM1 to specify critical domains for SERPINC1."""
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

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.00002
        var_data.thresholds.ba1_benign = 0.002
        var_data.thresholds.bs1_benign = 0.0002
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 is not applicable for SERPINC1."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to include VCEP-specific logic for SERPINC1."""
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
