"""
Predictor for Myeloid Malignancy VCEP.
Included gene: RUNX1 (HGNC:10471).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN008
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

#: VCEP specification for Myeloid Malignancy.
SPEC: VcepSpec = VcepSpec(
    identifier="GN008",
    version="2.0.0",
)

PM1_CLUSTER = {
    "HGNC:10471": {  # RUNX1
        "moderate": [  # Specific residues within the RHD (Runt Homology Domain)
            107,  # R107
            110,  # K110
            134,  # A134
            162,  # R162
            166,  # R166
            167,  # S167
            169,  # R169
            170,  # G170
            194,  # K194
            196,  # T196
            198,  # D198
            201,  # R201
            204,  # R204
        ],
        "supporting": list(range(89, 205)),  # Other residues within the RHD (89-204)
    }
}


class MyeloidMalignancyPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override PM1 to specify critical residues for Myeloid Malignancy VCEP."""
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check moderate level criteria
        if var_data.prot_pos in gene_cluster["moderate"]:
            comment = (
                f"Variant affects a critical residue within the Runt Homology Domain (RHD) "
                f"at position {var_data.prot_pos} in RUNX1. PM1 is met at the Moderate level."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=comment,
            )

        # Check supporting level criteria
        if var_data.prot_pos in gene_cluster["supporting"]:
            comment = (
                f"Variant affects a residue within the Runt Homology Domain (RHD) "
                f"at position {var_data.prot_pos} in RUNX1. PM1 is met at the Supporting level."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=comment,
            )

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary="Variant does not meet the PM1 criteria for RUNX1.",
        )

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGData) -> bool:
        """BP3 is not applicable for Myeloid Malignancy VCEP."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to return not applicable status for RUNX1."""
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
        """Change BP7 thresholds for Myeloid Malignancy VCEP."""
        var_data.thresholds.spliceAI_acceptor_gain = 0.2
        var_data.thresholds.spliceAI_acceptor_loss = 0.2
        var_data.thresholds.spliceAI_donor_gain = 0.2
        var_data.thresholds.spliceAI_donor_loss = 0.2
        var_data.thresholds.phyloP100 = 2.0
        return super().predict_bp7(seqvar, var_data)
