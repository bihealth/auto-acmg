"""
Predictor for Myeloid Malignancy VCEP.
Included gene: RUNX1 (HGNC:10471).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN008
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

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
        "supporting": [  # Other residues within the RHD (89-204)
            *range(89, 205),
        ],
    }
}


class MyeloidMalignancyPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Myeloid Malignancy.
        """
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
