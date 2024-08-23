"""
Predictor for Phenylketonuria (PKU) VCEP.
Included gene: PAH (HGNC:8582).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN006
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

PM1_CLUSTER_PKU = {
    "HGNC:8582": {  # PAH
        "residues": [
            # active site residues
            138,  # Tyr138
            158,  # Arg158
            245,  # Val245
            268,  # Tyr268
            278,  # Thr278
            279,  # Pro279
            289,  # Glu289
            300,  # Ala300
            315,  # Asp315
            331,  # Phe331
            345,  # Ala345
            346,  # Gly346
            349,  # Ser349
            377,  # Tyr377
        ]
        # substrate binding residues
        + list(range(46, 49))
        + list(range(63, 70))
        # cofactor binding residues
        + [
            285,  # His285
            290,  # His290
            330,  # Glu330
        ]
        + list(range(246, 267))
        + list(range(280, 284))
        + list(range(322, 327))
        + list(range(377, 380)),  # Additional ranges for cofactor binding
    }
}


class PKUPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Phenylketonuria (PKU).
        """
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER_PKU.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check moderate level criteria
        if var_data.prot_pos in gene_cluster.get("residues", []):
            comment = (
                f"Variant affects a residue at position {var_data.prot_pos} "
                f"in {var_data.hgnc_id}, which is critical for PKU."
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
            summary="Variant does not meet the PM1 criteria for PAH.",
        )
