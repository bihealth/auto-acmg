"""
Predictor for Leber Congenital Amaurosis/early onset Retinal Dystrophy VCEP.
Included gene: RPE65 (HGNC:10294)
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN120
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

PM1_CLUSTER = {
    "HGNC:10294": {
        "residues": [180, 182, 241, 313, 417, 527] + list(range(107, 126)),
    }
}


class LeberCongenitalAmaurosisPredictor(DefaultPredictor):
    """
    Predictor for Leber Congenital Amaurosis/early onset Retinal Dystrophy VCEP.
    """

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Leber Congenital Amaurosis/early
        onset Retinal Dystrophy.
        """
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        if var_data.prot_pos in gene_cluster["residues"]:
            comment = (
                f"Variant affects a residue in RPE65 at position {var_data.prot_pos}, "
                f"which is associated with Leber Congenital Amaurosis/early onset Retinal "
                "Dystrophy."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=comment,
            )

        # If no criteria are met, return Not Met
        comment = (
            f"Variant does not meet the PM1 criteria for Leber Congenital Amaurosis/early onset "
            "Retinal Dystrophy."
        )
        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary=comment,
        )
