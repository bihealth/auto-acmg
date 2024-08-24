"""
Predictor for Familial Hypercholesterolemia VCEP.
Included gene: LDLR (HGNC:6547).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN013
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

# fmt: off
PM1_CLUSTER_LDLR = {
    "HGNC:6547": {  # LDLR
        "residues": [
            27,  34,  39,  46,  52,  63,  68,  75,  82,  89,
            95,  104, 109, 116, 121, 128, 134, 143, 148, 155,
            160, 167, 173, 184, 197, 204, 209, 216, 231, 236,
            243, 248, 255, 261, 270, 276, 284, 289, 296, 302,
            313, 318, 325, 329, 338, 340, 352, 358, 364, 368,
            377, 379, 392, 667, 677, 681, 696, 698, 711,
        ]
    },
}



class FamilialHypercholesterolemiaPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override predict_pm1 to include VCEP-specific logic for LDLR."""
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER_LDLR.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

         # Check if the variant affects a critical residue
        if var_data.prot_pos in gene_cluster["residues"]:
            comment = (
                f"Variant affects a critical cysteine residue in LDLR at position "
                f"{var_data.prot_pos}, which is involved in disulfide bond formation. PM1 is met."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=comment,
            )

        # Check if the variant in the 4th exon
        affected_exon = self._get_affected_exon(var_data, seqvar)
        if affected_exon == 4:
            comment = f"Variant affects the 4th exon in LDLR. PM1 is met."
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=comment,
            )

        # If no criteria match
        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary=f"Variant does not meet the PM1 criteria for LDLR.",
        )
