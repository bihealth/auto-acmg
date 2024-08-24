"""
Predictor for Lysosomal Diseases VCEP.
Included gene: GAA (HGNC:4065).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN010
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

PM1_CLUSTER = {
    "HGNC:4065": {
        "residues": [
            282,  # D282
            376,  # W376
            404,  # D404
            405,  # L405
            441,  # I441
            481,  # W481
            516,  # W516
            518,  # D518
            519,  # M519
            600,  # R600
            613,  # W613
            616,  # D616
            618,  # W618
            649,  # F649
            650,  # L650
            674,  # H674
        ]
    }
}


class LysosomalDiseasesPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Lysosomal Diseases.
        """
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        if var_data.prot_pos in gene_cluster["residues"]:
            comment = (
                f"Variant affects a residue in GAA at position {var_data.prot_pos}, "
                f"which is associated with Lysosomal Diseases."
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
            summary="Variant does not meet the PM1 criteria for Lysosomal Diseases.",
        )
