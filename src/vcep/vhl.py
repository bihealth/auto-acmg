"""
Predictor for VHL VCEP.
Included gene: VHL (HGNC:12687).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN078
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

# fmt: off
PM1_CLUSTER = {
    # "NM_000551.3": {
    "HGNC:12687": {
        "moderate": [
            # Stebbins er al (PMID: 10205047)
            167, 162, 178, 98, 78, 86,
            # Chiorean et al (PMID: 35475554)
            65, 76, 78, 80, 86, 88, 96, 98, 112, 117,
            161, 162, 167, 170, 176,
            # Walsh et al (PMIDs: 29247016, 30311369)
            68, 74, 89, 111, 114, 115, 121, 135, 151, 158, 169
        ]
    }
}


class VHLPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for VHL.
        """
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.transcript_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check if variant is in the moderate cluster
        if var_data.prot_pos in gene_cluster["moderate"]:
            comment = (
                f"Variant affects a germline hotspot or key functional domain in VHL "
                f"at position {var_data.prot_pos}. PM1 is met at the Moderate level."
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
            summary="Variant does not meet the PM1 criteria for VHL.",
        )
