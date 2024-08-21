"""
Predictor for Cardiomyopathy VCEP.
Included gene:
MYH7 (HGNC:7577),
MYBPC3 (HGNC:7551),
TNNI3 (HGNC:11947),
TNNT2 (HGNC:11949),
TPM1 (HGNC:12010),
ACTC1 (HGNC:143),
MYL2 (HGNC:7583),
MYL3 (HGNC:7584).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN002
https://cspec.genome.network/cspec/ui/svi/doc/GN095
https://cspec.genome.network/cspec/ui/svi/doc/GN098
https://cspec.genome.network/cspec/ui/svi/doc/GN099
https://cspec.genome.network/cspec/ui/svi/doc/GN100
https://cspec.genome.network/cspec/ui/svi/doc/GN101
https://cspec.genome.network/cspec/ui/svi/doc/GN102
https://cspec.genome.network/cspec/ui/svi/doc/GN103
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

PM1_CLUSTER = {
    # MYH7
    "NM_000257.4": [(167, 931)],
    "ENST00000355349": [(167, 931)],
    # MYBPC3
    "NM_000256.3": [
        (485, 502),
        (1248, 1266),
    ],
    "ENST00000545968": [
        (485, 502),
        (1248, 1266),
    ],
    # TNNI3
    "NM_000363.5": [
        (141, 209),
    ],
    "ENST00000344887": [
        (141, 209),
    ],
    # TNNT2
    "ENST00000367318": [
        (79, 179),
    ],
    "NM_001276345.2": [
        (89, 189),
    ],
    "ENST00000656932.1": [
        (89, 189),
    ],
}


class CardiomyopathyPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override predict_pm1 to include VCEP-specific logic for Cardiomyopathy."""
        logger.info("Predict PM1")

        if var_data.transcript_id in PM1_CLUSTER:
            domains = PM1_CLUSTER[var_data.transcript_id]

            # Check if the variant falls within any of the specified domains
            for start_aa, end_aa in domains:
                if start_aa <= var_data.prot_pos <= end_aa:
                    comment = (
                        f"Variant falls within a critical domain for {var_data.transcript_id} "
                        f"between positions {start_aa}-{end_aa}. PM1 is met."
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
                summary=(
                    "Variant does not fall within any critical domain for the specified gene. "
                    "PM1 is not met."
                ),
            )

        # TPM1, ACTC1, MYL2, MYL3
        if var_data.hgnc_id in ["HGNC:12010", "HGNC:143", "HGNC:7583", "HGNC:7584"]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicModerate,
                summary="Not applicable for TPM1, ACTC1, MYL2, MYL3.",
            )

        return super().predict_pm1(seqvar, var_data)
