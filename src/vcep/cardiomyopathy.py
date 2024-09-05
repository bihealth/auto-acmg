"""
Predictor for Cardiomyopathy VCEP.
Included genes:
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

from typing import List, Tuple

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

#: VCEP spcifications for Cardiomyopathy.
SPECs: List[VcepSpec] = [
    VcepSpec(
        identifier="GN002",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN095",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN098",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN099",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN100",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN101",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN102",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN103",
        version="1.0.0",
    ),
]

PM1_CLUSTER = {
    # MYH7
    # "NM_000257.4": [(167, 931)],
    # "ENST00000355349": [(167, 931)],
    "HGNC:7577": [(167, 931)],
    # MYBPC3
    # "NM_000256.3": [
    #     (485, 502),
    #     (1248, 1266),
    # ],
    # "ENST00000545968": [
    #     (485, 502),
    #     (1248, 1266),
    # ],
    "HGNC:7551": [
        (485, 502),
        (1248, 1266),
    ],
    # TNNI3
    # "NM_000363.5": [(141, 209)],
    # "ENST00000344887": [(141, 209)],
    "HGNC:11947": [(141, 209)],
    # TNNT2
    # "ENST00000367318": [(79, 179)],
    # "NM_001276345.2": [(89, 189)],
    # "ENST00000656932.1": [(89, 189)],
    "HGNC:11949": [(89, 189)],
}


class CardiomyopathyPredictor(DefaultSeqVarPredictor):

    def predict_pvs1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """PVS1 is not applicable."""
        if var_data.hgnc_id in [
            "HGNC:7577",
            "HGNC:11947",
            "HGNC:11949",
            "HGNC:12010",
            "HGNC:143",
            "HGNC:7583",
            "HGNC:7584",
        ]:
            logger.info("Predict PVS1")
            return AutoACMGCriteria(
                name="PVS1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicVeryStrong,
                summary="PVS1 is not applicable for the gene.",
            )
        return super().predict_pvs1(seqvar, var_data)

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Specify PM1 domains for Cardiomyopathy."""
        logger.info("Predict PM1")

        if var_data.hgnc_id in PM1_CLUSTER:
            domains = PM1_CLUSTER[var_data.hgnc_id]

            # Check if the variant falls within any of the specified domains
            for start_aa, end_aa in domains:
                if start_aa <= var_data.prot_pos <= end_aa:
                    comment = (
                        f"Variant falls within a critical domain for {var_data.hgnc_id} "
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

    def _bs2_not_applicable(self, var_data: AutoACMGSeqVarData) -> bool:
        """BS2 is not applicable for Cardiomyopathy."""
        return True

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.00004
        var_data.thresholds.ba1_benign = 0.001
        if var_data.hgnc_id == "HGNC:7551":
            var_data.thresholds.bs1_benign = 0.0002
        else:
            var_data.thresholds.bs1_benign = 0.0001
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 is not applicable for Cardiomyopathy."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override PP2 and BP1 for Cardiomyopathy to return not applicable status."""
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

    def predict_pp3bp4(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override PP3 and BP4 for Cardiomyopathy to use REVEL thresholds."""
        var_data.thresholds.pp3bp4_strategy = "revel"
        var_data.thresholds.revel_pathogenic = 0.7
        var_data.thresholds.revel_benign = 0.4
        return super().predict_pp3bp4(seqvar, var_data)

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override donor and acceptor positions for Cardiomyopathy VCEP."""
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 4
        return super().predict_bp7(seqvar, var_data)
