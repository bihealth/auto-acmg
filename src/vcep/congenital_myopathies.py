"""
Predictor for Congenital Myopathies VCEP.
Included genes:
NEB (HGNC:7720),
ACTA1 (HGNC:129),
DNM2 (HGNC:2974),
MTM1 (HGNC:7448),
RYR1 (HGNC:10483).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN146
https://cspec.genome.network/cspec/ui/svi/doc/GN147
https://cspec.genome.network/cspec/ui/svi/doc/GN148
https://cspec.genome.network/cspec/ui/svi/doc/GN149
https://cspec.genome.network/cspec/ui/svi/doc/GN150
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

#: VCEP specifications for Congenital Myopathies.
SPECs: List[VcepSpec] = [
    VcepSpec(
        identifier="GN146",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN147",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN148",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN149",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN150",
        version="1.0.0",
    ),
]


PM1_CLUSTER = {
    # RYR1
    "HGNC:10483": [(4800, 4950)],
}


class CongenitalMyopathiesPredictor(DefaultSeqVarPredictor):

    def predict_pvs1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """PVS1 is not applicable."""
        if var_data.hgnc_id == "HGNC:2974":
            logger.info("Predict PVS1")
            return AutoACMGCriteria(
                name="PVS1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicVeryStrong,
                summary="PVS1 is not applicable for the gene.",
            )
        return super().predict_pvs1(seqvar, var_data)

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override PM1 to specify critical domains for congenital myopathies."""
        logger.info("Predict PM1")

        # NEB, ACTA1, DNM2, MTM1
        if var_data.hgnc_id in ["HGNC:7720", "HGNC:129", "HGNC:2974", "HGNC:7448"]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicModerate,
                summary="Not applicable for NEB, ACTA1, DNM2, MTM1.",
            )

        if var_data.hgnc_id in PM1_CLUSTER:
            for start, end in PM1_CLUSTER[var_data.hgnc_id]:
                if start <= var_data.prot_pos <= end:
                    comment = (
                        f"Variant falls within a critical domain for {var_data.hgnc_id} "
                        f"between positions {start}-{end}. PM1 is met."
                    )
                    return AutoACMGCriteria(
                        name="PM1",
                        prediction=AutoACMGPrediction.Met,
                        strength=AutoACMGStrength.PathogenicModerate,
                        summary=comment,
                    )
        else:
            return super().predict_pm1(seqvar, var_data)

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary="Variant does not fall within a critical domain.",
        )

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.0000001  # Practically should be absent
        if var_data.hgnc_id == "HGNC:7720":
            var_data.thresholds.ba1_benign = 0.00559
            var_data.thresholds.bs1_benign = 0.000237
        elif var_data.hgnc_id == "HGNC:129":
            var_data.thresholds.ba1_benign = 0.0000781
            var_data.thresholds.bs1_benign = 0.00000781
        elif var_data.hgnc_id == "HGNC:2974":
            var_data.thresholds.ba1_benign = 0.0000015
            var_data.thresholds.bs1_benign = 0.00000015
        elif var_data.hgnc_id == "HGNC:7448":
            var_data.thresholds.ba1_benign = 0.000016
            var_data.thresholds.bs1_benign = 0.0000016
        elif var_data.hgnc_id == "HGNC:10483":
            var_data.thresholds.ba1_benign = 0.0000486
            var_data.thresholds.bs1_benign = 0.00000486
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 is not applicable for Congenital Myopathies."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override PP2, BP1 for Congenital Myopathies to return not applicable status."""
        pp2 = False
        comment = "Not applicable for the gene."
        if var_data.hgnc_id in ["HGNC:129", "HGNC:2974"]:
            if self._is_missense(var_data):
                pp2 = True
                comment = f"PP2 is met for {var_data.hgnc_id} as the variant is a missense change."
            else:
                pp2 = False
                comment = (
                    f"PP2 is not met for {var_data.hgnc_id} as the variant is not a missense "
                    "change."
                )
        return (
            AutoACMGCriteria(
                name="PP2",
                prediction=AutoACMGPrediction.Met if pp2 else AutoACMGPrediction.NotMet,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=comment,
            ),
            AutoACMGCriteria(
                name="BP1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.BenignSupporting,
                summary="BP1 is not applicable for the gene.",
            ),
        )

    def predict_pp3bp4(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Use REVEL scores for PP3 and BP4."""
        var_data.thresholds.pp3bp4_strategy = "revel"
        var_data.thresholds.revel_pathogenic = 0.7
        var_data.thresholds.revel_benign = 0.15
        var_data.thresholds.spliceAI_acceptor_gain = 0.5
        var_data.thresholds.spliceAI_acceptor_loss = 0.5
        var_data.thresholds.spliceAI_donor_gain = 0.5
        var_data.thresholds.spliceAI_donor_loss = 0.5
        return super().predict_pp3bp4(seqvar, var_data)
