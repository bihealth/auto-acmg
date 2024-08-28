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

from typing import Tuple

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

PM1_CLUSTER = {
    # RYR1
    "HGNC:10483": [(4800, 4950)],
}


class CongenitalMyopathiesPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override PM1 prediction for congenital myopathies."""
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

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGData) -> bool:
        """Override BP3 for congenital myopathies."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override PP2 and BP1 for RASopathy."""
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
