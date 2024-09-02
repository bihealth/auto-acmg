"""
Predictor for Heriditary Breast, Ovarian and Pancreatic Cancer VCEP.
Included genes:
ATM (HGNC:795),
PALB2 (HGNC:26144).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN020
https://cspec.genome.network/cspec/ui/svi/doc/GN077
"""

from typing import List, Tuple

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.seqvar import SeqVar

#: VCEP specifications for Heriditary Breast, Ovarian and Pancreatic Cancer.
SPECs: List[VcepSpec] = [
    VcepSpec(
        identifier="GN020",
        version="1.3.0",
    ),
    VcepSpec(
        identifier="GN077",
        version="1.1.0",
    ),
]


class HBOPCPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Heriditary Breast, Ovarian and
        Pancreatic Cancer. ATM and PALB2 are not applicable for PM1.
        """
        logger.info("Predict PM1")

        if var_data.hgnc_id in [
            "HGNC:795",  # ATM
            "HGNC:26144",  # PALB2
        ]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=f"PM1 is not applicable for {var_data.hgnc_id}.",
            )

        return super().predict_pm1(seqvar, var_data)

    def _bs2_not_applicable(self, var_data: AutoACMGData) -> bool:
        """BS2 is not applicable for ATM gene."""
        if var_data.hgnc_id == "HGNC:795":
            return True
        return False

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        if var_data.hgnc_id == "HGNC:795":
            var_data.thresholds.pm2_pathogenic = 0.00001
            var_data.thresholds.ba1_benign = 0.005
            var_data.thresholds.bs1_benign = 0.0005
        elif var_data.hgnc_id == "HGNC:26144":
            var_data.thresholds.pm2_pathogenic = 0.000003  # 1/300,000
            var_data.thresholds.ba1_benign = 0.001
            var_data.thresholds.bs1_benign = 0.0001
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def predict_pm4bp3(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """
        Override predict_pm4bp3 to include VCEP-specific logic for CDH1. PM4 is not changed, but
        BP3 is not applicable.
        """
        logger.info("Predict PM4 and BP3")

        if var_data.hgnc_id == "HGNC:795":
            pred, comment = self.verify_pm4bp3(seqvar, var_data)
            if pred:
                pm4 = (
                    AutoACMGPrediction.Met
                    if pred.PM4
                    else (
                        AutoACMGPrediction.NotMet
                        if pred.PM4 is False
                        else AutoACMGPrediction.Failed
                    )
                )
            else:
                pm4 = AutoACMGPrediction.Failed
                comment = "PM4 could not be verified."
            return (
                AutoACMGCriteria(
                    name="PM4",
                    prediction=pm4,
                    strength=AutoACMGStrength.PathogenicSupporting,
                    summary=comment,
                ),
                AutoACMGCriteria(
                    name="BP3",
                    prediction=AutoACMGPrediction.NotApplicable,
                    strength=AutoACMGStrength.BenignSupporting,
                    summary="BP3 is not applicable for ATM.",
                ),
            )
        elif var_data.hgnc_id == "HGNC:26144":
            return (
                AutoACMGCriteria(
                    name="PM4",
                    prediction=AutoACMGPrediction.NotApplicable,
                    strength=AutoACMGStrength.PathogenicSupporting,
                    summary="PM4 is not applicable for PALB2.",
                ),
                AutoACMGCriteria(
                    name="BP3",
                    prediction=AutoACMGPrediction.NotApplicable,
                    strength=AutoACMGStrength.BenignSupporting,
                    summary="BP3 is not applicable for PALB2.",
                ),
            )
        return super().predict_pm4bp3(seqvar, var_data)

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """
        Override predict_pp2bp1 to include VCEP-specific logic for ATM and PALB2. Check if the
        variant is synonymous, missense or inframe indel and not in an important domain and not
        predicted to affect splicing. PP2 is not applicable.
        """
        logger.info("Predict PP2 and BP1")
        if var_data.hgnc_id == "HGNC:26144":
            if (
                "missense_variant" in var_data.consequence.mehari
                or var_data.consequence.cadd == "missense"
            ):
                bp1 = True
                comment = "Variant is missense. BP1 is met."
            else:
                bp1 = False
                comment = "Variant is not missense. BP1 is not met."
        else:
            bp1 = False
            comment = "BP1 is not applicable for the gene."
        return (
            AutoACMGCriteria(
                name="PP2",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary="PP2 is not applicable for the gene.",
            ),
            AutoACMGCriteria(
                name="BP1",
                prediction=AutoACMGPrediction.Met if bp1 else AutoACMGPrediction.NotMet,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=comment,
            ),
        )

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override donor and acceptor positions for ATM and PALB2."""
        if var_data.hgnc_id == "HGNC:26144":
            var_data.thresholds.bp7_donor = 7
            var_data.thresholds.bp7_acceptor = 21
        elif var_data.hgnc_id == "HGNC:795":
            var_data.thresholds.bp7_donor = 7
            var_data.thresholds.bp7_acceptor = 40
        return super().predict_bp7(seqvar, var_data)
