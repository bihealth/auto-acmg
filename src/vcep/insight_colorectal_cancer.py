"""
Predictor for InSIGHT Hereditary Colorectal Cancer/Polyposis VCEP.
Included genes:
APC (HGNC:583),
MLH1 (HGNC:7127),
MSH2 (HGNC:7325),
MSH6 (HGNC:7329),
PMS2 (HGNC:9122).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN089
https://cspec.genome.network/cspec/ui/svi/doc/GN115
https://cspec.genome.network/cspec/ui/svi/doc/GN137
https://cspec.genome.network/cspec/ui/svi/doc/GN138
https://cspec.genome.network/cspec/ui/svi/doc/GN139
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

#: VCEP specifications for InSIGHT Hereditary Colorectal Cancer/Polyposis.
SPECs: List[VcepSpec] = [
    VcepSpec(
        identifier="GN089",
        version="2.1.0",
    ),
    VcepSpec(
        identifier="GN115",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN137",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN138",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN139",
        version="1.0.0",
    ),
]


class InsightColorectalCancerPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override PM1 prediction for InSIGHT Hereditary Colorectal Cancer/Polyposis."""
        logger.info("Predict PM1")

        if var_data.hgnc_id in [
            "HGNC:583",  # APC
            "HGNC:7127",  # MLH1
            "HGNC:7325",  # MSH2
            "HGNC:7329",  # MSH6
            "HGNC:9122",  # PMS2
        ]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=f"PM1 is not applicable for {var_data.hgnc_id}.",
            )

        return super().predict_pm1(seqvar, var_data)

    def predict_pm4bp3(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        return (
            AutoACMGCriteria(
                name="PM4",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=(
                    "PM4 is not applicable for the InSIGHT Hereditary Colorectal Cancer/Polyposis "
                    "VCEP."
                ),
            ),
            AutoACMGCriteria(
                name="BP3",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.BenignSupporting,
                summary=(
                    "BP3 is not applicable for the InSIGHT Hereditary Colorectal Cancer/Polyposis "
                    "VCEP."
                ),
            ),
        )

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override PP2 and BP1 prediction for InSIGHT Hereditary Colorectal Cancer/Polyposis."""
        bp1 = False
        comment = "BP1 is not applicable for the gene."
        if var_data.hgnc_id == "HGNC:583":  # APC
            if self._is_missense(var_data):
                start_pos, end_pos = min(var_data.cds_start, var_data.cds_end), max(
                    var_data.cds_start, var_data.cds_end
                )
                _, benign_count, total_count = self._get_missense_vars(seqvar, start_pos, end_pos)
                benign_ratio = benign_count / total_count
                if benign_ratio > var_data.thresholds.pp2bp1_benign and not (
                    1021 <= seqvar.pos <= 1035
                ):
                    bp1 = True
                    comment = (
                        f"Benign ratio {benign_ratio} is greater than "
                        f"{var_data.thresholds.pp2bp1_benign} and not in β-catenin binding domain. "
                        "BP1 is met."
                    )
                else:
                    bp1 = False
                    comment = (
                        f"Benign ratio {benign_ratio} is less than "
                        f"{var_data.thresholds.pp2bp1_benign} or in β-catenin binding domain. "
                        "BP1 is not met."
                    )
            else:
                bp1 = False
                comment = "Variant is not a missense variant. BP1 is not met."
        return (
            AutoACMGCriteria(
                name="PP2",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=(
                    "PP2 is not applicable for the InSIGHT Hereditary Colorectal Cancer/Polyposis "
                    "VCEP."
                ),
            ),
            AutoACMGCriteria(
                name="BP1",
                prediction=AutoACMGPrediction.Met if bp1 else AutoACMGPrediction.NotMet,
                strength=AutoACMGStrength.BenignSupporting,
                summary=comment,
            ),
        )

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override donor and acceptor positions for InSIGHT Hereditary Colorectal Cancer/Polyposis
        VCEP.
        """
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
