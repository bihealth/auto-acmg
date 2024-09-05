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

from typing import List, Optional, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    PS1PM5,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

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


class InsightColorectalCancerPredictor(DefaultSeqVarPredictor):

    def verify_ps1pm5(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PS1PM5], str]:
        """Override PS1/PM5 for InSIGHT Hereditary Colorectal Cancer/Polyposis."""
        self.prediction_ps1pm5, self.comment_ps1pm5 = super().verify_ps1pm5(seqvar, var_data)
        var_data.thresholds.spliceAI_acceptor_gain = 0.2
        var_data.thresholds.spliceAI_acceptor_loss = 0.2
        var_data.thresholds.spliceAI_donor_gain = 0.2
        var_data.thresholds.spliceAI_donor_loss = 0.2
        if (
            self.prediction_ps1pm5
            and self._is_missense(var_data)
            and self._affect_splicing(var_data)
        ):
            self.prediction_ps1pm5.PS1 = False
            self.comment_ps1pm5 = "Variant affects splicing. PS1 is not applicable."
            self.prediction_ps1pm5.PM5 = False
            self.comment_ps1pm5 = "Variant affects splicing. PM5 is not applicable."
        return self.prediction_ps1pm5, self.comment_ps1pm5

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for InSIGHT Hereditary Colorectal
        Cancer/Polyposis. Use default logic for all genes except APC, MLH1, MSH2, MSH6, and PMS2.
        For them return not applicable status for PM1.
        """
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

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.00002  # 1/50,000
        if var_data.hgnc_id == "HGNC:583":  # APC
            var_data.thresholds.pm2_pathogenic = 0.000003
        if var_data.hgnc_id == "HGNC:7329":  # MSH6
            var_data.thresholds.ba1_benign = 0.0022
            var_data.thresholds.bs1_benign = 0.00022
        elif var_data.hgnc_id == "HGNC:9122":  # PMS2
            var_data.thresholds.ba1_benign = 0.0028
            var_data.thresholds.bs1_benign = 0.0001
        else:
            var_data.thresholds.ba1_benign = 0.001
            var_data.thresholds.bs1_benign = 0.0001
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def predict_pm4bp3(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """PM4 and BP3 are not applicable for InSIGHT Hereditary Colorectal Cancer/Polyposis."""
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
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """
        Override PP2 and BP1 for InSIGHT Hereditary Colorectal Cancer/Polyposis VCEP. Check benign
        ratio for APC missense variants and exclude β-catenin binding domain. PP2 is not applicable
        for the gene.
        """
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

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """
        Override donor and acceptor positions for InSIGHT Hereditary Colorectal Cancer/Polyposis
        VCEP.
        """
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
