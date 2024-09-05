"""
Predictor for Heriditary Breast, Ovarian and Pancreatic Cancer VCEP.
Included genes:
ATM (HGNC:795),
PALB2 (HGNC:26144).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN020
https://cspec.genome.network/cspec/ui/svi/doc/GN077
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
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException
from src.defs.seqvar import SeqVar
from src.seqvar.auto_ps1_pm5 import DNA_BASES
from src.seqvar.default_predictor import DefaultSeqVarPredictor

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


class HBOPCPredictor(DefaultSeqVarPredictor):

    def _is_nonsense(self, var_data: AutoACMGSeqVarData) -> bool:
        """
        Check if the variant is a nonsense/frameshift variant.

        Args:
            var_data (AutoACMGSeqVarData): The variant information.

        Returns:
            bool: True if the variant is a nonsense/frameshift variant.
        """
        # Check the nonsense
        if "nonsense" in var_data.consequence.cadd:
            return True
        if any("nonsense" in cons for cons in var_data.consequence.mehari):
            return True

        # Check the frameshift
        if "frameshift" in var_data.consequence.cadd:
            return True
        if any("frameshift" in cons for cons in var_data.consequence.mehari):
            return True

        # Stop gained
        if "stop_gained" in var_data.consequence.mehari:
            return True

        return False

    def verify_ps1pm5(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PS1PM5], str]:
        """Override PS1/PM5 for Hearing Loss."""
        self.prediction_ps1pm5, self.comment_ps1pm5 = super().verify_ps1pm5(seqvar, var_data)
        if var_data.hgnc_id == "HGNC:795" and self.prediction_ps1pm5:
            if self._is_missense(var_data):
                self.prediction_ps1pm5.PM5 = False
                if self._affect_splicing(var_data):
                    self.prediction_ps1pm5.PS1 = False
            elif self._is_splice_affecting(var_data):
                if (
                    self.prediction_ps1pm5.PM5
                    and var_data.prot_pos < 3047
                    and self.undergo_nmd(
                        var_data.tx_pos_utr, var_data.hgnc_id, var_data.strand, var_data.exons
                    )
                ):
                    self.comment_ps1pm5 += (
                        "Variant is splice-affecting variant and is upstream of p.3047 and "
                        "undergoes NMD. PM5 is met."
                    )
                    self.prediction_ps1pm5.PM5 = True
                else:
                    self.comment_ps1pm5 += (
                        "Variant is splice-affecting variant and doesn't fulfill PM5 or is "
                        "upstream of p.3047 or does not undergo NMD. PM5 is not met."
                    )
                    self.prediction_ps1pm5.PM5 = False
        elif var_data.hgnc_id == "HGNC:26144" and self.prediction_ps1pm5:
            if self._is_missense(var_data):
                self.prediction_ps1pm5.PS1 = False
                self.prediction_ps1pm5.PM5 = False
            elif self._is_splice_affecting(var_data):
                if (
                    self.prediction_ps1pm5.PM5
                    and var_data.prot_pos < 1183
                    and self.undergo_nmd(
                        var_data.tx_pos_utr, var_data.hgnc_id, var_data.strand, var_data.exons
                    )
                ):
                    self.comment_ps1pm5 += (
                        "Variant is splice-affecting variant and is upstream of p.1183 and "
                        "undergoes NMD. PM5 is met."
                    )
                    self.prediction_ps1pm5.PM5 = True
                else:
                    self.comment_ps1pm5 += (
                        "Variant is splice-affecting variant and doesn't fulfill PM5 or is "
                        "upstream of p.1183 or does not undergo NMD. PM5 is not met."
                    )
                    self.prediction_ps1pm5.PM5 = False

        # Check nonsense/frameshift variants for PM5
        try:
            if not self.prediction_ps1pm5 or not self._is_nonsense(var_data):
                raise AlgorithmError("Variant is not nonsense/frameshift variant.")

            if (
                var_data.hgnc_id == "HGNC:795"
                and var_data.prot_pos >= 3047
                or var_data.hgnc_id == "HGNC:26144"
                and var_data.prot_pos >= 1183
            ):
                raise AlgorithmError(
                    "Variant is nonsence/frameshift variant and is downstream of p.3047 or p.1183. "
                    "PM5 is not met."
                )

            primary_aa_change = self._parse_HGVSp(var_data.pHGVS)
            if not primary_aa_change:
                raise AlgorithmError("No valid primary amino acid change for PS1/PM5 prediction.")

            for alt_base in DNA_BASES:
                # Skip the same base insert
                if alt_base == seqvar.insert:
                    continue
                alt_seqvar = SeqVar(
                    genome_release=seqvar.genome_release,
                    chrom=seqvar.chrom,
                    pos=seqvar.pos,
                    delete=seqvar.delete,
                    insert=alt_base,
                )
                alt_info = self._get_var_info(alt_seqvar)

                if alt_info and alt_info.result.dbnsfp and alt_info.result.dbnsfp.HGVSp_VEP:
                    alt_aa_change = self._parse_HGVSp(alt_info.result.dbnsfp.HGVSp_VEP)
                    if alt_aa_change and self._is_pathogenic(alt_info.result):
                        if primary_aa_change != alt_aa_change:
                            self.comment_ps1pm5 += "PM5 is met as nonsense/frameshift variant."
                            self.prediction_ps1pm5.PM5 = True

        except AutoAcmgBaseException:
            pass
        return self.prediction_ps1pm5, self.comment_ps1pm5

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
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

    def _bs2_not_applicable(self, var_data: AutoACMGSeqVarData) -> bool:
        """BS2 is not applicable for ATM gene."""
        if var_data.hgnc_id == "HGNC:795":
            return True
        return False

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
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
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
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
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
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

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override donor and acceptor positions for ATM and PALB2."""
        if var_data.hgnc_id == "HGNC:26144":
            var_data.thresholds.bp7_donor = 7
            var_data.thresholds.bp7_acceptor = 21
        elif var_data.hgnc_id == "HGNC:795":
            var_data.thresholds.bp7_donor = 7
            var_data.thresholds.bp7_acceptor = 40
        return super().predict_bp7(seqvar, var_data)
