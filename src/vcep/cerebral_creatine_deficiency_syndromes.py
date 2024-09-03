"""
Predictor for Cerebral Creatine Deficiency Syndromes VCEP.
Included genes:
GATM (HGNC:4175),
GAMT (HGNC:4136),
SLC6A8 (HGNC:11055).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN025
https://cspec.genome.network/cspec/ui/svi/doc/GN026
https://cspec.genome.network/cspec/ui/svi/doc/GN027
"""

from typing import List, Tuple

from src.defs.auto_acmg import (
    AlleleCondition,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.exceptions import MissingDataError
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

#: VCEP specifications for Cerebral Creatine Deficiency Syndromes.
SPECs: List[VcepSpec] = [
    VcepSpec(
        identifier="GN025",
        version="1.1.0",
    ),
    VcepSpec(
        identifier="GN026",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN027",
        version="1.1.0",
    ),
]


class CerebralCreatineDeficiencySyndromesPredictor(DefaultSeqVarPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override predict_pm1 to return a not applicable status for PM1."""
        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotApplicable,
            strength=AutoACMGStrength.PathogenicModerate,
            summary="PM1 is not applicable for Cerebral Creatine Deficiency Syndromes.",
        )

    def _check_zyg(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """
        Check the zygosity of the sequence variant.

        BS2 only to be used when variant is observed in the homozygous state in a healthy adult.

        Args:
            variant_data: The variant data.

        Returns:
            True if the variant is recessive (homozygous), dominant (heterozygous), or X-linked
            (hemizygous) disorder.
        """
        allele_condition = self._get_allele_cond(seqvar)
        self.comment_pm2ba1bs1bs2 += f"Allele condition: {allele_condition.name}.\n"
        controls_af = self._get_control_af(var_data)
        any_af = self._get_any_af(var_data)
        af = controls_af or any_af
        if not af or not af.bySex:
            self.comment_pm2ba1bs1bs2 += "No controls allele data found in control data.\n"
            raise MissingDataError("No raw data found in control data.")

        if not af.bySex.overall:
            self.comment_pm2ba1bs1bs2 += "No allele data found for overall in control data.\n"
            raise MissingDataError("No allele data found for overall in control data.")
        ac = af.bySex.overall.ac if af.bySex.overall.ac else 0
        nhomalt = af.bySex.overall.nhomalt if af.bySex.overall.nhomalt else 0
        self.comment_pm2ba1bs1bs2 += f"Allele count: {ac}, Nhomalt: {nhomalt}.\n"
        if allele_condition == AlleleCondition.Recessive:
            if nhomalt > 5:
                self.comment_pm2ba1bs1bs2 += (
                    f"Nhomalt {nhomalt} > 5.\n"
                    "The variant is in a recessive (homozygous) disorder."
                )
                return True
        return False

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        if var_data.hgnc_id == "HGNC:4175":
            var_data.thresholds.pm2_pathogenic = 0.000055
            var_data.thresholds.ba1_benign = 0.0005
            var_data.thresholds.bs1_benign = 0.0001
        elif var_data.hgnc_id == "HGNC:4136":
            var_data.thresholds.pm2_pathogenic = 0.0004
            var_data.thresholds.ba1_benign = 0.003
            var_data.thresholds.bs1_benign = 0.001
        elif var_data.hgnc_id == "HGNC:11055":
            var_data.thresholds.pm2_pathogenic = 0.00002
            var_data.thresholds.ba1_benign = 0.002
            var_data.thresholds.bs1_benign = 0.0002
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 is not applicable for Cerebral Creatine Deficiency Syndromes."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to return a not applicable status for PP2 and BP1."""
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

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override donor and acceptor positions for Cerebral Creatine Deficiency Syndromes VCEP."""
        if var_data.hgnc_id == "HGNC:4136":
            var_data.thresholds.bp7_donor = 7
            var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
