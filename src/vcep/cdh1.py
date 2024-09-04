"""
Predictor for CDH1 VCEP.
Included gene: CDH1 (HGNC:1748).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN007
"""

from typing import Tuple

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

#: VCEP specification for CDH1.
SPEC: VcepSpec = VcepSpec(
    identifier="GN007",
    version="3.1.0",
)


class CDH1Predictor(DefaultSeqVarPredictor):

    def predict_ps1pm5(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_ps1pm5 to return a not applicable status for PS1 and evaluate PM5."""
        ps1 = AutoACMGCriteria(
            name="PS1",
            prediction=AutoACMGPrediction.NotApplicable,
            strength=AutoACMGStrength.PathogenicSupporting,
            summary="PS1 is not applicable for CDH1.",
        )
        if any(
            cons in ["frameshift_variant", "stop_gained"] for cons in var_data.consequence.mehari
        ) and self.undergo_nmd(
            var_data.tx_pos_utr, var_data.hgnc_id, var_data.strand, var_data.exons
        ):
            pm5 = AutoACMGCriteria(
                name="PM5",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary="Nonsense or frameshift variant predicted to undergo NMD. PM5 is met.",
            )
        else:
            pm5 = AutoACMGCriteria(
                name="PM5",
                prediction=AutoACMGPrediction.NotMet,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=(
                    "Consequence is not frameshift or nonsense or variant is not predicted to "
                    "undergo NMD. PM5 is not met."
                ),
            )
        return ps1, pm5

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override predict_pm1 to return a not applicable status for PM1."""
        logger.info("Predict PM1")

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotApplicable,
            strength=AutoACMGStrength.PathogenicSupporting,
            summary="PM1 is not applicable for CDH1.",
        )

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.00001
        var_data.thresholds.ba1_benign = 0.002
        var_data.thresholds.bs1_benign = 0.001
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def predict_pm4bp3(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pm4bp3 for CDH1 VCEP. PM4 is not changed, but BP3 is not applicable."""
        logger.info("Predict PM4 and BP3")
        pred, comment = super().verify_pm4bp3(seqvar, var_data)
        if pred:
            pm4 = (
                AutoACMGPrediction.Met
                if pred.PM4
                else (AutoACMGPrediction.NotMet if pred.PM4 is False else AutoACMGPrediction.Failed)
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
                summary="BP3 is not applicable for CDH1.",
            ),
        )

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
        """Override donor and acceptor positions for CDH1 VCEP."""
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
