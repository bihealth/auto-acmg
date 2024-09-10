"""
Predictor for PTEN VCEP.
Included gene: PTEN (HGNC:9588).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN003
"""

from typing import Optional, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    PM4BP3,
    PP3BP4,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.exceptions import AutoAcmgBaseException
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

#: VCEP specification for PTEN.
SPEC: VcepSpec = VcepSpec(
    identifier="GN003",
    version="3.1.0",
)

PM1_CLUSTER = {
    "HGNC:9588": {  # PTEN
        "residues":
        # catalytic motifs
        list(range(90, 95)) + list(range(123, 131)) + list(range(166, 169)),
    }
}


class PTENPredictor(DefaultSeqVarPredictor):
    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override predict_pm1 to specify critical residues for PTEN."""
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check if the variant falls within the specified catalytic motifs
        if var_data.prot_pos in gene_cluster.get("residues", []):
            comment = (
                f"Variant affects a critical residue within the catalytic motifs of PTEN "
                f"at position {var_data.prot_pos}. PM1 is met at the Moderate level."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Applicable,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=comment,
            )

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotApplicable,
            strength=AutoACMGStrength.PathogenicModerate,
            summary="Variant does not meet the PM1 criteria for PTEN.",
        )

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.00001
        var_data.thresholds.ba1_benign = 0.00056
        var_data.thresholds.bs1_benign = 0.000043
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 is not applicable for PTEN."""
        return True

    def verify_pm4bp3(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PM4BP3], str]:
        """
        Override predict_pm4bp3 to include VCEP-specific logic for PTEN. Include in-frame deletions
        and insertions that affect the catalytic motifs as PM4. Stop-loss variants are also
        considered as PM4.
        """
        self.prediction_pm4bp3 = PM4BP3()
        self.comment_pm4bp3 = ""
        try:
            # Stop-loss variants are considered as PM4
            if self._is_stop_loss(var_data):
                self.comment_pm4bp3 = "Variant consequence is stop-loss. PM4 is met."
                self.prediction_pm4bp3.PM4 = True
            # In-frame deletions/insertions
            elif self.is_inframe_delins(var_data):
                self.comment_pm4bp3 = "Variant consequence is in-frame deletion/insertion. "

                # Check if the variant affects the catalytic motifs
                if var_data.prot_pos in PM1_CLUSTER.get(var_data.hgnc_id, {}).get("residues", []):
                    self.comment_pm4bp3 += "Impacting catalytic motif. PM4 is met."
                    self.prediction_pm4bp3.PM4 = True
                else:
                    self.comment_pm4bp3 += "No impact on catalytic motif. PM4 is not met."
                    self.prediction_pm4bp3.PM4 = False
            else:
                self.comment_pm4bp3 = (
                    "Variant consequence is not stop-loss or in-frame deletion/insertion. "
                )
                self.prediction_pm4bp3.PM4 = False
                self.prediction_pm4bp3.BP3 = False
        except Exception as e:
            self.prediction_pm4bp3 = None
            self.comment_pm4bp3 = f"An error occured while predicting PM4 and BP3 criteria: {e}"

        return self.prediction_pm4bp3, self.comment_pm4bp3

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to return not applicable status for PTEN. PP2 is not changed."""
        logger.info("Predict PP2 and BP1")
        pred, comment = self.verify_pp2bp1(seqvar, var_data)
        if pred:
            pp2_pred = (
                AutoACMGPrediction.Applicable
                if pred.PP2
                else (
                    AutoACMGPrediction.NotApplicable
                    if pred.PP2 is False
                    else AutoACMGPrediction.Failed
                )
            )
            pp2_strength = pred.PP2_strength
        else:
            pp2_pred = AutoACMGPrediction.Failed
            pp2_strength = AutoACMGStrength.PathogenicSupporting
        return (
            AutoACMGCriteria(
                name="PP2",
                prediction=pp2_pred,
                strength=pp2_strength,
                summary=comment,
            ),
            AutoACMGCriteria(
                name="BP1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary="BP1 is not applicable for the gene.",
            ),
        )

    def verify_pp3bp4(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PP3BP4], str]:
        """Predict PP3 and BP4 criteria."""
        self.prediction_pp3bp4 = PP3BP4()
        self.comment_pp3bp4 = ""
        try:
            score = "revel"
            var_data.thresholds.revel_pathogenic = 0.7
            var_data.thresholds.revel_benign = 0.5
            self.prediction_pp3bp4.PP3 = self._is_pathogenic_score(
                var_data,
                (score, getattr(var_data.thresholds, f"{score}_pathogenic")),
            )
            self.prediction_pp3bp4.BP4 = self._is_benign_score(
                var_data,
                (score, getattr(var_data.thresholds, f"{score}_benign")),
            )

            var_data.thresholds.spliceAI_acceptor_gain = 0.5
            var_data.thresholds.spliceAI_acceptor_loss = 0.5
            var_data.thresholds.spliceAI_donor_gain = 0.5
            var_data.thresholds.spliceAI_donor_loss = 0.5
            self.prediction_pp3bp4.PP3 = self.prediction_pp3bp4.PP3 or self._affect_spliceAI(
                var_data
            )

            if self._is_synonymous_variant(var_data) or self._is_intronic(var_data):
                var_data.thresholds.spliceAI_acceptor_gain = 0.2
                var_data.thresholds.spliceAI_acceptor_loss = 0.2
                var_data.thresholds.spliceAI_donor_gain = 0.2
                var_data.thresholds.spliceAI_donor_loss = 0.2
                self.prediction_pp3bp4.BP4 = not self._affect_spliceAI(var_data)

        except AutoAcmgBaseException as e:
            self.comment_pp3bp4 = f"An error occurred during prediction. Error: {e}"
            self.prediction_pp3bp4 = None
        return self.prediction_pp3bp4, self.comment_pp3bp4

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override donor and acceptor positions for PTEN VCEP."""
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
