"""
Predictor for VHL VCEP.
Included gene: VHL (HGNC:12687).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN078
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

#: VCEP specification for VHL.
SPEC: VcepSpec = VcepSpec(
    identifier="GN078",
    version="1.0.0",
)

# fmt: off
PM1_CLUSTER = {
    # "NM_000551.3": {
    "HGNC:12687": {
        "residues": [
            # Stebbins er al (PMID: 10205047)
            167, 162, 178, 98, 78, 86,
            # Chiorean et al (PMID: 35475554)
            65, 76, 78, 80, 86, 88, 96, 98, 112, 117,
            161, 162, 167, 170, 176,
            # Walsh et al (PMIDs: 29247016, 30311369)
            68, 74, 89, 111, 114, 115, 121, 135, 151, 158, 169
        ]
    }
}

#: Important domains for PM4 in VHL
PM4_CLUSTER = [
    (63, 155),  # Beta (β) domain
    (156, 192), # Alpha (ɑ) domain
    (193, 204), # Second Beta (β) domain
]

#: Repeat regions for BP3 in VHL
BP3_REPEAT_REGIONS = [
    (14, 48)
]


class VHLPredictor(DefaultSeqVarPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override PM1 prediction to specify critical residues for VHL."""
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check if variant is in the moderate cluster
        if var_data.prot_pos in gene_cluster["residues"]:
            comment = (
                f"Variant affects a germline hotspot or key functional domain in VHL "
                f"at position {var_data.prot_pos}. PM1 is met at the Moderate level."
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
            summary="Variant does not meet the PM1 criteria for VHL.",
        )

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.00000156
        var_data.thresholds.ba1_benign = 0.000156
        var_data.thresholds.bs1_benign = 0.0000156
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    @staticmethod
    def _in_vhl_important_domain(var_data: AutoACMGSeqVarData) -> bool:
        """Check if the variant is in an important VHL domain."""
        for start, end in PM4_CLUSTER:
            if start <= var_data.prot_pos <= end:
                return True
        return False

    @staticmethod
    def _in_gxeex_repeat_region(var_data: AutoACMGSeqVarData) -> bool:
        """Check if the variant is in the GXEEX repeat region."""
        for start, end in BP3_REPEAT_REGIONS:
            if start <= var_data.prot_pos <= end:
                return True
        return False

    def verify_pm4bp3(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> Tuple[Optional[PM4BP3], str]:
        """
        Override PM4 and BP3 verification for VHL. PM4 is met if the variant is in an important
        domain. BP3 is met if the variant is in the GXEEX repeat region.
        """
        self.prediction_pm4bp3 = PM4BP3()
        self.comment_pm4bp3 = ""
        try:
            # Stop-loss variants are considered as PM4
            if self._is_stop_loss(var_data):
                self.comment_pm4bp3 = "Variant consequence is stop-loss. PM4 is met."
                self.prediction_pm4bp3.PM4 = True
                self.prediction_pm4bp3.BP3 = False
            # In-frame deletions/insertions
            elif self.is_inframe_delins(var_data):
                self.comment_pm4bp3 = "Variant consequence is in-frame deletion/insertion. "
                if self._in_vhl_important_domain(var_data):
                    self.comment_pm4bp3 += (
                        "Variant is in an important domain of VHL. PM4 is met."
                    )
                    self.prediction_pm4bp3.PM4 = True
                    self.prediction_pm4bp3.BP3 = False
                elif self._in_gxeex_repeat_region(var_data):
                    self.comment_pm4bp3 += (
                        "Variant is in the GXEEX repeat motif in the VHL gene. BP3 is met."
                    )
                    self.prediction_pm4bp3.PM4 = False
                    self.prediction_pm4bp3.BP3 = True
                else:
                    self.comment_pm4bp3 += (
                        "Variant is not in an important domain or in a repeat region. BP3 is met."
                    )
                    self.prediction_pm4bp3.PM4 = False
                    self.prediction_pm4bp3.BP3 = True
            else:
                self.comment_pm4bp3 = (
                    "Variant consequence is not an in-frame indel or stop-loss. PM4 and BP3 are not met."
                )
                self.prediction_pm4bp3.PM4 = False
                self.prediction_pm4bp3.BP3 = False

        except AutoAcmgBaseException as e:
            logger.error("Failed to predict PM4 and BP3 criteria. Error: {}", e)
            self.comment_pm4bp3 = f"An error occurred while predicting PM4 and BP3 criteria: {e}"
            self.prediction_pm4bp3 = None
        return self.prediction_pm4bp3, self.comment_pm4bp3

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to return not applicable status for VHL."""
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

    def verify_pp3bp4(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PP3BP4], str]:
        """Predict PP3 and BP4 criteria."""
        self.prediction_pp3bp4 = PP3BP4()
        self.comment_pp3bp4 = ""
        try:
            score = "revel"
            var_data.thresholds.revel_pathogenic = 0.664
            var_data.thresholds.revel_benign = 0.3
            self.prediction_pp3bp4.PP3 = self._is_pathogenic_score(
                var_data,
                (score, getattr(var_data.thresholds, f"{score}_pathogenic")),
            )

            var_data.thresholds.spliceAI_acceptor_gain = 0.5
            var_data.thresholds.spliceAI_acceptor_loss = 0.5
            var_data.thresholds.spliceAI_donor_gain = 0.5
            var_data.thresholds.spliceAI_donor_loss = 0.5
            self.prediction_pp3bp4.PP3 = self.prediction_pp3bp4.PP3 or self._affect_spliceAI(
                var_data
            )
            var_data.thresholds.spliceAI_acceptor_gain = 0.1
            var_data.thresholds.spliceAI_acceptor_loss = 0.1
            var_data.thresholds.spliceAI_donor_gain = 0.1
            var_data.thresholds.spliceAI_donor_loss = 0.1
            self.prediction_pp3bp4.BP4 = self._affect_spliceAI(var_data)

        except AutoAcmgBaseException as e:
            self.comment_pp3bp4 = f"An error occurred during prediction. Error: {e}"
            self.prediction_pp3bp4 = None
        return self.prediction_pp3bp4, self.comment_pp3bp4


    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Change the BP7 threshold for PhyloP."""
        var_data.thresholds.phyloP100 = 0.2
        return super().predict_bp7(seqvar, var_data)
