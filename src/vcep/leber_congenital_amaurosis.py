"""
Predictor for Leber Congenital Amaurosis/early onset Retinal Dystrophy VCEP.
Included gene: RPE65 (HGNC:10294)
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN120
"""

from typing import Optional, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    PM2BA1BS1BS2,
    PP3BP4,
    PS1PM5,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    GenomicStrand,
    VcepSpec,
)
from src.defs.exceptions import AutoAcmgBaseException
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

#: VCEP specifications for Leber Congenital Amaurosis/early onset Retinal Dystrophy.
SPEC: VcepSpec = VcepSpec(
    identifier="GN120",
    version="1.0.0",
)

PM1_CLUSTER = {
    "HGNC:10294": {
        "residues": [180, 182, 241, 313, 417, 527] + list(range(107, 126)),
    }
}


class LeberCongenitalAmaurosisPredictor(DefaultSeqVarPredictor):

    def verify_ps1pm5(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PS1PM5], str]:
        """Override PS1/PM5 for Leber Congenital Amaurosis/early onset Retinal Dystrophy."""
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
        """Override predict_pm1 to include domain information for RPE65."""
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        if var_data.prot_pos in gene_cluster["residues"]:
            comment = (
                f"Variant affects a residue in RPE65 at position {var_data.prot_pos}, "
                f"which is associated with Leber Congenital Amaurosis/early onset Retinal "
                "Dystrophy."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=comment,
            )

        # If no criteria are met, return Not Met
        comment = (
            f"Variant does not meet the PM1 criteria for Leber Congenital Amaurosis/early onset "
            "Retinal Dystrophy."
        )
        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary=comment,
        )

    def verify_pm2ba1bs1bs2(
        self,
        seqvar: SeqVar,
        var_data: AutoACMGSeqVarData,
    ) -> Tuple[Optional[PM2BA1BS1BS2], str]:
        """
        Predicts the PM2, BA1, BS1, BS2 criteria for the sequence variant.

        Note:
            Rules:
            PM2: Absent from controls allele frequency data.

            BA1: Allele frequency is greater than 0,008.

            BS1: Allele frequency is between 0.0008 and 0.008.

            BS2: No change from the reference implementation.

        Returns:
            BA1BS1BS2PM2: The prediction result.
        """
        self.prediction_pm2ba1bs1bs2 = PM2BA1BS1BS2()
        self.comment_pm2ba1bs1bs2 = ""
        var_data.thresholds.ba1_benign = 0.008
        var_data.thresholds.bs1_benign = 0.0008
        try:
            af = self._get_af(seqvar, var_data)
            if not af:
                self.comment_pm2ba1bs1bs2 = "No allele frequency data found. "
            elif self._ba1_exception(seqvar):
                self.comment_pm2ba1bs1bs2 = "The variant is in the exception list for BA1 criteria."
                self.prediction_pm2ba1bs1bs2.BA1 = False
                self.prediction_pm2ba1bs1bs2.BS1 = False
            elif af >= var_data.thresholds.ba1_benign:
                self.comment_pm2ba1bs1bs2 = "Allele frequency > 5%: BA1 is met. "
                self.prediction_pm2ba1bs1bs2.BA1 = True
            elif af >= var_data.thresholds.bs1_benign:
                self.comment_pm2ba1bs1bs2 = "Allele frequency > 1%: BS1 is met. "
                self.prediction_pm2ba1bs1bs2.BS1 = True
            elif not self._get_control_af(var_data):
                self.comment_pm2ba1bs1bs2 = "No control allele frequency data found. "
                self.prediction_pm2ba1bs1bs2.PM2 = True

            if not self.prediction_pm2ba1bs1bs2.BA1 and af and self._check_zyg(seqvar, var_data):
                self.comment_pm2ba1bs1bs2 += (
                    "The variant is in a recessive, dominant, or X-linked disorder: BS2 is met."
                )
                self.prediction_pm2ba1bs1bs2.BS2 = True

        except AutoAcmgBaseException as e:
            logger.error("Error occurred during PM2, BA1, BS1, BS2 prediction. Error: {}", e)
            self.comment_pm2ba1bs1bs2 = (
                f"An error occurred while predicting PM2, BA1, BS1, BS2 criteria: {e}"
            )
            self.prediction_pm2ba1bs1bs2 = None

        return self.prediction_pm2ba1bs1bs2, self.comment_pm2ba1bs1bs2

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 is not applicable for RPE65."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to return not applicable status for RPE65."""
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
            var_data.thresholds.revel_pathogenic = 0.644
            var_data.thresholds.revel_benign = 0.29
            self.prediction_pp3bp4.PP3 = self._is_pathogenic_score(
                var_data,
                (score, getattr(var_data.thresholds, f"{score}_pathogenic")),
            )
            self.prediction_pp3bp4.BP4 = self._is_benign_score(
                var_data,
                (score, getattr(var_data.thresholds, f"{score}_benign")),
            )

            var_data.thresholds.spliceAI_acceptor_gain = 0.2
            var_data.thresholds.spliceAI_acceptor_loss = 0.2
            var_data.thresholds.spliceAI_donor_gain = 0.2
            var_data.thresholds.spliceAI_donor_loss = 0.2
            self.prediction_pp3bp4.PP3 = self.prediction_pp3bp4.PP3 or self._affect_spliceAI(
                var_data
            )
            var_data.thresholds.spliceAI_acceptor_gain = 0.1
            var_data.thresholds.spliceAI_acceptor_loss = 0.1
            var_data.thresholds.spliceAI_donor_gain = 0.1
            var_data.thresholds.spliceAI_donor_loss = 0.1
            self.prediction_pp3bp4.BP4 = self.prediction_pp3bp4.BP4 and not self._affect_spliceAI(
                var_data
            )

        except AutoAcmgBaseException as e:
            self.comment_pp3bp4 = f"An error occurred during prediction. Error: {e}"
            self.prediction_pp3bp4 = None
        return self.prediction_pp3bp4, self.comment_pp3bp4

    def _is_bp7_exception(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """
        Add an exception for RPE65.

        Positions excluded:
        Synonymous substitutions at the first base of an exon
        Synonymous substitutions in the last 3 bases of an exon
        """
        for exon in var_data.exons:
            if var_data.strand == GenomicStrand.Plus:
                if exon.altStartI <= seqvar.pos <= exon.altStartI + 1:
                    return True
                if exon.altEndI - 3 <= seqvar.pos <= exon.altEndI:
                    return True
            elif var_data.strand == GenomicStrand.Minus:
                if exon.altStartI <= seqvar.pos <= exon.altStartI + 3:
                    return True
                if exon.altEndI - 1 <= seqvar.pos <= exon.altEndI:
                    return True
        return False

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override donor and acceptor positions for RPE65 gene."""
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
