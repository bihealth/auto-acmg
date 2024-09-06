"""
Predictor for ENIGMA BRCA1 and BRCA2 VCEP.
Included genes:
BRCA1 (HGNC:1100),
BRCA2 (HGNC:1101).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN092
https://cspec.genome.network/cspec/ui/svi/doc/GN097
"""

from typing import List, Optional, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    BP7,
    PM2BA1BS1BS2,
    PP3BP4,
    PS1PM5,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.exceptions import AutoAcmgBaseException
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

#: VCEP specifications for ENIGMA BRCA1 and BRCA2.
SPECs: List[VcepSpec] = [
    VcepSpec(
        identifier="GN092",
        version="1.1.0",
    ),
    VcepSpec(
        identifier="GN097",
        version="1.1.0",
    ),
]

BP7_IMPORTANT_DOMAINS = {
    "HGNC:1100": [
        (2, 102),  # RING domain
        (1391, 1423),  # coiled-coil aa
        (1650, 1857),  # BRCT domain
    ],
    "HGNC:1101": [
        (10, 40),  # PALB2 binding
        (2481, 3186),  # DNA binding
    ],
}


class ENIGMAPredictor(DefaultSeqVarPredictor):

    def verify_ps1pm5(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PS1PM5], str]:
        """Override PS1/PM5 for ENIGMA BRCA1 and BRCA2."""
        self.prediction_ps1pm5, self.comment_ps1pm5 = super().verify_ps1pm5(seqvar, var_data)
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
        """Override PM1 to return a not applicable status."""
        logger.info("Predict PM1")

        if var_data.hgnc_id in ["HGNC:1100", "HGNC:1101"]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicModerate,
                summary="PM1 is not applicable for BRCA1 and BRCA2.",
            )

        return super().predict_pm1(seqvar, var_data)

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

            BA1: Allele frequency is greater than 0.0926%.

            BS1: Allele frequency is between 0.0185% and 0.0926%.

            BS2: No change from the reference implementation.

        Returns:
            BA1BS1BS2PM2: The prediction result.
        """
        self.prediction_pm2ba1bs1bs2 = PM2BA1BS1BS2()
        self.comment_pm2ba1bs1bs2 = ""
        var_data.thresholds.ba1_benign = 0.001
        var_data.thresholds.bs1_benign = 0.00002
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

    def predict_pm4bp3(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pm4bp3 to return not applicable status for PM4 and BP3."""
        logger.info("Predict PM4 and BP3")
        return (
            AutoACMGCriteria(
                name="PM4",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicModerate,
                summary="PM4 is not applicable for the ENIGMA VCEP.",
            ),
            AutoACMGCriteria(
                name="BP3",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.BenignSupporting,
                summary="BP3 is not applicable for the ENIGMA VCEP.",
            ),
        )

    def _in_important_domain(self, var_data: AutoACMGSeqVarData) -> bool:
        """Check if the variant is in an important domain."""
        for start, end in BP7_IMPORTANT_DOMAINS.get(var_data.hgnc_id, []):
            if start <= var_data.prot_pos <= end:
                return True
        return False

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """
        Override predict_pp2bp1 to include VCEP-specific logic for ENIGMA BRCA1 and BRCA2. Check if
        the variant is synonymous, missense or inframe indel and not in an important domain and not
        predicted to affect splicing.
        """
        if (
            any(
                crit in var_data.consequence.mehari or var_data.consequence.cadd == crit
                for crit in [
                    "synonymous_variant",
                    "missense_variant",
                    "inframe_insertion",
                    "inframe_deletion",
                    "synonymous",
                    "missense",
                    "inframe",
                ]
            )
            and not self._in_important_domain(var_data)
            and not self._spliceai_impact(var_data)
        ):
            bp1 = True
            comment = (
                "Variant is synonymous, missense or inframe indel and not in an important domain "
                "and not predicted to affect splicing. BP1 is met."
            )
        else:
            bp1 = False
            comment = (
                "Variant is not synonymous, missense or inframe indel, or in an important domain, "
                "or predicted to affect splicing. BP1 is not met."
            )
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
                strength=AutoACMGStrength.PathogenicStrong,
                summary=comment,
            ),
        )

    def predict_pp3bp4(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Predict PP3 and BP4 criteria based on ENIGMA VCEP specific rules."""
        logger.info("Predict PP3 and BP4")
        pp3_met = False
        bp4_met = False
        comments = []
        self.prediction_pp3bp4 = PP3BP4()

        # Evaluate missense changes
        if self._is_missense_variant(var_data) or self._is_inframe_indel(var_data):
            bayesDel_noAF_score = var_data.scores.dbnsfp.bayesDel_noAF
            var_data.thresholds.bayesDel_noAF_pathogenic = 0.28
            var_data.thresholds.bayesDel_noAF_benign = 0.15
            if self._in_important_domain(var_data) and self._is_pathogenic_score(
                var_data, ("bayesDel_noAF", var_data.thresholds.bayesDel_noAF_pathogenic)
            ):
                pp3_met = True
                comments.append(
                    f"BayesDel_noAF score {bayesDel_noAF_score} > {var_data.thresholds.bayesDel_noAF_pathogenic}, PP3 met."
                )
            if (
                self._in_important_domain(var_data)
                and self._is_benign_score(
                    var_data, ("bayesDel_noAF", var_data.thresholds.bayesDel_noAF_benign)
                )
                and not self._affect_spliceAI(var_data)
            ):
                bp4_met = True
                comments.append(
                    f"BayesDel_noAF score {bayesDel_noAF_score} < {var_data.thresholds.bayesDel_noAF_benign}, BP4 met."
                )

        # Evaluate splice changes for pp3
        if (
            self._is_missense_variant(var_data)
            or self._is_inframe_indel(var_data)
            or self._is_intronic(var_data)
            or self._is_synonymous_variant(var_data)
        ):
            var_data.thresholds.spliceAI_acceptor_gain = 0.2
            var_data.thresholds.spliceAI_acceptor_loss = 0.2
            var_data.thresholds.spliceAI_donor_gain = 0.2
            var_data.thresholds.spliceAI_donor_loss = 0.2
            if self._affect_spliceAI(var_data):
                pp3_met = True
                comments.append("SpliceAI ≥0.2, PP3 met.")

        # Evaluate splice changes for bp4
        if self._is_intronic(var_data) or self._is_synonymous_variant(var_data):
            var_data.thresholds.spliceAI_acceptor_gain = 0.1
            var_data.thresholds.spliceAI_acceptor_loss = 0.1
            var_data.thresholds.spliceAI_donor_gain = 0.1
            var_data.thresholds.spliceAI_donor_loss = 0.1
            if not self._affect_spliceAI(var_data):
                bp4_met = True
                comments.append("SpliceAI ≥0.2, BP4 met.")

        # Set criteria results
        pp3_result = AutoACMGCriteria(
            name="PP3",
            prediction=AutoACMGPrediction.Met if pp3_met else AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicSupporting,
            summary=" | ".join(comments) if pp3_met else "PP3 criteria not met.",
        )
        bp4_result = AutoACMGCriteria(
            name="BP4",
            prediction=AutoACMGPrediction.Met if bp4_met else AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.BenignSupporting,
            summary=" | ".join(comments) if bp4_met else "BP4 criteria not met.",
        )

        return (pp3_result, bp4_result)

    def verify_bp7(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> Tuple[Optional[BP7], str]:
        """
        Override verify BP7 criterion for ENIGMA BRCA1 and BRCA2. Check if the variant is synonymous
        and in an important domain or intronic and not conserved.
        """
        # Change the donor/acceptor thresholds for BRCA1 and BRCA2
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21

        self.prediction_bp7 = BP7()
        self.comment_bp7 = ""
        try:
            if (self._is_synonymous(var_data) and self._in_important_domain(var_data)) or (
                self._is_intronic(var_data)
                and not self._affect_canonical_ss(seqvar, var_data)
                and not self._is_conserved(var_data)
            ):
                self.comment_bp7 += (
                    "Synonymous variant is in an important domain or intronic variant is not "
                    "conserved. BP7 is met."
                )
                self.prediction_bp7.BP7 = True
            else:
                self.comment_bp7 += (
                    "Variant is not synonymous or not in an important domain, or intronic variant "
                    "is conserved. BP7 is not met."
                )
                self.prediction_bp7.BP7 = False

        except AutoAcmgBaseException as e:
            logger.error("Failed to predict BP7 criterion. Error: {}", e)
            self.comment_bp7 = f"Failed to predict BP7 criterion. Error: {e}"
            self.prediction_bp7 = None

        return self.prediction_bp7, self.comment_bp7
