"""
Predictor for Myeloid Malignancy VCEP.
Included gene: RUNX1 (HGNC:10471).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN008
"""

from typing import Optional, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    PM2BA1BS1BS2,
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

#: VCEP specification for Myeloid Malignancy.
SPEC: VcepSpec = VcepSpec(
    identifier="GN008",
    version="2.0.0",
)

PM1_CLUSTER = {
    "HGNC:10471": {  # RUNX1
        "moderate": [  # Specific residues within the RHD (Runt Homology Domain)
            107,  # R107
            110,  # K110
            134,  # A134
            162,  # R162
            166,  # R166
            167,  # S167
            169,  # R169
            170,  # G170
            194,  # K194
            196,  # T196
            198,  # D198
            201,  # R201
            204,  # R204
        ],
        "supporting": list(range(89, 205)),  # Other residues within the RHD (89-204)
    }
}


class MyeloidMalignancyPredictor(DefaultSeqVarPredictor):

    def _is_allowed_nonsense(self, var_data: AutoACMGSeqVarData) -> bool:
        """
        Check if the variant is a nonsense/frameshift variant and is downstream of c.98.

        Args:
            var_data (AutoACMGSeqVarData): The variant information.

        Returns:
            bool: True if the variant is a nonsense/frameshift variant and is downstream of c.98.
        """
        allowed = False
        # Check the nonsense
        if "nonsense" in var_data.consequence.cadd:
            allowed = True
        if any("nonsense" in cons for cons in var_data.consequence.mehari):
            allowed = True

        # Check the frameshift
        if "frameshift" in var_data.consequence.cadd:
            allowed = True
        if any("frameshift" in cons for cons in var_data.consequence.mehari):
            allowed = True

        # Check the protein position (downstream of c.98)
        if allowed and var_data.cds_pos > 98:
            allowed = True
        else:
            allowed = False
        return allowed

    def verify_ps1pm5(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PS1PM5], str]:
        """Override PS1/PM5 for InSIGHT Hereditary Colorectal Cancer/Polyposis."""
        self.prediction_ps1pm5, self.comment_ps1pm5 = super().verify_ps1pm5(seqvar, var_data)
        # Check nonsense/frameshift variants for PM5
        try:
            if not self.prediction_ps1pm5 or not self._is_allowed_nonsense(var_data):
                raise AlgorithmError("Variant is not nonsense/frameshift variant.")

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
        """Override PM1 to specify critical residues for Myeloid Malignancy VCEP."""
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check moderate level criteria
        if var_data.prot_pos in gene_cluster["moderate"]:
            comment = (
                f"Variant affects a critical residue within the Runt Homology Domain (RHD) "
                f"at position {var_data.prot_pos} in RUNX1. PM1 is met at the Moderate level."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=comment,
            )

        # Check supporting level criteria
        if var_data.prot_pos in gene_cluster["supporting"]:
            comment = (
                f"Variant affects a residue within the Runt Homology Domain (RHD) "
                f"at position {var_data.prot_pos} in RUNX1. PM1 is met at the Supporting level."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=comment,
            )

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary="Variant does not meet the PM1 criteria for RUNX1.",
        )

    def _bs2_not_applicable(self, var_data: AutoACMGSeqVarData) -> bool:
        """BS2 is not applicable for Myeloid Malignancy VCEP."""
        return True

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

            BA1: Allele frequency is greater than 0.15%.

            BS1: Allele frequency is between 0.015% and 0.15%.

            BS2: No change from the reference implementation.

        Returns:
            BA1BS1BS2PM2: The prediction result.
        """
        self.prediction_pm2ba1bs1bs2 = PM2BA1BS1BS2()
        self.comment_pm2ba1bs1bs2 = ""
        var_data.thresholds.ba1_benign = 0.0015
        var_data.thresholds.bs1_benign = 0.00015
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
        """BP3 is not applicable for Myeloid Malignancy VCEP."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to return not applicable status for RUNX1."""
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
        """Change BP7 thresholds for Myeloid Malignancy VCEP."""
        var_data.thresholds.spliceAI_acceptor_gain = 0.2
        var_data.thresholds.spliceAI_acceptor_loss = 0.2
        var_data.thresholds.spliceAI_donor_gain = 0.2
        var_data.thresholds.spliceAI_donor_loss = 0.2
        var_data.thresholds.phyloP100 = 2.0
        return super().predict_bp7(seqvar, var_data)
