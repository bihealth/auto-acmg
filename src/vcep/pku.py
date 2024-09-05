"""
Predictor for Phenylketonuria (PKU) VCEP.
Included gene: PAH (HGNC:8582).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN006
"""

from typing import Optional, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    PP3BP4,
    AlleleCondition,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    GenomicStrand,
    VcepSpec,
)
from src.defs.exceptions import AutoAcmgBaseException, MissingDataError
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

#: VCEP specification for Phenylketonuria.
SPEC: VcepSpec = VcepSpec(
    identifier="GN006",
    version="2.0.0",
)

PM1_CLUSTER_PKU = {
    "HGNC:8582": {  # PAH
        "residues": [
            # active site residues
            138,  # Tyr138
            158,  # Arg158
            245,  # Val245
            268,  # Tyr268
            278,  # Thr278
            279,  # Pro279
            289,  # Glu289
            300,  # Ala300
            315,  # Asp315
            331,  # Phe331
            345,  # Ala345
            346,  # Gly346
            349,  # Ser349
            377,  # Tyr377
        ]
        # substrate binding residues
        + list(range(46, 49))
        + list(range(63, 70))
        # cofactor binding residues
        + [
            285,  # His285
            290,  # His290
            330,  # Glu330
        ]
        + list(range(246, 267))
        + list(range(280, 284))
        + list(range(322, 327))
        + list(range(377, 380)),  # Additional ranges for cofactor binding
    }
}


class PKUPredictor(DefaultSeqVarPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override predict_pm1 to include domain information for PAH."""
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER_PKU.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check moderate level criteria
        if var_data.prot_pos in gene_cluster.get("residues", []):
            comment = (
                f"Variant affects a residue at position {var_data.prot_pos} "
                f"in {var_data.hgnc_id}, which is critical for PKU."
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
            summary="Variant does not meet the PM1 criteria for PAH.",
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
        var_data.thresholds.pm2_pathogenic = 0.000002
        var_data.thresholds.ba1_benign = 0.015
        var_data.thresholds.bs1_benign = 0.002
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 is not applicable for PKU."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to return not applicable status for PAH."""
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
            self.prediction_pp3bp4.BP4 = self.prediction_pp3bp4.BP4 and not self._affect_spliceAI(
                var_data
            )

        except AutoAcmgBaseException as e:
            self.comment_pp3bp4 = f"An error occurred during prediction. Error: {e}"
            self.prediction_pp3bp4 = None
        return self.prediction_pp3bp4, self.comment_pp3bp4

    def _is_bp7_exception(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """
        Add an exception for Phenylketonuria.

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
        """Override donor and acceptor positions for PAH gene."""
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
