"""
Predictor for ACADVL VCEP.
Included gene: ACADVL (HGNC:92).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN021
"""

from typing import Tuple

from loguru import logger

from src.defs.auto_acmg import (
    PP3BP4,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGSeqVarResult,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

#: VCEP specification for ACADVL.
SPEC: VcepSpec = VcepSpec(
    identifier="GN021",
    version="1.0.0",
)

PM1_CLUSTER = [
    (214, 223),  # Nucleotide and substrate binding
    (249, 251),  # Nucleotide and substrate binding
    (460, 466),  # Nucleotide and substrate binding
    (481, 516),  # Membrane binding
    (1, 40),  # Mitochondrial signal peptide
]


class ACADVLPredictor(DefaultSeqVarPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override predict_pm1 to return a not applicable status for PM1."""
        logger.info("Predict PM1")
        # Check if variant falls within critical regions
        for start, end in PM1_CLUSTER:
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant falls within a critical region for ACADVL between positions "
                    f"{start}-{end}. PM1 is met."
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
            summary="Variant does not fall within any critical region for ACADVL. PM1 is not met.",
        )

    def _bs2_not_applicable(self, var_data: AutoACMGSeqVarData) -> bool:
        """BS2 is not applicable for ACADVL."""
        return True

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.001
        var_data.thresholds.ba1_benign = 0.007
        var_data.thresholds.bs1_benign = 0.0035
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 is not applicable for ACADVL."""
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

    def predict_pp3bp4(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Predict PP3 and BP4 criteria based on ACADVL VCEP specific rules."""
        logger.info("Predict PP3 and BP4.")

        pp3_met = False
        bp4_met = False
        comments = []
        self.prediction_pp3bp4 = PP3BP4()

        # Evaluate missense changes
        if self._is_missense_variant(var_data):
            revel_score = var_data.scores.dbnsfp.revel
            if revel_score:
                if revel_score > 0.75:
                    pp3_met = True
                    comments.append(f"REVEL score {revel_score} > 0.75, PP3 met.")
                if revel_score < 0.5:
                    bp4_met = True
                    comments.append(f"REVEL score {revel_score} < 0.5, BP4 met.")

        # Evaluate in-frame deletions and insertions
        if self._is_inframe_indel(var_data):
            provean = var_data.scores.dbnsfp.provean
            mutation_taster = var_data.scores.dbnsfp.mutationTaster
            # Assume thresholds for PROVEAN and Mutation Taster are set properly elsewhere
            if provean and mutation_taster:  # Need actual logic to determine thresholds
                pp3_met = pp3_met or (provean < -2.5 and mutation_taster > 0.5)
                bp4_met = bp4_met or (provean > -2.5 and mutation_taster < 0.5)

        # Evaluate splice variants
        if self._is_splice_variant(var_data):
            self.comment_pp3bp4 = "Variant is a splice variant."
            self.prediction_pp3bp4.PP3 = self._is_pathogenic_splicing(var_data)
            self.prediction_pp3bp4.BP4 = self._is_benign_splicing(var_data)
            self.comment_pp3bp4 += (
                f"Ada score: {var_data.scores.dbscsnv.ada}, "
                f"Ada threshold: {var_data.thresholds.ada}. "
                f"RF score: {var_data.scores.dbscsnv.rf}, "
                f"RF threshold: {var_data.thresholds.rf}. "
            )

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
