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

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.seqvar import SeqVar

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


class CerebralCreatineDeficiencySyndromesPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override predict_pm1 to return a not applicable status for PM1."""
        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotApplicable,
            strength=AutoACMGStrength.PathogenicModerate,
            summary="PM1 is not applicable for Cerebral Creatine Deficiency Syndromes.",
        )

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGData) -> bool:
        """BP3 is not applicable for Cerebral Creatine Deficiency Syndromes."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
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
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Predict PP3 and BP4 criteria for Cerebral Creatine Deficiency Syndromes based on VCEP specific rules."""
        logger.info("Predict PP3 and BP4 for Cerebral Creatine Deficiency Syndromes.")

        pp3_met = False
        bp4_met = False
        comments_pp3 = []
        comments_bp4 = []

        # Check for REVEL score implications
        revel_score = var_data.scores.dbnsfp.revel
        if revel_score:
            if revel_score >= 0.75:
                pp3_met = True
                comments_pp3.append(f"REVEL score {revel_score} >= 0.75, meeting PP3.")
            if revel_score <= 0.15:
                bp4_met = True
                comments_bp4.append(f"REVEL score {revel_score} <= 0.15, meeting BP4.")

        # Check for in-frame indels
        if self._is_inframe_indel(var_data):
            provean_score = var_data.scores.dbnsfp.provean
            mutation_taster_score = var_data.scores.dbnsfp.mutationTaster
            if provean_score and mutation_taster_score:
                # Checking arbitrary deleterious and benign thresholds for demo purposes
                if provean_score < -2.5 and mutation_taster_score > 0.5:
                    pp3_met = True
                    comments_pp3.append(
                        "In-frame indel predicted deleterious by PROVEAN and MutationTaster."
                    )
                if provean_score > -1.5 and mutation_taster_score < 0.5:
                    bp4_met = True
                    comments_bp4.append(
                        "In-frame indel predicted benign by PROVEAN and MutationTaster."
                    )

        # Check for splicing implications using SpliceAI
        splice_impact = self._affect_spliceAI(var_data)
        if splice_impact:
            pp3_met = True
            comments_pp3.append("Splicing predictions indicate an impact, meeting PP3.")
        else:
            bp4_met = True
            comments_bp4.append("No significant splicing impact predicted, meeting BP4.")

        # Compile PP3 and BP4 results
        pp3_result = AutoACMGCriteria(
            name="PP3",
            prediction=AutoACMGPrediction.Met if pp3_met else AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicSupporting,
            summary=" | ".join(comments_pp3),
        )
        bp4_result = AutoACMGCriteria(
            name="BP4",
            prediction=AutoACMGPrediction.Met if bp4_met else AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.BenignSupporting,
            summary=" | ".join(comments_bp4),
        )

        return (pp3_result, bp4_result)

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override donor and acceptor positions for Cerebral Creatine Deficiency Syndromes VCEP."""
        if var_data.hgnc_id == "HGNC:4136":
            var_data.thresholds.bp7_donor = 7
            var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
