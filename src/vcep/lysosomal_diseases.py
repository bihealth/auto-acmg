"""
Predictor for Lysosomal Diseases VCEP.
Included gene: GAA (HGNC:4065).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN010
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

#: VCEP specification for Lysosomal Diseases.
SPEC: VcepSpec = VcepSpec(
    identifier="GN010",
    version="2.0.0",
)

PM1_CLUSTER = {
    "HGNC:4065": {
        "residues": [
            282,  # D282
            376,  # W376
            404,  # D404
            405,  # L405
            441,  # I441
            481,  # W481
            516,  # W516
            518,  # D518
            519,  # M519
            600,  # R600
            613,  # W613
            616,  # D616
            618,  # W618
            649,  # F649
            650,  # L650
            674,  # H674
        ]
    }
}


class LysosomalDiseasesPredictor(DefaultSeqVarPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override predict_pm1 to specify critical domains for Lysosomal Diseases."""
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        if var_data.prot_pos in gene_cluster["residues"]:
            comment = (
                f"Variant affects a residue in GAA at position {var_data.prot_pos}, "
                f"which is associated with Lysosomal Diseases."
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
            summary="Variant does not meet the PM1 criteria for Lysosomal Diseases.",
        )

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.001
        var_data.thresholds.ba1_benign = 0.01
        var_data.thresholds.bs1_benign = 0.005
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 is not applicable for GAA."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to return not applicable status for GAA."""
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
