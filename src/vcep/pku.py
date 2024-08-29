"""
Predictor for Phenylketonuria (PKU) VCEP.
Included gene: PAH (HGNC:8582).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN006
"""

from typing import Tuple

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
    GenomicStrand,
    VcepSpec,
)
from src.defs.seqvar import SeqVar

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


class PKUPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Phenylketonuria (PKU).
        """
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

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGData) -> bool:
        """Override BP3 for Phenylketonuria."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to include VCEP-specific logic for PAH."""
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

    def _is_bp7_exception(self, seqvar: SeqVar, var_data: AutoACMGData) -> bool:
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

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override donor and acceptor positions for PAH gene."""
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
