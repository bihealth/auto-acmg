"""
Predictor for Leber Congenital Amaurosis/early onset Retinal Dystrophy VCEP.
Included gene: RPE65 (HGNC:10294)
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN120
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
)
from src.defs.seqvar import SeqVar

PM1_CLUSTER = {
    "HGNC:10294": {
        "residues": [180, 182, 241, 313, 417, 527] + list(range(107, 126)),
    }
}


class LeberCongenitalAmaurosisPredictor(DefaultPredictor):
    """
    Predictor for Leber Congenital Amaurosis/early onset Retinal Dystrophy VCEP.
    """

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Leber Congenital Amaurosis/early
        onset Retinal Dystrophy.
        """
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

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGData) -> bool:
        """Override BP3 for Leber Congenital Amaurosis."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to include VCEP-specific logic for RPE65."""
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

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override donor and acceptor positions for RPE65 gene."""
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
