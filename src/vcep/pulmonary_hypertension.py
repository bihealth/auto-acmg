"""
Predictor for Pulmonary Hypertension VCEP.
Included gene: BMPR2 (HGNC:1078).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN125
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

#: VCEP specification for Pulmonary Hypertension.
SPEC: VcepSpec = VcepSpec(
    identifier="GN125",
    version="1.1.0",
)

# fmt: off
PM1_CLUSTER_BMPR2 = {
    "HGNC:1078": {  # BMPR2
        "strong": [
            # Extracellular domain critical residues
            34, 60, 66, 84, 94, 99, 116, 117, 118, 123,
            # Kinase domain critical residues
            210, 212, 230, 245, 333, 338, 351, 353, 386, 405, 410, 491,
            # KD heterodimerization critical residues
            485, 486, 487, 488, 489, 490, 491, 492,
        ],
        "moderate":
            list(range(33, 132))  # Extracellular domain: 33-131
            + list(range(203, 505)),  # Kinase domain: 203-504
        "non_critical": [
            42, 47, 82, 102, 107, 182, 186, 503, 899  # Non-critical residues
        ],
    }
}


class PulmonaryHypertensionPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Pulmonary Hypertension.
        """
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER_BMPR2.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check if the variant falls within the strong level critical residues
        if var_data.prot_pos in gene_cluster["strong"]:
            comment = (
                f"Variant affects a critical residue in BMPR2 at position {var_data.prot_pos}. "
                f"PM1 is met at the Strong level."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicStrong,
                summary=comment,
            )

        # Check if the variant falls within the moderate level critical residues
        if (
            var_data.prot_pos in gene_cluster["moderate"] and
            not var_data.prot_pos in gene_cluster["non_critical"]
        ):
            comment = (
                f"Variant affects a residue in BMPR2 at position {var_data.prot_pos} "
                f"within the extracellular or kinase domain. PM1 is met at the Moderate level."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=comment,
            )

        # If the variant falls in a non-critical residue
        if var_data.prot_pos in gene_cluster["non_critical"]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotMet,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=(
                    f"Variant affects a residue at position {var_data.prot_pos} in BMPR2, "
                    f"which is demonstrated to be non-critical for kinase activity."
                ),
            )

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary="Variant does not meet the PM1 criteria for BMPR2.",
        )

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to include VCEP-specific logic for BMPR2."""
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
        Add an exception for Pulmonary Hypertension.

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
