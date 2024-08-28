"""
Predictor for DICER1 and miRNA-Processing Gene VCEP.
Included gene: DICER1 (HGNC:17098).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN024
"""

from typing import Dict, List, Tuple

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

PM1_CLUSTER: Dict[str, Dict[str, Dict[str, List]]] = {
    # DICER1
    "HGNC:17098": {
        "moderate": {
            "residues": [1344, 1705, 1709, 1713, 1809, 1810, 1813],  # Metal ion-binding residues
        },
        "supporting": {
            "domains": [(1682, 1846)],  # RNase IIIb domain excluding the metal ion-binding residues
        },
    },
}


class DICER1Predictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override predict_pm1 to include VCEP-specific logic for DICER1."""
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotMet,
                strength=AutoACMGStrength.PathogenicModerate,
                summary="Variant does not affect a critical domain for DICER1.",
            )

        # Check moderate level criteria for metal ion-binding residues
        if var_data.prot_pos in gene_cluster.get("moderate", {}).get("residues", []):
            comment = (
                f"Variant affects a metal ion-binding residue in DICER1 "
                f"(position {var_data.prot_pos}). PM1 is met at the Moderate level."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=comment,
            )

        # Check supporting level criteria for RNase IIIb domain excluding metal ion-binding residues
        for start, end in gene_cluster.get("supporting", {}).get("domains", []):
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant affects a residue in the RNase IIIb domain of DICER1 "
                    f"(position {var_data.prot_pos}) outside the key metal ion-binding residues. "
                    f"PM1 is met at the Supporting level."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicSupporting,
                    summary=comment,
                )

        # If no criteria match
        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary="Variant does not affect a critical domain for DICER1.",
        )

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGData) -> bool:
        """Override BP3 to be not applicable for DICER1."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override PP2 and BP1 to be not applicable for DICER1."""
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

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override donor and acceptor positions for DICER1 VCEP."""
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
