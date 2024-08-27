"""
Predictor for Coagulation Factor Deficiency VCEP.
Included genes:
F8 (HGNC:3546),
F9 (HGNC:3551).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN071
https://cspec.genome.network/cspec/ui/svi/doc/GN080
"""

from typing import Dict, List

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.exceptions import AlgorithmError
from src.defs.seqvar import SeqVar

PM1_CLUSTER: Dict[str, Dict[str, Dict[str, List]]] = {
    "HGNC:3546": {  # F8
        "strong": {
            "residues": [
                (391, 392),
                (759, 760),
                (1701, 1705),
                (1708, 1709),
                (1683, 1683),
                (1689, 1689),
                (737, 737),
                (742, 742),
            ]
        },
        "moderate": {
            "residues": [(1667, 1667), (1332, 1332)],  # Residues affecting secretion
            "domains": [(2267, 2304)],  # FXa-binding residues, excluding Ser2283
            "excluded_residues": [(2283, 2283)],  # Excluded residue in FXa-binding region
        },
    },
    "HGNC:3551": {  # F9
        "strong": {
            "residues": [
                (267, 267),
                (315, 315),
                (411, 411),  # Catalytic residues
                (191, 192),
                (226, 227),  # Activation residues
            ]
        },
        "moderate": {"exons": [3, 4, 5]},  # Moderate level for variants within exons 3, 4, or 5
    },
}


class CoagulationFactorDeficiencyPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Coagulation Factor Deficiency.
        """
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            raise AlgorithmError(
                f"PM1 criteria of CoagulationFactorDeficiencyPredictor VCEP "
                f"is not defined for {var_data.hgnc_id}."
            )

        # Check strong level criteria
        for start, end in gene_cluster.get("strong", {}).get("residues", []):
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant affects a critical residue for {var_data.hgnc_id} "
                    f"between positions {start}-{end}. PM1 is met at the Strong level."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicStrong,
                    summary=comment,
                )

        # Check moderate level criteria for residues
        for start, end in gene_cluster.get("moderate", {}).get("residues", []):
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant affects a moderate-level residue for {var_data.hgnc_id} "
                    f"between positions {start}-{end}. PM1 is met at the Moderate level."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicModerate,
                    summary=comment,
                )

        # Check moderate level criteria for domains, excluding specific residues
        for start, end in gene_cluster.get("moderate", {}).get("domains", []):
            if start <= var_data.prot_pos <= end:
                excluded_residues = gene_cluster.get("moderate", {}).get("excluded_residues", [])
                if not any(start <= var_data.prot_pos <= end for start, end in excluded_residues):
                    comment = (
                        f"Variant falls within a moderate-level region for {var_data.hgnc_id} "
                        f"between positions {start}-{end}. PM1 is met at the Moderate level."
                    )
                    return AutoACMGCriteria(
                        name="PM1",
                        prediction=AutoACMGPrediction.Met,
                        strength=AutoACMGStrength.PathogenicModerate,
                        summary=comment,
                    )

        # Check moderate level criteria for exons in F9
        affected_exon = self._get_affected_exon(var_data, seqvar)
        if var_data.hgnc_id == "HGNC:3551" and affected_exon in [3, 4, 5]:
            comment = (
                f"Variant falls within exon {affected_exon} for {var_data.hgnc_id}. "
                f"PM1 is met at the Moderate level."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=comment,
            )

        # If no criteria match
        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary=f"Variant does not meet the PM1 criteria for {var_data.hgnc_id}.",
        )

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Change the spliceAI and phyloP threshold for BP7."""
        if var_data.hgnc_id == "HGNC:3546":
            var_data.thresholds.spliceAI_acceptor_gain = 0.05
            var_data.thresholds.spliceAI_acceptor_loss = 0.05
            var_data.thresholds.spliceAI_donor_gain = 0.05
            var_data.thresholds.spliceAI_donor_loss = 0.05
        elif var_data.hgnc_id == "HGNC:3551":
            var_data.thresholds.spliceAI_acceptor_gain = 0.01
            var_data.thresholds.spliceAI_acceptor_loss = 0.01
            var_data.thresholds.spliceAI_donor_gain = 0.01
            var_data.thresholds.spliceAI_donor_loss = 0.01
        var_data.thresholds.phyloP100 = 0.1
        return super().predict_bp7(seqvar, var_data)
