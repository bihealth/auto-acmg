"""
Predictor for Monogenic Diabetes VCEP.
Included genes:
HNF1A (HGNC:11621),
HNF4A (HGNC:5024),
GCK (HGNC:4195).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN017
https://cspec.genome.network/cspec/ui/svi/doc/GN085
https://cspec.genome.network/cspec/ui/svi/doc/GN086
"""

from typing import Dict, List, Tuple, Union

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

# fmt: off
PM1_CLUSTER: Dict[str, Dict[str, Dict[str, List[Union[int, Tuple[int, int]]]]]] = {
    "HGNC:11621": {  # HNF1A
        "moderate": {
            "residues": [
            130, 131, 132, 143, 144, 145, 146, 147, 148, 149,
            155, 156, 157, 158, 203, 204, 205, 206, 263, 264,
            265, 270, 271, 272, 273
        ]},
        "supporting": {
            "domains": [
                (1, 32),  # Dimerization domain
                (107, 174),  # Subset of DNA binding domains
                (201, 280),  # Subset of DNA binding domains
            ],
            "promoter_regions": [
                (-195, -187),  # AP1 binding site
                (-227, -209),  # Overlapping HNF3 & NF-Y sites
                (-259, -238),  # HNF1A binding site
                (-288, -276),  # HNF4A binding site
            ]
        }
    },
    "HGNC:5024": {  # HNF4A
        "moderate": {
            "residues": [
                43, 49, 50, 51, 56, 57, 59, 63, 64, 67, 70, 72,
                75, 87, 88, 91, 94, 109, 112, 113,
                38, 41, 55, 58, 74, 80, 90, 93   # zinc finger residues
            ],
        },
        "supporting": {
            "domains": [
                (37, 113),  # DNA binding domain
                (180, 220),  # Ligand binding domain
                (300, 350),  # Ligand binding domain
            ],
            "promoter_regions": [
                (-151, -132),  # HNF6/OC2 and IPF1 binding sites
                (-181, -169),  # HNF1A/HNF1B binding sites
            ]
        }
    },
    "HGNC:4195": {  # GCK
        "moderate": {
            "residues": [
                151, 153, 168, 169, 204, 205, 225, 231, 254, 258,
                287, 290, 78, 82, 85, 228, 229, 295, 296, 331, 333,
                336, 410, 416
            ]
        }
    }
}


class MonogenicDiabetesPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Monogenic Diabetes.
        """
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check moderate level criteria
        if var_data.prot_pos in gene_cluster.get("moderate", {}).get("residues", []):
            comment = (
                f"Variant affects a residue at position {var_data.prot_pos} "
                f"in {var_data.hgnc_id}, which is critical for Monogenic Diabetes."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=comment,
            )

        # Check supporting level criteria
        if "supporting" in gene_cluster:
            # Check if the variant is in the specified domains
            for start, end in gene_cluster.get("supporting", {}).get("domains", []):   # type: ignore
                if start <= var_data.prot_pos <= end:
                    comment = (
                        f"Variant falls within a critical region in {var_data.hgnc_id} between "
                        f"positions {start}-{end}. PM1 is met at the Supporting level."
                    )
                    return AutoACMGCriteria(
                        name="PM1",
                        prediction=AutoACMGPrediction.Met,
                        strength=AutoACMGStrength.PathogenicSupporting,
                        summary=comment,
                    )

            # Check if the variant is in the promoter regions
            for start, end in gene_cluster.get("supporting", {}).get("promoter_regions", []):   # type: ignore
                if start <= var_data.cds_pos <= end:
                    comment = (
                        "Variant falls within a critical promoter region in "
                        f"{var_data.hgnc_id} between "
                        f"positions {start}-{end}. PM1 is met at the Supporting level."
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
            summary=f"Variant does not meet the PM1 criteria for {var_data.hgnc_id}.",
        )
