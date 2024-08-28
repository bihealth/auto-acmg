"""
Predictor for Severe Combined Immunodeficiency Disease VCEP.
Included genes:
FOXN1 (HGNC:12765),
ADA (HGNC:186),
DCLRE1C (HGNC:17642),
IL7R (HGNC:6024),
JAK3 (HGNC:6193),
RAG1 (HGNC:9831),
RAG2 (HGNC:9832),
IL2RG (HGNC:6010).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN113
https://cspec.genome.network/cspec/ui/svi/doc/GN114
https://cspec.genome.network/cspec/ui/svi/doc/GN116
https://cspec.genome.network/cspec/ui/svi/doc/GN119
https://cspec.genome.network/cspec/ui/svi/doc/GN121
https://cspec.genome.network/cspec/ui/svi/doc/GN123
https://cspec.genome.network/cspec/ui/svi/doc/GN124
https://cspec.genome.network/cspec/ui/svi/doc/GN129
"""

from typing import Dict, List, Tuple, Union

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

# fmt: off
PM1_CLUSTER_SCID: Dict[str, Dict[str, List[Union[int, Tuple[int, int]]]]] = {
    "HGNC:12765": {  # FOXN1
        "domains": [(270, 367)],  # DNA binding forkhead domain
    },
    "HGNC:6193": {  # JAK3
        "domains": [(651, 759)],  # JH2 domain residues R651W and C759R
    },
    "HGNC:9831": {  # RAG1
        "domains": [(394, 460), (461, 517)],  # NBD domain, DDBD domain
        "supporting_domains": [(387, 1011)],  # Core domain (supporting strength)
    },
    "HGNC:9832": {  # RAG2
        "domains": [(414, 487)],  # PHD domain
        "supporting_domains": [(1, 383)],  # Core domain (supporting strength)
    },
    "HGNC:6010": {  # IL2RG
        "strong_residues": [
            62, 72, 102, 115,  # Conserved cysteine residues
            224, 226, 691, 285,  # CpG dinucleotides
            237, 238, 239, 240, 241,  # WSxWS motif
        ],
        "strong_domains": [(263, 283)],  # Transmembrane domain
    },
}


class SCIDPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Severe Combined Immunodeficiency Disease.
        """
        logger.info("Predict PM1")

        if var_data.hgnc_id in [
            "HGNC:186",  # ADA
            "HGNC:17642",  # DCLRE1C
            "HGNC:6024",  # IL7R
        ]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=f"PM1 is not applicable for {var_data.hgnc_id}.",
            )

        gene_info = PM1_CLUSTER_SCID.get(var_data.hgnc_id, None)
        if not gene_info:
            return super().predict_pm1(seqvar, var_data)

        # Check if the variant falls within a critical region
        for start, end in gene_info.get("domains", []):   # type: ignore
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant affects a critical residue in {var_data.hgnc_id} "
                    f"within positions {start}-{end}. PM1 is met."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicModerate,
                    summary=comment,
                )

        # Check if the variant falls within a supporting region
        for start, end in gene_info.get("supporting_domains", []):   # type: ignore
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant affects a residue in the supporting region of {var_data.hgnc_id} "
                    f"within positions {start}-{end}. PM1 is met at a supporting level."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicSupporting,
                    summary=comment,
                )

        # Check if the variant affects a strong residue
        if var_data.prot_pos in gene_info.get("strong_residues", []):
            comment = (
                f"Variant affects a strong residue in {var_data.hgnc_id} "
                f"at position {var_data.prot_pos}. PM1 is met."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicStrong,
                summary=comment,
            )

        # Check if the variant falls within a strong domain
        for start, end in gene_info.get("strong_domains", []):   # type: ignore
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant affects a residue in the strong domain of {var_data.hgnc_id} "
                    f"within positions {start}-{end}. PM1 is met at a strong level."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicStrong,
                    summary=comment,
                )

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary=f"Variant does not meet the PM1 criteria for {var_data.hgnc_id}.",
        )

    def _is_conserved(self, var_data: AutoACMGData) -> bool:
        """
        Override the default _is_conserved method to ignore this check for SCID genes.
        """
        if var_data.hgnc_id == "HGNC:12765":
            return super()._is_conserved(var_data)
        return False  # Ignore conservation check for SCID genes

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override donor and acceptor positions for Severe Combined Immunodeficiency Disease genes
        VCEP.
        """
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
