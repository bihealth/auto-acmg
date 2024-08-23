"""
Predictor for Rett and Angelman-like Disorders VCEP.
Included genes:
TCF4 (HGNC:11634),
SLC9A6 (HGNC:11079),
CDKL5 (HGNC:11411),
FOXG1 (HGNC:3811),
MECP2 (HGNC:6990),
UBE3A (HGNC:12496).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN032
https://cspec.genome.network/cspec/ui/svi/doc/GN033
https://cspec.genome.network/cspec/ui/svi/doc/GN034
https://cspec.genome.network/cspec/ui/svi/doc/GN035
https://cspec.genome.network/cspec/ui/svi/doc/GN036
https://cspec.genome.network/cspec/ui/svi/doc/GN037
"""

from typing import Dict, List, Tuple, Union

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

# fmt: off
PM1_CLUSTER: Dict[str, Dict[str, List[Tuple[int, int]]]] = {
    "HGNC:11634": {  # TCF4
        "regions": [(564, 617)],  # Basic Helix-Loop-Helix (bHLH) domain
    },
    "HGNC:11411": {  # CDKL5
        "regions": [(19, 43), (169, 171)],  # ATP binding region, TEY phosphorylation site
    },
    "HGNC:3811": {  # FOXG1
        "regions": [(181, 275)],  # Forkhead domain
    },
    "HGNC:6990": {  # MECP2
        "regions": [(90, 162), (302, 306)],  # Methyl-DNA binding (MBD), Transcriptional repression domain (TRD)
    },
    "HGNC:12496": {  # UBE3A
        "regions": [(820, 820)],  # 3’ cysteine binding site
    },
}


class RettAngelmanPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Rett and Angelman-like Disorders.
        """
        logger.info("Predict PM1")

        # PM1 is not applicable for SLC9A6
        if var_data.hgnc_id == "HGNC:11079":
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=f"PM1 is not applicable for SLC9A6.",
            )

        gene_info = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_info:
            return super().predict_pm1(seqvar, var_data)

        # Check if the variant falls within a critical region
        for start, end in gene_info.get("regions", []):
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

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary=f"Variant does not meet the PM1 criteria for {var_data.hgnc_id}.",
        )
