"""
Predictor for Hearing Loss VCEP.
Included genes:
CDH23 (HGNC:13733),
COCH (HGNC:2180),
GJB2 (HGNC:4284),
KCNQ4 (HGNC:6298),
MYO6 (HGNC:7605),
MYO7A (HGNC:7606),
SLC26A4 (HGNC:8818),
TECTA (HGNC:11720),
USH2A (HGNC:12601),
MYO15A (HGNC:7594),
OTOF (HGNC:8515).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN005
https://cspec.genome.network/cspec/ui/svi/doc/GN023
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
from src.defs.exceptions import AlgorithmError
from src.defs.seqvar import SeqVar

#: VCEP specifications for Hearing Loss.
SPECs: List[VcepSpec] = [
    VcepSpec(
        identifier="GN005",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN023",
        version="1.0.0",
    ),
]

PM1_CLUSTER = {
    "HGNC:6298": {
        "domains": [
            (271, 292),  # Pore-forming intramembrane region
        ]
    }
}


class HearingLossPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override predict_pm1 to include domains for KCNQ4. For other genes - not applicable."""
        logger.info("Predict PM1")

        if var_data.hgnc_id in [
            "HGNC:13733",  # CDH23
            "HGNC:2180",  # COCH
            "HGNC:4284",  # GJB2
            "HGNC:7605",  # MYO6
            "HGNC:7606",  # MYO7A
            "HGNC:8818",  # SLC26A4
            "HGNC:11720",  # TECTA
            "HGNC:12601",  # USH2A
            "HGNC:7594",  # MYO15A
            "HGNC:8515",  # OTOF
        ]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=f"PM1 is not applicable for {var_data.hgnc_id}.",
            )

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check strong level criteria
        for start, end in gene_cluster["domains"]:
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant falls within the critical pore-forming intramembrane region "
                    f"of KCNQ4 between positions {start}-{end}. PM1 is met at the Strong level."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicStrong,
                    summary=comment,
                )

        # If no criteria match
        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary="Variant does not meet the PM1 criteria for KCNQ4.",
        )

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override PP2 and BP1 for Hearing Loss to return not applicable status."""
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
