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

from typing import List, Optional, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    PS1PM5,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.exceptions import AlgorithmError
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

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


class HearingLossPredictor(DefaultSeqVarPredictor):

    def verify_ps1pm5(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PS1PM5], str]:
        """Override PS1/PM5 for Hearing Loss."""
        self.prediction_ps1pm5, self.comment_ps1pm5 = super().verify_ps1pm5(seqvar, var_data)
        if (
            self.prediction_ps1pm5
            and self._is_missense(var_data)
            and self._affect_splicing(var_data)
        ):
            self.prediction_ps1pm5.PS1 = False
            self.comment_ps1pm5 = "Variant affects splicing. PS1 is not applicable."
            self.prediction_ps1pm5.PM5 = False
            self.comment_ps1pm5 = "Variant affects splicing. PM5 is not applicable."
        return self.prediction_ps1pm5, self.comment_ps1pm5

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
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

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.00007
        var_data.thresholds.ba1_benign = 0.001
        var_data.thresholds.bs1_benign = 0.0007
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
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
