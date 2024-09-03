"""
Predictor for Hereditary Hemorrhagic Telangiectasia VCEP.
Included genes:
ACVRL1 (HGNC:175),
ENG (HGNC:3349).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN135
https://cspec.genome.network/cspec/ui/svi/doc/GN136
"""

from typing import Dict, List, Tuple, Union

from loguru import logger

from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

#: VCEP specifications for Hereditary Hemorrhagic Telangiectasia.
SPECs: List[VcepSpec] = [
    VcepSpec(
        identifier="GN135",
        version="1.1.0",
    ),
    VcepSpec(
        identifier="GN136",
        version="1.1.0",
    ),
]

# fmt: off
PM1_CLUSTER: Dict[str, Dict[str, List[Union[Tuple[int, int], int]]]] = {
    "HGNC:175": {
        "domains": [
            (209, 216),  # glycine-rich loop
            (229, 229),  # phosphate anchor
            (242, 242),  # C-helix E pairing the phosphate anchor
            (329, 335),  # catalytic loop
            (348, 351),  # metal-binding loop
            # BMP10 interaction cluster (individual residues)
            40, 54, 56, 57, 58, 59, 66, 71, 72, 73, 75, 76, 78, 79, 80, 82, 83, 84, 85, 87,
        ]
    },
    "HGNC:3349": {
        "domains": [
            278, 282,  # BMP9 binding sites
            207, 363, 382, 412, 549,  # Pathogenic or likely pathogenic cysteine residues
            350, 394, 516, 582,  # Cysteine residues critical for ENG function
        ]
    },
}


class HHTPredictor(DefaultSeqVarPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to include VCEP-specific logic for Hereditary Hemorrhagic
        Telangiectasia. Include critical domains for ACVRL1 and ENG.
        """
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check moderate level criteria
        moderate_criteria = gene_cluster.get("domains", [])
        if var_data.prot_pos in moderate_criteria or any(
            isinstance(region, tuple) and region[0] <= var_data.prot_pos <= region[1]
            for region in moderate_criteria
        ):
            comment = (
                f"Variant falls within a critical residue for {var_data.hgnc_id}. "
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

    def _bs2_not_applicable(self, var_data: AutoACMGSeqVarData) -> bool:
        """BS2 is not applicable for ACVRL1 and ENG."""
        return True

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.00004
        var_data.thresholds.ba1_benign = 0.01
        var_data.thresholds.bs1_benign = 0.0008
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 is not applicable for Hereditary Hemorrhagic Telangiectasia."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pp2bp1 to return not applicable status for ACVRL1 and ENG."""
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
