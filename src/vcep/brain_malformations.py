"""
Predicting brain malformations VCEP.
Included genes: AKT3 (HGNC:393), MTOR (HGNC:3942), PIK3CA (HGNC:8975), PIK3R2 (HGNC:8980)
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN018
"""

from typing import Optional, Tuple

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import (
    BP7,
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
)
from src.defs.seqvar import SeqVar

PM1_CLUSTER = {
    "HGNC:393": [  # AKT3
        (5, 109),  # Pleckstrin homology domain
        (151, 388),  # Catalytic kinase domain
        (425, 475),  # C-terminal Protein Kinase
    ],
    "HGNC:3942": [  # MTOR
        (1382, 1982),  # Kinase domain
        (2015, 2114),  # FKBP-rapamycin-binding (FRB) domain
    ],
    "HGNC:8975": [  # PIK3CA
        (173, 292),  # Kinase Ras-binding domain
        (322, 483),  # Kinase domains
        (797, 1068),  # Kinase domains
    ],
    "HGNC:8980": [  # PIK3R2
        (328, 716),  # SH2, sequence homology 2 domain
        (31, 108),  # Adaptor binding domain (PI3K ABD)
    ],
}


class BrainMalformationsPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        "Override predict_pm1 to include VCEP-specific logic for brain malformations VCEP."
        logger.info("Predict PM1")

        if var_data.hgnc_id in PM1_CLUSTER:
            domains = PM1_CLUSTER[var_data.hgnc_id]

            # Check if the variant falls within any of the specified domains
            for start_aa, end_aa in domains:
                if start_aa <= var_data.prot_pos <= end_aa:
                    comment = (
                        f"Variant falls within a critical domain for {var_data.hgnc_id} "
                        f"between positions {start_aa}-{end_aa}. PM1 is met."
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
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=(
                    "Variant does not fall within any critical domain for the specified gene. "
                    "PM1 is not met."
                ),
            )
        else:
            return super().predict_pm1(seqvar, var_data)

    def predict_pm4bp3(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override predict_pm4bp3 to include VCEP-specific logic for brain malformations VCEP."""
        logger.info("Predict PM4 and BP3")
        return (
            AutoACMGCriteria(
                name="PM4",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicModerate,
                summary="PM4 is not applicable for the brain malformations VCEP.",
            ),
            AutoACMGCriteria(
                name="BP3",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.BenignSupporting,
                summary="BP3 is not applicable for the brain malformations VCEP.",
            ),
        )

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Change the PhyloP100 score threshold for BP7."""
        var_data.thresholds.phyloP100 = 0.1
        return super().predict_bp7(seqvar, var_data)
