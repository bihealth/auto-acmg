"""
Predictor for Malignant Hyperthermia Susceptibility VCEP.
Included gene: RYR1 (HGNC:10483).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN012
"""

from typing import Optional, Tuple

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import (
    PM2BA1BS1BS2,
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.exceptions import AutoAcmgBaseException
from src.defs.seqvar import SeqVar

#: VCEP specification for Malignant Hyperthermia Susceptibility.
SPEC: VcepSpec = VcepSpec(
    identifier="GN012",
    version="2.0.0",
)

PM1_CLUSTER = {
    "HGNC:10483": {
        "moderate": {
            "domains": [
                (1, 552),  # N-terminal region
                (2101, 2458),  # Central region
            ],
        },
        "supporting": {
            "domains": [
                (1, 552),  # N-terminal region (if PS1/PM5 applicable)
                (2101, 2458),  # Central region (if PS1/PM5 applicable)
                (4631, 4991),  # C-terminal region
            ],
        },
    }
}


class MalignantHyperthermiaPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override PM1 to specify critical domains for Malignant Hyperthermia Susceptibility."""
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check moderate level domains
        for start, end in gene_cluster.get("moderate", {}).get("domains", []):
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant falls within a critical region in RYR1 between positions {start}-{end}. "
                    f"PM1 is met at the Moderate level."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicModerate,
                    summary=comment,
                )

        # Check supporting level domains
        for start, end in gene_cluster.get("supporting", {}).get("domains", []):
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant falls within a critical region in RYR1 between positions {start}-{end}. "
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
            summary="Variant does not meet the PM1 criteria for Malignant Hyperthermia Susceptibility.",
        )

    def verify_pm2ba1bs1bs2(
        self,
        seqvar: SeqVar,
        var_data: AutoACMGData,
    ) -> Tuple[Optional[PM2BA1BS1BS2], str]:
        """
        Predicts the PM2, BA1, BS1, BS2 criteria for the sequence variant.

        Note:
            Rules:
            PM2: Not applicable for Malignant Hyperthermia Susceptibility.

            BA1: Allele frequency is greater than 0,0038.

            BS1: Allele frequency is between 0.0008 and 0.0038.

            BS2: No change from the reference implementation.

        Returns:
            BA1BS1BS2PM2: The prediction result.
        """
        self.prediction_pm2ba1bs1bs2 = PM2BA1BS1BS2()
        var_data.thresholds.ba1_benign = 0.0038
        var_data.thresholds.bs1_benign = 0.0008
        try:
            af = self._get_af(seqvar, var_data)
            if not af:
                self.comment_pm2ba1bs1bs2 = "No allele frequency data found. "
            elif self._ba1_exception(seqvar):
                self.comment_pm2ba1bs1bs2 = "The variant is in the exception list for BA1 criteria."
                self.prediction_pm2ba1bs1bs2.BA1 = False
                self.prediction_pm2ba1bs1bs2.BS1 = False
            elif af >= var_data.thresholds.ba1_benign:
                self.comment_pm2ba1bs1bs2 = "Allele frequency > 5%: BA1 is met. "
                self.prediction_pm2ba1bs1bs2.BA1 = True
            elif af >= var_data.thresholds.bs1_benign:
                self.comment_pm2ba1bs1bs2 = "Allele frequency > 1%: BS1 is met. "
                self.prediction_pm2ba1bs1bs2.BS1 = True

            self.comment_pm2ba1bs1bs2 += (
                "PM2 is not applicable for Malignant Hyperthermia Susceptibility. "
            )

            if not self.prediction_pm2ba1bs1bs2.BA1 and af and self._check_zyg(seqvar, var_data):
                self.comment_pm2ba1bs1bs2 += (
                    "The variant is in a recessive, dominant, or X-linked disorder: BS2 is met."
                )
                self.prediction_pm2ba1bs1bs2.BS2 = True

        except AutoAcmgBaseException as e:
            logger.error("Error occurred during PM2, BA1, BS1, BS2 prediction. Error: {}", e)
            self.comment_pm2ba1bs1bs2 = (
                f"An error occurred while predicting PM2, BA1, BS1, BS2 criteria: {e}"
            )
            self.prediction_pm2ba1bs1bs2 = None

        return self.prediction_pm2ba1bs1bs2, self.comment_pm2ba1bs1bs2

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """
        Override PP2 and BP1 for Malignant Hyperthermia Susceptibility to return not applicable
        status.
        """
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
