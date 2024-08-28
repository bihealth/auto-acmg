"""
Predictor for Glaucoma VCEP.
Included gene: MYOC (HGNC:7610).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN019
"""

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.exceptions import MissingDataError
from src.defs.seqvar import SeqVar


class GlaucomaPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override predict_pm1 to include VCEP-specific logic for Glaucoma."""
        logger.info("Predict PM1")

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotApplicable,
            strength=AutoACMGStrength.PathogenicSupporting,
            summary="PM1 is not applicable for MYOC.",
        )

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGData) -> bool:
        """Override _bp3_not_applicable for Glaucoma VCEP."""
        return True

    def _is_conserved(self, var_data: AutoACMGData) -> bool:
        """
        Predict if the variant is conserved.

        Check if the variant is conserved using the GERP score.

        Args:
            variant_info: The variant information.

        Returns:
            bool: True if the variant is conserved, False otherwise.

        Raises:
            MissingDataError: If the GERP score is missing
        """
        if not var_data.scores.cadd.gerp:
            raise MissingDataError("GERP score is missing.")
        return var_data.scores.cadd.gerp >= var_data.thresholds.gerp

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Change BP7 thresholds for Glaucoma VCEP."""
        var_data.thresholds.spliceAI_acceptor_gain = 0.2
        var_data.thresholds.spliceAI_acceptor_loss = 0.2
        var_data.thresholds.spliceAI_donor_gain = 0.2
        var_data.thresholds.spliceAI_donor_loss = 0.2
        return super().predict_bp7(seqvar, var_data)
