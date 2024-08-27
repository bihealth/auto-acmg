"""
Predictor for ENIGMA BRCA1 and BRCA2 VCEP.
Included genes:
BRCA1 (HGNC:1100),
BRCA2 (HGNC:1101).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN092
https://cspec.genome.network/cspec/ui/svi/doc/GN097
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
from src.defs.exceptions import AutoAcmgBaseException
from src.defs.seqvar import SeqVar

BP7_IMPORTANT_DOMAINS = {
    "HGNC:1100": [
        (2, 102),  # RING domain
        (1391, 1423),  # coiled-coil aa
        (1650, 1857),  # BRCT domain
    ],
    "HGNC:1101": [
        (10, 40),  # PALB2 binding
        (2481, 3186),  # DNA binding
    ],
}


class ENIGMAPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override PM1 prediction for ENIGMA BRCA1 and BRCA2."""
        logger.info("Predict PM1")

        if var_data.hgnc_id in ["HGNC:1100", "HGNC:1101"]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicModerate,
                summary="PM1 is not applicable for BRCA1 and BRCA2.",
            )

        return super().predict_pm1(seqvar, var_data)

    def _in_important_domain(self, var_data: AutoACMGData) -> bool:
        """Check if the variant is in an important domain."""
        for start, end in BP7_IMPORTANT_DOMAINS.get(var_data.hgnc_id, []):
            if start <= var_data.prot_pos <= end:
                return True
        return False

    def verify_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> Tuple[Optional[BP7], str]:
        """Override verify BP7 criterion for ENIGMA BRCA1 and BRCA2."""
        # Change the donor/acceptor thresholds for BRCA1 and BRCA2
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21

        self.prediction_bp7 = BP7()
        self.comment_bp7 = ""
        try:
            if (self._is_synonymous(var_data) and self._in_important_domain(var_data)) or (
                self._is_intronic(var_data)
                and not self._affect_canonical_ss(seqvar, var_data)
                and not self._is_conserved(var_data)
            ):
                self.comment_bp7 += (
                    "Synonymous variant is in an important domain or intronic variant is not "
                    "conserved. BP7 is met."
                )
                self.prediction_bp7.BP7 = True
            else:
                self.comment_bp7 += (
                    "Variant is not synonymous or not in an important domain, or intronic variant "
                    "is conserved. BP7 is not met."
                )
                self.prediction_bp7.BP7 = False

        except AutoAcmgBaseException as e:
            logger.error("Failed to predict BP7 criterion. Error: {}", e)
            self.comment_bp7 = f"Failed to predict BP7 criterion. Error: {e}"
            self.prediction_bp7 = None

        return self.prediction_bp7, self.comment_bp7
