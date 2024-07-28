from typing import Optional

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.criteria.auto_bp7 import AutoBP7
from src.criteria.auto_pm1 import AutoPM1
from src.criteria.auto_pm2_ba1_bs1_bs2 import AutoPM2BA1BS1BS2
from src.criteria.auto_pm4_bp3 import AutoPM4BP3
from src.criteria.auto_pp2_bp1 import AutoPP2BP1
from src.criteria.auto_pp3_bp4 import AutoPP3BP4
from src.criteria.auto_ps1_pm5 import AutoPS1PM5
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.auto_acmg import ACMGPrediction, ACMGResult
from src.defs.exceptions import AutoAcmgBaseException
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


class AutoACMGCriteria:
    """Predict ACMG criteria for sequence variant."""

    def __init__(self, seqvar: SeqVar, *, config: Optional[Config] = None):
        #: Configuration to use.
        self.config: Config = config or Config()
        # Attributes to be set
        self.seqvar: SeqVar = seqvar
        self.annonars_client: AnnonarsClient = AnnonarsClient(
            api_base_url=self.config.api_base_url_annonars
        )
        self.prediction: Optional[ACMGResult] = None

    def _get_variant_info(self, seqvar: SeqVar) -> Optional[AnnonarsVariantResponse]:
        """Get variant information from Annonars.

        Returns:
            AnnonarsVariantResponse: Annonars response.
        """
        try:
            logger.debug("Getting variant information for {}.", seqvar)
            return self.annonars_client.get_variant_info(seqvar)
        except AutoAcmgBaseException as e:
            logger.error("Failed to get variant information. Error: {}", e)
            return None

    def predict(self) -> Optional[ACMGResult]:
        """Predict ACMG criteria for sequence variant."""
        self.prediction = ACMGResult()

        variant_info = self._get_variant_info(self.seqvar)
        if not variant_info:
            logger.error("Failed to get variant information for {}.", self.seqvar)
            return None

        # # PS1 and PM5
        # try:
        #     logger.info("Predicting PS1 and PM5 criteria.")
        #     ps1pm5 = AutoPS1PM5(self.seqvar, variant_info.result, config=self.config)
        #     ps1_pm5_prediction, ps1_pm5_comment = ps1pm5.predict()
        #     if not ps1_pm5_prediction:
        #         logger.error("Failed to predict PS1&PM5 criteria.")
        #         self.prediction.PS1.comment = ps1_pm5_comment
        #         self.prediction.PM5.comment = ps1_pm5_comment
        #     else:
        #         self.prediction.PS1.prediction = (
        #             ACMGPrediction.Met if ps1_pm5_prediction.PS1 else ACMGPrediction.Unmet
        #         )
        #         self.prediction.PM5.prediction = (
        #             ACMGPrediction.Met if ps1_pm5_prediction.PM5 else ACMGPrediction.Unmet
        #         )
        #         self.prediction.PS1.comment = ps1_pm5_comment
        #         self.prediction.PM5.comment = ps1_pm5_comment

        # except AutoAcmgBaseException as e:
        #     logger.error("Failed to predict PS1 and PM5 criteria. Error: {}", e)

        # # PM4 and BP3
        # try:
        #     logger.info("Predicting PM4 and BP3 criteria.")
        #     pm4bp3 = AutoPM4BP3(self.seqvar, variant_info.result, config=self.config)
        #     pm4_bp3_prediction, pm4_bp3_comment = pm4bp3.predict()
        #     if not pm4_bp3_prediction:
        #         logger.error("Failed to predict PM4&BP3 criteria.")
        #         self.prediction.PM4.comment = pm4_bp3_comment
        #         self.prediction.BP3.comment = pm4_bp3_comment
        #     else:
        #         self.prediction.PM4.prediction = (
        #             ACMGPrediction.Met if pm4_bp3_prediction.PM4 else ACMGPrediction.Unmet
        #         )
        #         self.prediction.BP3.prediction = (
        #             ACMGPrediction.Met if pm4_bp3_prediction.BP3 else ACMGPrediction.Unmet
        #         )
        #         self.prediction.PM4.comment = pm4_bp3_comment
        #         self.prediction.BP3.comment = pm4_bp3_comment

        # except AutoAcmgBaseException as e:
        #     logger.error("Failed to predict PM4 and BP3 criteria. Error: {}", e)

        # # BA1, BS1, BS2, PM2
        # try:
        #     logger.info("Predicting BA1, BS1, BS2, and PM2 criteria.")
        #     ba1bs1bs2pm2 = AutoPM2BA1BS1BS2(self.seqvar, variant_info.result, config=self.config)
        #     ba1bs1bs2pm2_prediction, ba1bs1bs2pm2_comment = ba1bs1bs2pm2.predict()
        #     if not ba1bs1bs2pm2_prediction:
        #         logger.error("Failed to predict BA1, BS1, BS2, and PM2 criteria.")
        #         self.prediction.BA1.comment = ba1bs1bs2pm2_comment
        #         self.prediction.BS1.comment = ba1bs1bs2pm2_comment
        #         self.prediction.BS2.comment = ba1bs1bs2pm2_comment
        #         self.prediction.PM2.comment = ba1bs1bs2pm2_comment
        #     else:
        #         self.prediction.BA1.prediction = (
        #             ACMGPrediction.Met if ba1bs1bs2pm2_prediction.BA1 else ACMGPrediction.Unmet
        #         )
        #         self.prediction.BS1.prediction = (
        #             ACMGPrediction.Met if ba1bs1bs2pm2_prediction.BS1 else ACMGPrediction.Unmet
        #         )
        #         self.prediction.BS2.prediction = (
        #             ACMGPrediction.Met if ba1bs1bs2pm2_prediction.BS2 else ACMGPrediction.Unmet
        #         )
        #         self.prediction.PM2.prediction = (
        #             ACMGPrediction.Met if ba1bs1bs2pm2_prediction.PM2 else ACMGPrediction.Unmet
        #         )
        #         self.prediction.BA1.comment = ba1bs1bs2pm2_comment
        #         self.prediction.BS1.comment = ba1bs1bs2pm2_comment
        #         self.prediction.BS2.comment = ba1bs1bs2pm2_comment
        #         self.prediction.PM2.comment = ba1bs1bs2pm2_comment

        # except AutoAcmgBaseException as e:
        #     logger.error("Failed to predict BA1, BS1, BS2, and PM2 criteria. Error: {}", e)

        # # PM1
        # try:
        #     logger.info("Predicting PM1 criteria.")
        #     pm1 = AutoPM1(self.seqvar, variant_info.result, config=self.config)
        #     pm1_prediction, pm1_comment = pm1.predict()
        #     if not pm1_prediction:
        #         logger.error("Failed to predict PM1 criteria.")
        #         self.prediction.PM1.comment = pm1_comment
        #     else:
        #         self.prediction.PM1.prediction = (
        #             ACMGPrediction.Met if pm1_prediction.PM1 else ACMGPrediction.Unmet
        #         )
        #         self.prediction.PM1.comment = pm1_comment

        # except AutoAcmgBaseException as e:
        #     logger.error("Failed to predict PM1 criteria. Error: {}", e)

        # # PP2 and BP1
        # try:
        #     logger.info("Predicting PP2 and BP1 criteria.")
        #     pp2bp1 = AutoPP2BP1(self.seqvar, variant_info.result, config=self.config)
        #     pp2_bp1_prediction, pp2_bp1_comment = pp2bp1.predict()
        #     if not pp2_bp1_prediction:
        #         logger.error("Failed to predict PP2 and BP1 criteria.")
        #         self.prediction.PP2.comment = pp2_bp1_comment
        #     else:
        #         self.prediction.PP2.prediction = (
        #             ACMGPrediction.Met if pp2_bp1_prediction.PP2 else ACMGPrediction.Unmet
        #         )
        #         self.prediction.BP1.prediction = (
        #             ACMGPrediction.Met if pp2_bp1_prediction.BP1 else ACMGPrediction.Unmet
        #         )
        #         self.prediction.PP2.comment = pp2_bp1_comment
        #         self.prediction.BP1.comment = pp2_bp1_comment

        # except AutoAcmgBaseException as e:
        #     logger.error("Failed to predict PP2 and BP1 criteria. Error: {}", e)

        # # BP7
        # try:
        #     logger.info("Predicting BP7 criteria.")
        #     bp7 = AutoBP7(self.seqvar, variant_info.result, config=self.config)
        #     bp7_prediction, bp7_comment = bp7.predict()
        #     if not bp7_prediction:
        #         logger.error("Failed to predict BP7 criteria.")
        #         self.prediction.BP7.comment = bp7_comment
        #     else:
        #         self.prediction.BP7.prediction = (
        #             ACMGPrediction.Met if bp7_prediction.BP7 else ACMGPrediction.Unmet
        #         )
        #         self.prediction.BP7.comment = bp7_comment

        # except AutoAcmgBaseException as e:
        #     logger.error("Failed to predict BP7 criteria. Error: {}", e)

        # PP3 and BP4
        try:
            logger.info("Predicting PP3 and BP4 criteria.")
            pp3_bp4 = AutoPP3BP4(self.seqvar, variant_info.result, config=self.config)
            pp3_bp4_prediction, pp3_bp4_comment = pp3_bp4.predict()
            if not pp3_bp4_prediction:
                logger.error("Failed to predict PP3 and BP4 criteria.")
                self.prediction.PP3.comment = pp3_bp4_comment
                self.prediction.BP4.comment = pp3_bp4_comment
            else:
                self.prediction.PP3.prediction = (
                    ACMGPrediction.Met if pp3_bp4_prediction.PP3 else ACMGPrediction.Unmet
                )
                self.prediction.BP4.prediction = (
                    ACMGPrediction.Met if pp3_bp4_prediction.BP4 else ACMGPrediction.Unmet
                )
                self.prediction.PP3.comment = pp3_bp4_comment
                self.prediction.BP4.comment = pp3_bp4_comment

        except AutoAcmgBaseException as e:
            logger.error("Failed to predict PP3 and BP4 criteria. Error: {}", e)

        return self.prediction
