from typing import Optional

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.criteria.auto_ba1_bs1_bs2_pm2 import AutoBA1BS1BS2PM2
from src.criteria.auto_pm4_bp3 import AutoPM4BP3
from src.criteria.auto_ps1_pm5 import AutoPS1PM5
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.auto_acmg import ACMGCriteria
from src.defs.exceptions import AutoAcmgBaseException
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


class AutoACMGCriteria:
    """Predict ACMG criteria for sequence variant."""

    def __init__(
        self, seqvar: SeqVar, genome_release: GenomeRelease, *, config: Optional[Config] = None
    ):
        #: Configuration to use.
        self.config = config or Config()
        # Attributes to be set
        self.seqvar = seqvar
        self.genome_release = genome_release
        self.annonars_client = AnnonarsClient(api_base_url=self.config.api_base_url_annonars)
        self.prediction: Optional[ACMGCriteria] = None

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

    def predict(self) -> Optional[ACMGCriteria]:
        """Predict ACMG criteria for sequence variant."""
        self.prediction = ACMGCriteria()

        variant_info = self._get_variant_info(self.seqvar)
        if not variant_info:
            logger.error("Failed to get variant information for {}.", self.seqvar)
            return None

        # PS1 and PM5
        try:
            ps1pm5 = AutoPS1PM5(
                self.seqvar, self.genome_release, variant_info.result, config=self.config
            )
            ps1_pm5_prediction = ps1pm5.predict()
            if not ps1_pm5_prediction:
                logger.error("Failed to predict PS1&PM5 criteria.")
            else:
                self.prediction.PS1 = ps1_pm5_prediction.PS1
                self.prediction.PM5 = ps1_pm5_prediction.PM5
        except AutoAcmgBaseException as e:
            logger.error("Failed to predict PS1 and PM5 criteria. Error: {}", e)

        # PM4 and BP3
        try:
            pm4bp3 = AutoPM4BP3(
                self.seqvar, self.genome_release, variant_info.result, config=self.config
            )
            pm4_bp3_prediction = pm4bp3.predict()
            if not pm4_bp3_prediction:
                logger.error("Failed to predict PM4&BP3 criteria.")
            else:
                self.prediction.PM4 = pm4_bp3_prediction.PM4
                self.prediction.BP3 = pm4_bp3_prediction.BP3
        except AutoAcmgBaseException as e:
            logger.error("Failed to predict PM4 and BP3 criteria. Error: {}", e)

        # BA1, BS1, BS2, PM2
        try:
            ba1bs1bs2pm2 = AutoBA1BS1BS2PM2(
                self.seqvar, self.genome_release, variant_info.result, config=self.config
            )
            ba1bs1bs2pm2_prediction = ba1bs1bs2pm2.predict()
            if not ba1bs1bs2pm2_prediction:
                logger.error("Failed to predict BA1, BS1, BS2, and PM2 criteria.")
            else:
                self.prediction.BA1 = ba1bs1bs2pm2_prediction.BA1
                self.prediction.BS1 = ba1bs1bs2pm2_prediction.BS1
                self.prediction.BS2 = ba1bs1bs2pm2_prediction.BS2
                self.prediction.PM2 = ba1bs1bs2pm2_prediction.PM2
        except AutoAcmgBaseException as e:
            logger.error("Failed to predict BA1, BS1, BS2, and PM2 criteria. Error: {}", e)

        return self.prediction
