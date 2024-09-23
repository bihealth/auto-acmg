from typing import Optional

from loguru import logger

from src.api.reev.annonars import AnnonarsClient
from src.core.config import settings
from src.defs.auto_acmg import AutoACMGStrucVarResult
from src.defs.strucvar import StrucVar
from src.strucvar.auto_pvs1 import AutoPVS1


class DefaultStrucVarPredictor(AutoPVS1):
    def __init__(self, strucvar: StrucVar, result: AutoACMGStrucVarResult):
        #: Structural variant to predict.
        self.strucvar = strucvar
        #: Annonars client.
        self.annonars_client: AnnonarsClient = AnnonarsClient(
            api_base_url=settings.AUTO_ACMG_API_ANNONARS_URL
        )
        #: Prediction result.
        self.result = result

    def predict(self) -> Optional[AutoACMGStrucVarResult]:
        """Predict ACMG criteria for the structural variant."""
        # Currently only PVS1 prediction
        logger.warning("Currently only PVS1 prediction is implemented.")

        # PVS1
        self.result.criteria.pvs1 = self.predict_pvs1(self.strucvar, self.result.data)

        logger.info("StrucVar prediction finished.")
        return self.result
