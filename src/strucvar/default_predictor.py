from typing import Optional

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.auto_acmg import AutoACMGStrucVarResult
from src.defs.strucvar import StrucVar
from src.strucvar.auto_pvs1 import AutoPVS1


class DefaultStrucVarPredictor(AutoPVS1):
    def __init__(self, strucvar: StrucVar, result: AutoACMGStrucVarResult, config: Config):
        #: Configuration to use.
        self.config = config or Config()
        #: Structural variant to predict.
        self.strucvar = strucvar
        #: Annonars client.
        self.annonars_client: AnnonarsClient = AnnonarsClient(
            api_base_url=self.config.api_base_url_annonars
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
