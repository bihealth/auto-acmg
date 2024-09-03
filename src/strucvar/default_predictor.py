from typing import Optional

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.auto_acmg import AutoACMGStrucVarResult
from src.defs.strucvar import StrucVar


class DefaultStrucVarPredictor:

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
        logger.info("Predicting ACMG criteria for structural variant.")
        return None
