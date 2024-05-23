"""Implementation of BA1, BS1, BS2, PS4, PM2 prediction for sequence variants."""

import re
from typing import Optional

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.annonars_variant import AnnonarsVariantResponse, VariantResult
from src.defs.auto_acmg import BA1BS1BS2PS4PM2
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


class AutoBA1BS1BS2PS4PM2:
    """Predicts BA1, BS1, BS2, PS4, PM2 criteria for sequence variants."""

    def __init__(
        self, seqvar: SeqVar, genome_release: GenomeRelease, *, config: Optional[Config] = None
    ):
        #: Configuration to use.
        self.config = config or Config()
        # Attributes to be set
        self.seqvar = seqvar
        self.genome_release = genome_release
        self.annonars_client = AnnonarsClient(api_base_url=self.config.api_base_url_annonars)
        self.prediction: BA1BS1BS2PS4PM2 | None = None

    def _get_variant_info(self, seqvar: SeqVar) -> Optional[AnnonarsVariantResponse]:
        """Get variant information from Annonars.

        Returns:
            dict: Annonars response.
        """
        try:
            logger.debug("Getting variant information for {}.", seqvar)
            return self.annonars_client.get_variant_info(seqvar)
        except AutoAcmgBaseException as e:
            logger.error("Failed to get variant information. Error: {}", e)
            return None

    def predict(self) -> Optional[BA1BS1BS2PS4PM2]:
        """
        Predicts the BA1, BS1, BS2, PS4, PM5 criteria for the sequence variant.

        Returns:
            BA1BS1BS2PS4PM2: The prediction result.
        """

        # Return the prediction result
        return self.prediction
