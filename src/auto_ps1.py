"""Implementation of PS1 prediction for sequence variants."""

from typing import Optional

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


class AutoPS1:
    """Predicts PS1 criteria for sequence variants."""

    def __init__(self, seqvar: SeqVar, genome_release: GenomeRelease):
        self.seqvar = seqvar
        self.genome_release = genome_release

    def get_variant_info(self, seqvar: SeqVar) -> Optional[AnnonarsVariantResponse]:
        """Get variant information from Annonars.

        Returns:
            dict: Annonars response.
        """
        try:
            logger.debug("Getting variant information for {}.", seqvar)
            annonars_client = AnnonarsClient()
            return annonars_client.get_variant_info(seqvar)
        except Exception as e:
            logger.error("Failed to get variant information. Error: {}", e)
            return None

    def predict(self) -> Optional[bool]:
        """Predicts PS1 criteria for sequence variants.

        Returns:
            bool: True if the variant meets PS1 criteria, False if it does not, or None if
            prediction failed.
        """
        try:
            # Predict PS1 criteria
            variant_info = self.get_variant_info(self.seqvar)
            if variant_info:
                return True
            else:
                logger.error("Failed to predict PS1 criteria.")
            return None
        except Exception as e:
            logger.error("Failed to predict PS1 criteria. Error: {}", e)
            return None
