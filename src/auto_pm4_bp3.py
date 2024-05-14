"""Implementation of PM4 and BP3 rules for sequence variants."""

from typing import Optional

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.auto_acmg import PM4BP3
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


class AutoPM4BP3:
    """Predicts PM4 and BP3 criteria for sequence variants."""

    def __init__(self, seqvar: SeqVar, genome_release: GenomeRelease):
        self.seqvar = seqvar
        self.genome_release = genome_release
        self.annonars_client = AnnonarsClient()
        self.prediction: PM4BP3 | None = None

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

    def _in_repeat_region(self, variant_info: AnnonarsVariantResponse) -> bool:
        """Check if the variant is in a repeat region.

        Args:
            variant_info (AnnonarsVariantResponse): The variant information.

        Returns:
            bool: True if the variant is in a repeat region, False otherwise.
        """
        return False

    def _in_conserved_domain(self, variant_info: AnnonarsVariantResponse) -> bool:
        """Check if the variant is in a conserved domain.

        Args:
            variant_info (AnnonarsVariantResponse): The variant information.

        Returns:
            bool: True if the variant is in a conserved domain, False otherwise.
        """
        return False

    def predict(self) -> Optional[PM4BP3]:
        """Predicts PM4 and BP3 criteria for the provided sequence variant.

        Note:
            Rule:
            PM4: Protein length changes due to in-frame deletions/insertions in a non-repeat region or stop-loss
            variants.
            BP3: In-frame deletions/insertions in a repetitive region without a known function.
            Implementation:
            - If the variant is a stop-loss variant, PM4 is True and BP3 is False.
            - If the variant is an in-frame deletion/insertion:
            - If the variant is not in a repeat region and in a conserved domain, PM4 is True and BP3 is False.
            - If the variant is in a repeat region and not in a conserved domain, PM4 is False and BP3 is True.

        Returns:
            PM4BP3: PM4 and BP3 prediction.
        """
        try:
            # Initialize the prediction result
            self.prediction = PM4BP3()

            primary_info = self._get_variant_info(self.seqvar)
            if not primary_info:
                raise AlgorithmError("No valid primary information.")

            # Stop-loss variants are considered as PM4
            if primary_info.result.cadd and primary_info.result.cadd.ConsDetail == "stop_loss":
                self.prediction.PM4 = True
                self.prediction.BP3 = False
            elif primary_info.result.cadd and primary_info.result.cadd.ConsDetail in [
                "inframe_deletion",
                "inframe_insertion",
            ]:
                if not self._in_repeat_region(primary_info) and self._in_conserved_domain(primary_info):
                    self.prediction.PM4 = True
                    self.prediction.BP3 = False
                elif self._in_repeat_region(primary_info) and not self._in_conserved_domain(primary_info):
                    self.prediction.PM4 = False
                    self.prediction.BP3 = True

        except AutoAcmgBaseException as e:
            logger.error("Failed to predict PM4 and BP3 criteria. Error: {}", e)
            self.prediction = None

        # Return the prediction
        return self.prediction
