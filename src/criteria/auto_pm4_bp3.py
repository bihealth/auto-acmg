"""Implementation of PM4 and BP3 rules for sequence variants."""

from typing import Optional, Tuple

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.annonars_variant import AnnonarsVariantResponse, VariantResult
from src.defs.auto_acmg import PM4BP3
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


class AutoPM4BP3:
    """Predicts PM4 and BP3 criteria for sequence variants."""

    def __init__(
        self,
        seqvar: SeqVar,
        genome_release: GenomeRelease,
        variant_info: VariantResult,
        *,
        config: Optional[Config] = None,
    ):
        #: Configuration to use.
        self.config = config or Config()
        #: Sequence variant to predict.
        self.seqvar = seqvar
        #: Genome release.
        self.genome_release = genome_release
        #: Variant information.
        self.variant_info = variant_info
        #: Annonars client.
        self.annonars_client = AnnonarsClient(api_base_url=self.config.api_base_url_annonars)
        #: Prediction result.
        self.prediction: PM4BP3 | None = None
        #: Comment to store the prediction explanation.
        self.comment: str = ""

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

    def _in_repeat_region(self, variant_info: VariantResult) -> bool:
        """Check if the variant is in a repeat region.

        Args:
            variant_info (VariantResult): The variant information.

        Returns:
            bool: True if the variant is in a repeat region, False otherwise.
        """
        return False

    def _in_conserved_domain(self, variant_info: VariantResult) -> bool:
        """Check if the variant is in a conserved domain.

        Args:
            variant_info (VariantResult): The variant information.

        Returns:
            bool: True if the variant is in a conserved domain, False otherwise.
        """
        return False

    def predict(self) -> Tuple[Optional[PM4BP3], str]:
        """Predicts PM4 and BP3 criteria for the provided sequence variant.

        Implementation of the rule:
        - If the variant is a stop-loss variant, PM4 is True and BP3 is False.
        - If the variant is an in-frame deletion/insertion:
        - If the variant is not in a repeat region and in a conserved domain, PM4 is True and BP3 is
        False.
        - If the variant is in a repeat region and not in a conserved domain, PM4 is False and BP3
        is True.

        Note:
            Rules:
            PM4: Protein length changes due to in-frame deletions/insertions in a non-repeat region
            or stop-loss variants.
            BP3: In-frame deletions/insertions in a repetitive region without a known function.

        Returns:
            PM4BP3: PM4 and BP3 prediction.
        """
        try:
            # Initialize the prediction result
            self.prediction = PM4BP3()

            if not self.variant_info:
                raise MissingDataError("No valid primary information.")

            self.comment = "Predicting???"
            # Stop-loss variants are considered as PM4
            if self.variant_info.cadd and self.variant_info.cadd.ConsDetail == "stop_loss":
                self.comment += f"Variant consequence is stop-loss. PM4 is met"
                self.prediction.PM4 = True
                self.prediction.BP3 = False
            elif self.variant_info.cadd and self.variant_info.cadd.ConsDetail in [
                "inframe_deletion",
                "inframe_insertion",
            ]:
                self.comment += f"Variant consequence is in-frame deletion/insertion."
                if not self._in_repeat_region(self.variant_info) and self._in_conserved_domain(
                    self.variant_info
                ):
                    self.comment += (
                        f"Variant is not in a repeat region and in a conserved domain."
                        f"PM4 is met."
                    )
                    self.prediction.PM4 = True
                    self.prediction.BP3 = False
                elif self._in_repeat_region(self.variant_info) and not self._in_conserved_domain(
                    self.variant_info
                ):
                    self.comment += (
                        f"Variant is in a repeat region and not in a conserved domain."
                        f"BP3 is met."
                    )
                    self.prediction.PM4 = False
                    self.prediction.BP3 = True

        except AutoAcmgBaseException as e:
            logger.error("Failed to predict PM4 and BP3 criteria. Error: {}", e)
            self.comment += f"An error occured while predicting PM4 and BP3 criteria: {e}"
            self.prediction = None

        # Return the prediction and the explanation
        return self.prediction, self.comment
