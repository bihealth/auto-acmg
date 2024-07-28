"""Implementation of PM4 and BP3 rules for sequence variants."""

import os
from typing import Optional, Tuple

import tabix  # type: ignore
from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config, settings
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
        variant_info: VariantResult,
        *,
        config: Optional[Config] = None,
    ):
        #: Configuration to use.
        self.config: Config = config or Config()
        #: Sequence variant to predict.
        self.seqvar: SeqVar = seqvar
        #: Variant information.
        self.variant_info: VariantResult = variant_info
        #: Prediction result.
        self.prediction: Optional[PM4BP3] = None
        #: Comment to store the prediction explanation.
        self.comment: str = ""

    def _in_repeat_region(self, seqvar: SeqVar) -> bool:
        """Check if the variant is in a repeat region.

        Check if the variant is in a repeat region using the RepeatMasker track.

        Args:
            seqvar (SeqVar): Sequence variant.

        Returns:
            bool: True if the variant is in a repeat region, False otherwise.
        """
        self.comment += "Check if the variant is in a repeat region.\n"
        try:
            # Find path to the lib file
            if seqvar.genome_release == GenomeRelease.GRCh37:
                path = os.path.join(settings.PATH_TO_ROOT, "lib", "rmsk", "grch37", "rmsk.bed.gz")
            else:
                path = os.path.join(settings.PATH_TO_ROOT, "lib", "rmsk", "grch38", "rmsk.bed.gz")
            tb = tabix.open(path)
            records = tb.query(f"chr{seqvar.chrom}", seqvar.pos - 1, seqvar.pos)
            # Check if iterator is not empty
            if any(True for _ in records):
                self.comment += "Variant is in a repeat region.\n"
                return True
            else:
                self.comment += "Variant is not in a repeat region.\n"
                return False
        except tabix.TabixError as e:
            logger.error("Failed to check if the variant is in a repeat region. Error: {}", e)
            raise AlgorithmError("Failed to check if the variant is in a repeat region.") from e

    def predict(self) -> Tuple[Optional[PM4BP3], str]:
        """Predicts PM4 and BP3 criteria for the provided sequence variant.

        Implementation of the rule:
        - If the variant is a stop-loss variant, PM4 is True and BP3 is False.
        - If the variant is an in-frame deletion/insertion:
        - If the variant is not in a repeat region, PM4 is True and BP3 is False.
        - If the variant is in a repeat region, PM4 is False and BP3 is True.
        - Otherwise, PM4 and BP3 are False.

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
            self.comment = "Check consequences of the variant for PM4.\n"
            # Stop-loss variants are considered as PM4
            if self.variant_info.cadd and self.variant_info.cadd.ConsDetail == "stop_loss":
                self.comment += f"Variant consequence is stop-loss. PM4 is met."
                self.prediction.PM4 = True
            # In-frame deletions/insertions
            elif self.variant_info.cadd and self.variant_info.cadd.ConsDetail in [
                "inframe_deletion",
                "inframe_insertion",
            ]:
                self.comment += f"Variant consequence is in-frame deletion/insertion.\n"
                if not self._in_repeat_region(self.seqvar):
                    self.comment += (
                        "Variant is not in a repeat region or a conserved domain. PM4 is met."
                    )
                    self.prediction.PM4 = True
                else:
                    self.comment += (
                        "Variant is in a repeat region or not in a conserved domain. BP3 is met."
                    )
                    self.prediction.BP3 = True
            else:
                self.comment += (
                    "Variant consequence is not indel or stop-loss. PM4 and BP3 are not met."
                )
                self._in_repeat_region(self.seqvar)
                self.prediction.PM4 = False
                self.prediction.BP3 = False

        except AutoAcmgBaseException as e:
            logger.error("Failed to predict PM4 and BP3 criteria. Error: {}", e)
            self.comment += f"An error occured while predicting PM4 and BP3 criteria: {e}"
            self.prediction = None

        # Return the prediction and the explanation
        return self.prediction, self.comment
