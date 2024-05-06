"""Implementations of the PVS1 algorithm."""

from typing import Optional

from loguru import logger

from src.auto_ps1 import AutoPS1
from src.defs.autopvs1 import (
    PVS1Prediction,
    PVS1PredictionPathMapping,
    PVS1PredictionSeqVarPath,
    PVS1PredictionStrucVarPath,
)
from src.defs.exceptions import ParseError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar, SeqVarResolver
from src.defs.strucvar import StrucVar, StrucVarResolver
from src.pvs1.auto_pvs1 import AutoPVS1


class AutoACMG:
    """Class for predicting ACMG criteria.

    This class handles both sequence variants and structural variants to determine their potential
    impact under the various criteria of the ACMG guidelines for variant classification. Currently
    it only implements the PVS1 criterion (not finished yet).

    Attributes:
        variant_name (str): The name or identifier of the variant being analyzed.
        genome_release (GenomeRelease): The genome release version, defaults to GRCh38.
    """

    def __init__(self, variant_name: str, genome_release: GenomeRelease = GenomeRelease.GRCh38):
        """Initializes the AutoACMG with the specified variant and genome release.

        Args:
            variant_name: The name or identifier of the variant.
            genome_release (Optional): The genome release version, such as GRCh38 or GRCh37.
        """
        self.variant_name = variant_name
        self.genome_release = genome_release
        logger.debug(
            "AutoACMG initialized with variant: {} and genome release: {}",
            variant_name,
            genome_release,
        )

    def resolve_variant(self) -> SeqVar | StrucVar | None:
        """Attempts to resolve the specified variant as either a sequence or structural variant.

        This method first tries to resolve the variant as a sequence variant. If it fails, it then
        attempts to resolve it as a structural variant.

        Returns:
            SeqVar, StrucVar, or None: The resolved variant object or None if resolution fails.

        Raises:
            Exception: Specific exceptions are caught and logged, but generic exceptions may be
            raised if both resolutions fail.
        """
        logger.debug("Resolving variant: {}", self.variant_name)
        try:
            seqvar_resolver = SeqVarResolver()
            seqvar: SeqVar = seqvar_resolver.resolve_seqvar(self.variant_name, self.genome_release)
            logger.debug("Resolved sequence variant: {}", seqvar)
            return seqvar
        except ParseError:
            logger.exception("Failed to resolve sequence variant, trying structural variant.")
            try:
                strucvar_resolver = StrucVarResolver()
                strucvar: StrucVar = strucvar_resolver.resolve_strucvar(
                    self.variant_name, self.genome_release
                )
                logger.debug("Resolved structural variant: {}", strucvar)
                return strucvar
            except ParseError as e:
                logger.error("Failed to resolve structural variant: {}", e)
                return None
        except Exception as e:
            logger.error("An unexpected error occurred: {}", e)
            return None

    def predict(self):
        """Runs the prediction algorithm to assess the PVS1 criteria for the resolved variant.

        This method resolves the variant and then, based on the type of variant, predicts its
        classification according to the PVS1 criteria. It handles both sequence and structural variants.

        Raises:
            Exception: Handles general exceptions that may occur during prediction and logs them.
        """
        logger.info("Predicting ACMG criteria for variant: {}", self.variant_name)
        variant = self.resolve_variant()
        if not variant:
            logger.error("Failed to resolve variant: {}", self.variant_name)
            return

        if isinstance(variant, SeqVar):
            logger.info(
                "Classifying ACMG criteria for sequence variant {}, genome release: {}.",
                variant.user_repr,
                self.genome_release.name,
            )
            self.seqvar: SeqVar = variant
            self.seqvar_pvs1_prediction: PVS1Prediction = PVS1Prediction.NotSet
            self.seqvar_pvs1_prediction_path: PVS1PredictionSeqVarPath = (
                PVS1PredictionSeqVarPath.NotSet
            )
            self.seqvar_ps1: Optional[bool] = None

            # PVS1
            try:
                logger.info("Predicting PVS1.")
                pvs1 = AutoPVS1(self.seqvar, self.genome_release)
                seqvar_prediction, seqvar_prediction_path = pvs1.predict()
                if seqvar_prediction is None or seqvar_prediction_path is None:
                    logger.error("Failed to predict PVS1 criteria.")
                else:
                    # Double check if the prediction path is indeed for sequence variant
                    assert isinstance(seqvar_prediction_path, PVS1PredictionSeqVarPath)
                    self.seqvar_pvs1_prediction = seqvar_prediction
                    self.seqvar_pvs1_prediction_path = seqvar_prediction_path
                    logger.info(
                        "PVS1 prediction for {}: {}.\n" "The prediction path is:\n{}.",
                        self.seqvar.user_repr,
                        self.seqvar_pvs1_prediction.name,
                        PVS1PredictionPathMapping[self.seqvar_pvs1_prediction_path],
                    )
            except Exception as e:
                logger.error("Failed to predict PVS1 criteria. Error: {}", e)

            # PS1
            try:
                logger.info("Predicting PS1.")
                ps1 = AutoPS1(self.seqvar, self.genome_release)
                self.seqvar_ps1 = ps1.predict()
                if self.seqvar_ps1 is None:
                    logger.error("Failed to predict PS1 criteria.")
                else:
                    logger.info(
                        "PS1 prediction for {}: {}.", self.seqvar.user_repr, self.seqvar_ps1
                    )
            except Exception as e:
                logger.error("Failed to predict PS1 criteria. Error: {}", e)

        elif isinstance(variant, StrucVar):
            logger.info(
                "Classifying ACMG criteria for structural variant {}, genome release: {}.",
                variant.user_repr,
                self.genome_release.name,
            )
            # PVS1
            try:
                logger.info("Predicting PVS1.")
                self.strucvar: StrucVar = variant
                self.strucvar_prediction: PVS1Prediction = PVS1Prediction.NotSet  # type: ignore
                self.strucvar_prediction_path: PVS1PredictionStrucVarPath = PVS1PredictionStrucVarPath.NotSet  # type: ignore

                pvs1 = AutoPVS1(self.strucvar, self.genome_release)
                strucvar_prediction, strucvar_prediction_path = pvs1.predict()
                if strucvar_prediction is None or strucvar_prediction_path is None:
                    logger.error("Failed to predict PVS1 criteria.")
                    return
                else:
                    # Double check if the prediction path is indeed for structural variant
                    assert isinstance(strucvar_prediction_path, PVS1PredictionStrucVarPath)
                    self.strucvar_prediction = strucvar_prediction
                    self.strucvar_prediction_path = strucvar_prediction_path
                    logger.info(
                        "PVS1 prediction for {}: {}.\n" "The prediction path is:\n{}.",
                        self.strucvar.user_repr,
                        self.strucvar_prediction.name,
                        PVS1PredictionPathMapping[self.strucvar_prediction_path],
                    )
            except Exception as e:
                logger.error("Failed to predict PVS1 criteria. Error: {}", e)
                return
        logger.info("ACMG criteria prediction completed.")
