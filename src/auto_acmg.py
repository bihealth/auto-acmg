"""Implementations of the PVS1 algorithm."""

from typing import List, Optional

from loguru import logger

from src.core.config import Config
from src.criteria.auto_acmg import AutoACMGCriteria
from src.defs.auto_acmg import ACMGPrediction, AutoACMGResult, CriteriaPrediction
from src.defs.auto_pvs1 import PVS1Prediction, PVS1PredictionPathMapping, PVS1PredictionSeqVarPath
from src.defs.exceptions import AutoAcmgBaseException, ParseError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar, SeqVarResolver
from src.defs.strucvar import StrucVar, StrucVarResolver
from src.pvs1.auto_pvs1 import AutoPVS1

#: Pathogenic PVS1 predictions of sequence variants.
PVS1_POSITIVE_SEQVAR_PREDICTIONS = [
    PVS1Prediction.PVS1,
    PVS1Prediction.PVS1_Strong,
    PVS1Prediction.PVS1_Moderate,
    PVS1Prediction.PVS1_Supporting,
]

#: Criteria, which have no automatic prediction yet.
NOT_IMPLEMENTED_CRITERIA: List[str] = [
    "PS2",  # De novo (both parents not tested)
    "PS3",  # Well-established in vitro or in vivo functional studies
    "PS4",  # Increased prevalence in affected individuals vs. controls
    "PM3",  # Recessive disorder (zygosity unknown)
    "PM6",  # Assumed de novo (both parents not tested)
    "PP1",  # Cosegregation with disease in multiple affected family members
    "PP4",  # Patient's phenotype or family history specificity
    "BS3",  # Well-established in vitro or in vivo functional studies
    "BS4",  # Lack of segregation in affected family members
    "BP2",  # Observed in healthy individuals (zygosity unknown)
    "BP5",  # Case with an alternate molecular basis for disease
]


class AutoACMG:
    """Class for predicting ACMG criteria.

    This class handles both sequence variants and structural variants to determine their potential
    impact under the various criteria of the ACMG guidelines for variant classification. Currently
    it only implements the PVS1 criterion (not finished yet).

    Attributes:
        variant_name (str): The name or identifier of the variant being analyzed.
        genome_release (GenomeRelease): The genome release version, defaults to GRCh38.
    """

    def __init__(
        self,
        variant_name: str,
        genome_release: GenomeRelease = GenomeRelease.GRCh38,
        *,
        config: Optional[Config] = None,
    ):
        """Initializes the AutoACMG with the specified variant and genome release.

        Args:
            variant_name: The name or identifier of the variant.
            genome_release (Optional): The genome release version, such as GRCh38 or GRCh37.
        """
        #: Configuration to use.
        self.config = config or Config()
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
            seqvar_resolver = SeqVarResolver(config=self.config)
            seqvar: SeqVar = seqvar_resolver.resolve_seqvar(self.variant_name, self.genome_release)
            logger.debug("Resolved sequence variant: {}", seqvar)
            return seqvar
        except ParseError:
            logger.exception("Failed to resolve sequence variant, trying structural variant.")
            try:
                strucvar_resolver = StrucVarResolver(config=self.config)
                strucvar: StrucVar = strucvar_resolver.resolve_strucvar(
                    self.variant_name, self.genome_release
                )
                logger.debug("Resolved structural variant: {}", strucvar)
                return strucvar
            except ParseError as e:
                logger.error("Failed to resolve structural variant: {}", e)
                return None
        except AutoAcmgBaseException as e:
            logger.error("An unexpected error occurred: {}", e)
            return None

    def predict(self) -> Optional[AutoACMGResult]:
        """
        Predict ACMG criteria for the specified variant.

        This method first resolves the variant, then predicts the PVS1 criterion for sequence
        variants and other ACMG criteria.

        Note:
            The method can resolve both sequence and structural variants, but currently only
            sequence variants are supported for ACMG criteria prediction.

        Raises:
            Exception: Specific exceptions are caught and logged, but generic exceptions may be
            raised if the prediction fails.
        """
        logger.info("Predicting ACMG criteria for variant: {}", self.variant_name)
        variant = self.resolve_variant()
        if not variant:
            logger.error("Failed to resolve variant: {}", self.variant_name)
            return None

        if isinstance(variant, SeqVar):
            logger.info(
                "Classifying ACMG criteria for sequence variant {}, genome release: {}.",
                variant.user_repr,
                self.genome_release.name,
            )
            self.seqvar: SeqVar = variant
            self.prediction: AutoACMGResult = AutoACMGResult()

            # PP5 and BP6 criteria
            self.prediction.pp5.prediction = ACMGPrediction.NotSet
            self.prediction.pp5.comment = "PP5 prediction is deprecated."
            self.prediction.bp6.prediction = ACMGPrediction.NotSet
            self.prediction.bp6.comment = "BP6 prediction is deprecated."
            logger.warning("Note, that PP5 and BP6 criteria are depricated and not predicted.")

            # Not implemented criteria
            for crit in NOT_IMPLEMENTED_CRITERIA:
                setattr(
                    getattr(self.prediction, crit.lower()),
                    "prediction",
                    ACMGPrediction.NotSet,
                )
                setattr(
                    getattr(self.prediction, crit.lower()),
                    "comment",
                    f"{crit} prediction is not implemented.",
                )
            logger.warning(
                "Some criteria are not implemented yet: {}",
                NOT_IMPLEMENTED_CRITERIA,
            )

            # PVS1
            try:
                logger.info("Predicting PVS1.")
                pvs1 = AutoPVS1(self.seqvar, self.genome_release, config=self.config)
                seqvar_prediction, seqvar_prediction_path = pvs1.predict()
                if seqvar_prediction is None or seqvar_prediction_path is None:
                    raise AutoAcmgBaseException(
                        "PVS1 prediction failed: prediction or prediction path is None."
                    )
                else:
                    if seqvar_prediction == PVS1Prediction.NotSet:
                        raise AutoAcmgBaseException("PVS1 prediction failed: prediction NotSet.")
                    self.prediction.pvs1.prediction = (
                        ACMGPrediction.Positive
                        if seqvar_prediction in PVS1_POSITIVE_SEQVAR_PREDICTIONS
                        else ACMGPrediction.Negative
                    )
                    self.prediction.pvs1.comment = (
                        f"PVS1 strength: {seqvar_prediction.name}. "
                        f"PVS1 prediction path: {PVS1PredictionPathMapping[seqvar_prediction_path]}."
                    )
            except AutoAcmgBaseException as e:
                self.prediction.pvs1.prediction = ACMGPrediction.NotSet
                self.prediction.pvs1.comment = "PVS1 prediction failed."
                logger.error("Failed to predict PVS1 criteria. Error: {}", e)

            # Other criteria
            try:
                logger.info("Predicting other ACMG criteria.")
                auto_criteria = AutoACMGCriteria(
                    self.seqvar, self.genome_release, config=self.config
                )
                criteria_preciction = auto_criteria.predict()
                if criteria_preciction is None:
                    raise AutoAcmgBaseException("Other ACMG criteria prediction failed.")
                else:
                    for criteria, prediction in criteria_preciction.model_dump().items():
                        setattr(
                            getattr(self.prediction, criteria.lower()),
                            "prediction",
                            ACMGPrediction.Positive if prediction else ACMGPrediction.Negative,
                        )
            except AutoAcmgBaseException as e:
                logger.error("Failed to predict other ACMG criteria. Error: {}", e)

            logger.info("ACMG criteria prediction completed.")
            return self.prediction

        elif isinstance(variant, StrucVar):
            logger.info("Classification of structural variants is not implemented yet.")
            # logger.info("Currently only PVS1 prediction is implemented for structural variants!")
            # logger.info(
            #     "Classifying ACMG criteria for structural variant {}, genome release: {}.",
            #     variant.user_repr,
            #     self.genome_release.name,
            # )
            # # PVS1
            # try:
            #     logger.info("Predicting PVS1.")
            #     self.strucvar: StrucVar = variant
            #     self.strucvar_prediction: PVS1Prediction = PVS1Prediction.NotSet  # type: ignore
            #     self.strucvar_prediction_path: PVS1PredictionStrucVarPath = PVS1PredictionStrucVarPath.NotSet  # type: ignore

            #     pvs1 = AutoPVS1(self.strucvar, self.genome_release)
            #     strucvar_prediction, strucvar_prediction_path = pvs1.predict()
            #     if strucvar_prediction is None or strucvar_prediction_path is None:
            #         logger.error("Failed to predict PVS1 criteria.")
            #         return
            #     else:
            #         # Double check if the prediction path is indeed for structural variant
            #         assert isinstance(strucvar_prediction_path, PVS1PredictionStrucVarPath)
            #         self.strucvar_prediction = strucvar_prediction
            #         self.strucvar_prediction_path = strucvar_prediction_path
            #         logger.info(
            #             "PVS1 prediction for {}: {}.\n" "The prediction path is:\n{}.",
            #             self.strucvar.user_repr,
            #             self.strucvar_prediction.name,
            #             PVS1PredictionPathMapping[self.strucvar_prediction_path],
            #         )
            # except AutoAcmgBaseException as e:
            #     logger.error("Failed to predict PVS1 criteria. Error: {}", e)
            #     return

        return None
