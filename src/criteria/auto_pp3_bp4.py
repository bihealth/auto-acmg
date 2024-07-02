"""Implementation of the PP3 and BP4 criteria."""

from typing import List, Optional, Tuple, Union

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.annonars_variant import VariantResult
from src.defs.auto_acmg import PP3BP4, MissenseScores
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


class AutoPP3BP4:
    """Class for automatic PP3 and BP4 prediction."""

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
        #: Annonars client.
        self.annonars_client: AnnonarsClient = AnnonarsClient(
            api_base_url=self.config.api_base_url_annonars
        )
        #: Prediction result.
        self.prediction: Optional[PP3BP4] = None
        #: Comment to store the prediction explanation.
        self.comment: str = ""

    def _convert_score_val(self, score_value: Optional[Union[str, float, int]]) -> Optional[float]:
        """
        Convert score value to float.

        Since the score values can be represented as strings (with ";" as separator), we pick the
        maximum value that is not empty ("."). If the value is already numeric, we return it as is.

        Args:
            score_value (Optional[Union[str, float, int]]): Score value to convert.

        Returns:
            Optional[float]: Converted score value.

        Raises:
            AlgorithmError: If the score value cannot be converted to float.
        """
        if score_value is None:
            return None
        if isinstance(score_value, (float, int)):
            return float(score_value)
        try:
            return max(float(score) for score in score_value.split(";") if score != ".")
        except ValueError as e:
            logger.error("Failed to convert score value to float. Error: {}", e)
            raise AlgorithmError("Failed to convert score value to float.") from e

    def _is_pathogenic_score(self, variant_info: VariantResult) -> bool:
        """
        Check if any of the pathogenic scores meet the threshold.

        Go through the Missense scores and check if any of the pathogenic scores meet the threshold.

        Args:
            variant_info (VariantResult): Variant information.

        Returns:
            bool: True if the variant is pathogenic, False otherwise.

        Raises:
            MissingDataError: If the variant information is missing.
        """
        if not variant_info.dbnsfp:
            logger.error("Missing dbNSFP data.")
            raise MissingDataError("Missing dbNSFP data.")
        self.comment += "Checking for pathogenic scores: \n"
        logger.debug("Checking for pathogenic scores.")
        for score in MissenseScores:
            self.comment += f"Checking {score.name} score: "
            logger.debug("Checking {} score.", score.name)
            score_value = self._convert_score_val(getattr(variant_info.dbnsfp, score.name, None))
            if score_value is not None:
                if score_value >= score.pathogenic_threshold:
                    self.comment += f"{score_value} >= {score.pathogenic_threshold}. =>\n"
                    logger.debug(
                        "Pathogenic score({}): {} >= {}",
                        score.name,
                        score_value,
                        score.pathogenic_threshold,
                    )
                    return True
                else:
                    self.comment += f"{score_value} < {score.pathogenic_threshold}.\n"
                    logger.debug(
                        "Pathogenic score({}): {} < {}",
                        score.name,
                        score_value,
                        score.pathogenic_threshold,
                    )
            else:
                self.comment += "Score not found. \n"
                logger.debug("Score not found.")
        return False

    def _is_benign_score(self, variant_info: VariantResult) -> bool:
        """
        Check if any of the benign scores meet the threshold.

        Go through the Missense scores and check if any of the benign scores meet the threshold.

        Args:
            variant_info (VariantResult): Variant information.

        Returns:
            bool: True if the variant is benign, False otherwise.

        Raises:
            MissingDataError: If the variant information is missing.
        """
        if not variant_info.dbnsfp:
            logger.error("Missing dbNSFP data.")
            raise MissingDataError("Missing dbNSFP data.")
        self.comment += "Checking for benign scores: \n"
        logger.debug("Checking for benign scores.")
        for score in MissenseScores:
            self.comment += f"Checking {score.name} score: "
            logger.debug("Checking {} score.", score.name)
            score_value = self._convert_score_val(getattr(variant_info.dbnsfp, score.name, None))
            if score_value is not None and score.benign_threshold is not None:
                if score_value <= score.benign_threshold:
                    self.comment += f"{score_value} <= {score.benign_threshold}. =>\n"
                    logger.debug(
                        "Benign score({}): {} <= {}",
                        score.name,
                        score_value,
                        score.benign_threshold,
                    )
                    return True
                else:
                    self.comment += f"{score_value} > {score.benign_threshold}.\n"
                    logger.debug(
                        "Benign score({}): {} > {}", score.name, score_value, score.benign_threshold
                    )
            else:
                self.comment += "Score not found. \n"
                logger.debug("Score not found.")
        return False

    def _is_pathogenic_spliceai(self, variant_info: VariantResult) -> bool:
        """
        Check if any of the pathogenic scores meet the threshold.

        Note:
            The threshold is set to 0.2 according to doi:10.1101/2023.02.24.23286431.

        Args:
            variant_info (VariantResult): Variant information.

        Returns:
            bool: True if the variant is pathogenic, False otherwise.

        Raises:
            MissingDataError: If the variant information is missing.
        """
        if (
            not variant_info.gnomad_exomes
            or not variant_info.gnomad_exomes.effectInfo
            or not variant_info.gnomad_exomes.effectInfo.spliceaiDsMax
        ):
            logger.error("Missing GnomAD exomes data.")
            raise MissingDataError("Missing GnomAD exomes data.")
        self.comment += "Checking for pathogenic SpliceAI score: \n"
        logger.debug("Checking for pathogenic SpliceAI score.")
        if variant_info.gnomad_exomes.effectInfo.spliceaiDsMax >= 0.2:
            self.comment += (
                "Pathogenic SpliceAI score: "
                f"{variant_info.gnomad_exomes.effectInfo.spliceaiDsMax} >= 0.2. =>\n"
            )
            logger.debug(
                "Pathogenic SpliceAI score: {}", variant_info.gnomad_exomes.effectInfo.spliceaiDsMax
            )
            return True
        else:
            self.comment += (
                "Benign SpliceAI score: "
                f"{variant_info.gnomad_exomes.effectInfo.spliceaiDsMax} < 0.2.\n"
            )
            logger.debug(
                "Benign SpliceAI score: {}", variant_info.gnomad_exomes.effectInfo.spliceaiDsMax
            )
            return False

    def _is_benign_spliceai(self, variant_info: VariantResult) -> bool:
        """
        Check if any of the pathogenic scores meet the threshold.

        Note:
            The threshold is set to 0.1 according to doi:10.1101/2023.02.24.23286431.

        Args:
            variant_info (VariantResult): Variant information.

        Returns:
            bool: True if the variant is benign, False otherwise.

        Raises:
            MissingDataError: If the variant information is missing.
        """
        if (
            not variant_info.gnomad_exomes
            or not variant_info.gnomad_exomes.effectInfo
            or not variant_info.gnomad_exomes.effectInfo.spliceaiDsMax
        ):
            logger.error("Missing GnomAD exomes data.")
            raise MissingDataError("Missing GnomAD exomes data.")
        self.comment += "Checking for benign SpliceAI score: \n"
        logger.debug("Checking for benign SpliceAI score.")
        if variant_info.gnomad_exomes.effectInfo.spliceaiDsMax <= 0.1:
            self.comment += (
                "Benign SpliceAI score: "
                f"{variant_info.gnomad_exomes.effectInfo.spliceaiDsMax} <= 0.1. =>\n"
            )
            logger.debug(
                "Benign SpliceAI score: {}", variant_info.gnomad_exomes.effectInfo.spliceaiDsMax
            )
            return True
        else:
            self.comment += (
                "Pathogenic SpliceAI score: "
                f"{variant_info.gnomad_exomes.effectInfo.spliceaiDsMax} > 0.1.\n"
            )
            logger.debug(
                "Pathogenic SpliceAI score: {}", variant_info.gnomad_exomes.effectInfo.spliceaiDsMax
            )
            return False

    def predict(self) -> Tuple[Optional[PP3BP4], str]:
        """Predict PP3 and BP4 criteria."""
        self.prediction = PP3BP4()

        try:
            if self.seqvar.chrom == "MT":
                self.comment = "Variant is in mitochondrial DNA. PP3 and BP4 criteria are not met."
                logger.debug("Variant is in mitochondrial DNA. PP3 and BP4 criteria are not met.")
                self.prediction.PP3 = False
                self.prediction.BP4 = False
                return self.prediction, self.comment

            if not self.variant_info:
                logger.error("Missing variant data.")
                raise MissingDataError("Missing variant data.")

            # Evaluate PP3 and BP4 criteria
            self.comment = "Checking Scores. => \n"
            logger.debug("Checking Scores.")
            is_pathogenic = self._is_pathogenic_score(
                self.variant_info
            ) or self._is_pathogenic_spliceai(self.variant_info)
            is_benign = self._is_benign_score(self.variant_info) or self._is_benign_spliceai(
                self.variant_info
            )
            self.comment += (
                f"Result: PP3 is {'met' if is_pathogenic else 'not met'}, "
                f"BP4 is {'met' if is_benign else 'not met'}. => \n"
            )
            logger.debug(
                "Result: PP3 is {}, BP4 is {}.",
                "met" if is_pathogenic else "not met",
                "met" if is_benign else "not met",
            )
            self.prediction.PP3 = is_pathogenic
            self.prediction.BP4 = is_benign

        except AutoAcmgBaseException as e:
            self.comment += f"An error occurred during prediction. Error: {e}"
            logger.error("Failed to predict PP3 and BP4 criteria. Error: {}", e)
            self.prediction = None

        return self.prediction, self.comment
