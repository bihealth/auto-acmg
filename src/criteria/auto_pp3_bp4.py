"""Implementation of the PP3 and BP4 criteria."""

from typing import List, Optional, Union

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.annonars_variant import VariantResult
from src.defs.auto_acmg import PP3BP4, MissenseScores
from src.defs.exceptions import AutoAcmgBaseException, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


class AutoPP3BP4:
    """Class for automatic PP3 and BP4 prediction."""

    def __init__(
        self,
        seqvar: SeqVar,
        genome_release: GenomeRelease,
        variant_info: VariantResult,
        *,
        config: Config,
    ):
        """Initialize the class."""
        self.seqvar = seqvar
        self.genome_release = genome_release
        self.variant_info = variant_info
        self.config = config
        self.annonars_client = AnnonarsClient(api_base_url=config.api_base_url_annonars)
        self.prediction: PP3BP4 | None = None

    @staticmethod
    def _convert_score_value(score_value: Optional[Union[str, float, int]]) -> Optional[float]:
        """Convert score value to float."""
        if score_value is None:
            return None
        if isinstance(score_value, (float, int)):
            return score_value
        try:
            return max(float(score) for score in score_value.split(";") if score != ".")
        except ValueError as e:
            logger.error("Failed to convert score value to float. Error: {}", e)
            raise

    def _is_pathogenic_score(self, variant_info: VariantResult) -> bool:
        """Check if any of the pathogenic scores meet the threshold."""
        if not variant_info.dbnsfp:
            logger.error("Missing dbNSFP data.")
            raise MissingDataError("Missing dbNSFP data.")
        for score in MissenseScores:
            score_value = self._convert_score_value(getattr(variant_info.dbnsfp, score.name, None))
            if score_value is not None and score.pathogenic_threshold is not None:
                if score_value >= score.pathogenic_threshold:
                    logger.debug(
                        "Pathogenic score({}): {} >= {}",
                        score.name,
                        score_value,
                        score.pathogenic_threshold,
                    )
                    return True
        return False

    def _is_benign_score(self, variant_info: VariantResult) -> bool:
        """Check if any of the benign scores meet the threshold."""
        if not variant_info.dbnsfp:
            logger.error("Missing dbNSFP data.")
            raise MissingDataError("Missing dbNSFP data.")
        for score in MissenseScores:
            score_value = self._convert_score_value(getattr(variant_info.dbnsfp, score.name, None))
            if score_value is not None and score.benign_threshold is not None:
                if score_value <= score.benign_threshold:
                    logger.debug(
                        "Benign score({}): {} <= {}",
                        score.name,
                        score_value,
                        score.benign_threshold,
                    )
                    return True
        return False

    @staticmethod
    def _is_pathogenic_spliceai(variant_info: VariantResult) -> bool:
        """Check if any of the pathogenic scores meet the threshold."""
        if not variant_info.cadd:
            logger.error("Missing CADD data.")
            raise MissingDataError("Missing CADD data.")
        spliceai_scores = [
            variant_info.cadd.SpliceAI_acc_gain,
            variant_info.cadd.SpliceAI_acc_loss,
            variant_info.cadd.SpliceAI_don_gain,
            variant_info.cadd.SpliceAI_don_loss,
        ]
        for score in spliceai_scores:
            if score is not None and score >= 0.5:
                logger.debug("Pathogenic SpliceAI score: {}", score)
                return True
        return False

    @staticmethod
    def _is_benign_spliceai(variant_info: VariantResult) -> bool:
        """Check if any of the pathogenic scores meet the threshold."""
        if not variant_info.cadd:
            logger.error("Missing CADD data.")
            raise MissingDataError("Missing CADD data.")
        spliceai_scores = [
            variant_info.cadd.SpliceAI_acc_gain,
            variant_info.cadd.SpliceAI_acc_loss,
            variant_info.cadd.SpliceAI_don_gain,
            variant_info.cadd.SpliceAI_don_loss,
        ]
        for score in spliceai_scores:
            if score is not None and score <= 0.1:
                logger.debug("Benign SpliceAI score: {}", score)
                return True
        return False

    def predict(self) -> Optional[PP3BP4]:
        """Predict PP3 and BP4 criteria."""
        self.prediction = PP3BP4()

        try:
            if self.seqvar.chrom == "MT":
                self.prediction.PP3 = False
                self.prediction.BP4 = False
                return self.prediction

            if not self.variant_info:
                logger.error("Missing variant data.")
                raise MissingDataError("Missing variant data.")

            # Evaluate PP3 and BP4 criteria
            is_pathogenic = self._is_pathogenic_score(
                self.variant_info
            )  # or self._is_pathogenic_spliceai(self.variant_info)
            is_benign = self._is_benign_score(
                self.variant_info
            )  # or self._is_benign_spliceai(self.variant_info)
            self.prediction.PP3 = is_pathogenic
            self.prediction.BP4 = is_benign

        except AutoAcmgBaseException as e:
            logger.error("Failed to predict PP3 and BP4 criteria. Error: {}", e)
            self.prediction = None

        return self.prediction
