"""Implementation of the PP3 and BP4 criteria."""

from typing import List, Optional

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
    def _is_pathogenic_score(variant_info: VariantResult) -> bool:
        """Check if any of the pathogenic scores meet the threshold."""
        for score in MissenseScores:
            score_value = getattr(variant_info.cadd, score.name, None)
            if score_value is not None and score.pathogenic_threshold is not None:
                if score_value >= score.pathogenic_threshold:
                    return True
        return False

    @staticmethod
    def _is_benign_score(variant_info: VariantResult) -> bool:
        """Check if any of the benign scores meet the threshold."""
        for score in MissenseScores:
            score_value = getattr(variant_info.cadd, score.name, None)
            if score_value is not None and score.benign_threshold is not None:
                if score_value <= score.benign_threshold:
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

            if not self.variant_info or not self.variant_info.cadd or not self.variant_info.dbnsfp:
                logger.error("Missing CADD or DBNSFP data.")
                raise MissingDataError("Missing CADD or DBNSFP data.")

            # Evaluate PP3 and BP4 criteria
            is_pathogenic = self._is_pathogenic_score(self.variant_info)
            is_benign = self._is_benign_score(self.variant_info)
            self.prediction.PP3 = is_pathogenic
            self.prediction.BP4 = is_benign

        except AutoAcmgBaseException as e:
            logger.error("Failed to predict PP3 and BP4 criteria. Error: {}", e)
            self.prediction = None

        return self.prediction
