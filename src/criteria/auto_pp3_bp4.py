"""Implementation of the PP3 and BP4 criteria."""

from typing import List, Optional

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.annonars_variant import VariantResult
from src.defs.auto_acmg import PP3BP4
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
    def _get_best_pathogenic_score(variant_info: VariantResult) -> Optional[float]:
        """Get the best pathogenic score available."""
        # Avoid linting errors
        assert variant_info.cadd is not None and variant_info.dbnsfp is not None
        scores = [
            # variant_info.mutpred2,
            variant_info.cadd.PolyPhenVal,
            # variant_info.revel,
            # variant_info.bayesdel,
            # variant_info.vest4,
            # variant_info.phylop,
        ]
        scores = [score for score in scores if score is not None]
        return max(scores) if scores else None  # type: ignore

    @staticmethod
    def _get_best_benign_score(variant_info: VariantResult) -> Optional[float]:
        """Get the best benign score available."""
        # Avoid linting errors
        assert variant_info.cadd is not None and variant_info.dbnsfp is not None
        scores = [
            # variant_info.dbnsfp.REVEL_rankscore,
            # variant_info.dbnsfp.MutPred_rankscore,
            variant_info.cadd.PolyPhenVal,
            # variant_info.bayesdel,
            # variant_info.vest4,
            # variant_info.phylop,
        ]
        scores = [score for score in scores if score is not None]
        return min(scores) if scores else None  # type: ignore

    def _predict_spliceai(self, variant_info: VariantResult) -> Optional[float]:
        """Predict splice site alterations using SpliceAI."""
        # TODO: Implement this method.
        return None

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

            best_pathogenic_score = self._get_best_pathogenic_score(self.variant_info)
            best_benign_score = self._get_best_benign_score(self.variant_info)
            spliceai_score = self._predict_spliceai(self.variant_info)

            # Evaluate PP3
            if best_pathogenic_score and best_pathogenic_score > 0.8:
                self.prediction.PP3 = True
            elif spliceai_score and spliceai_score > 0.8:
                self.prediction.PP3 = True
            else:
                self.prediction.PP3 = False

            # Evaluate BP4
            if best_benign_score and best_benign_score < 0.2:
                self.prediction.BP4 = True
            elif spliceai_score and spliceai_score < 0.2:
                self.prediction.BP4 = True
            else:
                self.prediction.BP4 = False

        except AutoAcmgBaseException as e:
            logger.error("Failed to predict PP3 and BP4 criteria. Error: {}", e)
            self.prediction = None

        return self.prediction
