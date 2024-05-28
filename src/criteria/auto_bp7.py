"""Implementation of BP7 criteria."""

from typing import Optional

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.annonars_variant import VariantResult
from src.defs.auto_acmg import BP7
from src.defs.exceptions import AutoAcmgBaseException, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.utils import SeqVarTranscriptsHelper


class AutoBP7:
    """Class for automatic BP7 prediction."""

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
        self.prediction: BP7 | None = None

    def _get_variant_info(self, seqvar: SeqVar) -> Optional[VariantResult]:
        """Get variant information from Annonars.

        Returns:
            VariantResult: Annonars response.
        """
        try:
            logger.debug("Getting variant information for {}.", seqvar)
            return self.annonars_client.get_variant_info(seqvar).result
        except AutoAcmgBaseException as e:
            logger.error("Failed to get variant information. Error: {}", e)
            return None

    def _check_proximity_to_pathogenic_variants(self, seqvar: SeqVar) -> bool:
        """Check for pathogenic variants +/- 2bp of the position in ClinVar."""
        response = self.annonars_client.get_variant_from_range(
            seqvar, seqvar.pos - 2, seqvar.pos + 2
        )
        if response and response.result.clinvar:
            pathogenic_variants = [
                v
                for v in response.result.clinvar
                if v.referenceAssertions
                and v.referenceAssertions[0].clinicalSignificance
                in [
                    "CLINICAL_SIGNIFICANCE_LIKELY_PATHOGENIC",
                    "CLINICAL_SIGNIFICANCE_PATHOGENIC",
                ]
            ]
            return len(pathogenic_variants) > 0
        return False

    def _check_proximity_to_splice_site(self, seqvar: SeqVar) -> bool:
        """Check if the variant is closer than 2bp to a splice site."""
        # TODO: Implement this method.
        return False

    def _predict_spliceai(self, seqvar: SeqVar) -> bool:
        """Predict splice site alterations using SpliceAI."""
        # TODO: Implement this method.
        return False

    def predict(self) -> Optional[BP7]:
        """Predict BP7 criterion."""
        self.prediction = BP7()
        try:
            if self.seqvar.chrom == "MT":
                # skipped according to McCormick et al. (2020).
                self.prediction.BP7 = False
                return self.prediction

            if self._check_proximity_to_pathogenic_variants(self.seqvar):
                self.prediction.BP7 = False
                return self.prediction
            if self._check_proximity_to_splice_site(self.seqvar):
                self.prediction.BP7 = False
                return self.prediction
            if self._predict_spliceai(self.seqvar):
                self.prediction.BP7 = True
            else:
                self.prediction.BP7 = False
        except AutoAcmgBaseException as e:
            logger.error("Failed to predict BP7 criterion. Error: {}", e)
            self.prediction = None

        return self.prediction
