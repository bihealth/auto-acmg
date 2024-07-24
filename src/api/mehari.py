"""Mehari API client."""

from typing import Optional

import requests
from loguru import logger
from pydantic import ValidationError

from src.core.cache import Cache
from src.core.config import settings
from src.defs.exceptions import MehariException
from src.defs.genome_builds import GenomeRelease
from src.defs.mehari import GeneTranscripts, TranscriptsSeqVar
from src.defs.seqvar import SeqVar

#: Mehari API base URL
MEHARI_API_BASE_URL = f"{settings.API_REEV_URL}/mehari"


class MehariClient:
    def __init__(self, *, api_base_url: Optional[str] = None):
        #: Mehari API base URL
        self.api_base_url = api_base_url or MEHARI_API_BASE_URL
        #: Persistent cache for API responses
        self.cache = Cache()

    def get_seqvar_transcripts(self, seqvar: SeqVar) -> TranscriptsSeqVar:
        """
        Get transcripts for a sequence variant.

        :param seqvar: Sequence variant
        :type seqvar: SeqVar
        :return: Transcripts
        :rtype: TranscriptsSeqVar | None
        """
        url = (
            f"{self.api_base_url}/seqvars/csq?"
            f"genome_release={seqvar.genome_release.name.lower()}"
            f"&chromosome={seqvar.chrom}"
            f"&position={seqvar.pos}"
            f"&reference={seqvar.delete}"
            f"&alternative={seqvar.insert}"
        )
        logger.debug("GET request to: {}", url)

        cached_response = self.cache.get(url)
        if cached_response:
            try:
                return TranscriptsSeqVar.model_validate(cached_response)
            except ValidationError as e:
                logger.exception("Validation failed for cached data: {}", e)
                raise MehariException("Cached data is invalid") from e

        response = requests.get(url)
        try:
            response.raise_for_status()
            response_data = response.json()
            self.cache.add(url, response_data)
            return TranscriptsSeqVar.model_validate(response_data)
        except requests.RequestException as e:
            logger.exception("Request failed: {}", e)
            raise MehariException("Request failed") from e
        except ValidationError as e:
            logger.exception("Validation failed: {}", e)
            raise MehariException("Mehari API returned invalid data") from e

    def get_gene_transcripts(self, hgnc_id: str, genome_build: GenomeRelease) -> GeneTranscripts:
        """ "
        Get transcripts for a gene.

        :param hgnc_id: HGNC gene ID
        :type hgnc_id: str
        :param genome_build: Genome build
        :type genome_build: GenomeRelease
        :return: Transcripts
        :rtype: GeneTranscripts | None
        """
        genome_build_mapping = {
            GenomeRelease.GRCh37: "GENOME_BUILD_GRCH37",
            GenomeRelease.GRCh38: "GENOME_BUILD_GRCH38",
        }
        url = (
            f"{self.api_base_url}/genes/txs?"
            f"hgncId={hgnc_id}"
            f"&genomeBuild={genome_build_mapping[genome_build]}"
        )
        logger.debug("GET request to: {}", url)

        cached_response = self.cache.get(url)
        if cached_response:
            try:
                return GeneTranscripts.model_validate(cached_response)
            except ValidationError as e:
                logger.exception("Validation failed for cached data: {}", e)
                raise MehariException("Cached data is invalid") from e

        response = requests.get(url)
        try:
            response.raise_for_status()
            response_data = response.json()
            self.cache.add(url, response_data)
            return GeneTranscripts.model_validate(response_data)
        except requests.RequestException as e:
            logger.exception("Request failed: {}", e)
            raise MehariException("Request failed") from e
        except ValidationError as e:
            logger.exception("Validation failed: {}", e)
            raise MehariException("Mehari API returned invalid data") from e
