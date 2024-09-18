"""Mehari API client."""

from typing import Optional

import httpx
from loguru import logger
from pydantic import ValidationError

from src.core.cache import Cache
from src.core.config import settings
from src.defs.exceptions import MehariException
from src.defs.genome_builds import GenomeRelease
from src.defs.mehari import GeneTranscripts, TranscriptsSeqVar, TranscriptsStrucVar
from src.defs.seqvar import SeqVar
from src.defs.strucvar import StrucVar

#: Mehari API base URL
MEHARI_API_BASE_URL = settings.AUTO_ACMG_API_MEHARI_URL or f"{settings.API_REEV_URL}/mehari"


class MehariClient:
    def __init__(self, *, api_base_url: Optional[str] = None):
        #: Mehari API base URL
        self.api_base_url = api_base_url or MEHARI_API_BASE_URL
        #: HTTPX client
        self.client = httpx.Client()
        #: Persistent cache for API responses
        self.cache = Cache()

    def get_seqvar_transcripts(self, seqvar: SeqVar) -> TranscriptsSeqVar:
        """
        Get transcripts for a sequence variant.

        :param seqvar: Sequence variant
        :type seqvar: SeqVar
        :return: Transcripts
        :rtype: TranscriptsSeqVar
        :raises MehariException: if the request failed
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

        response = self.client.get(url)
        if response.status_code != 200:
            logger.error("Request failed: {}", response.text)
            raise MehariException(
                f"Request failed. Status code: {response.status_code}, Text: {response.text}"
            )
        try:
            response.raise_for_status()
            response_data = response.json()
            self.cache.add(url, response_data)
            return TranscriptsSeqVar.model_validate(response_data)
        except ValidationError as e:
            logger.exception("Validation failed: {}", e)
            raise MehariException("Mehari API returned invalid data") from e

    def get_strucvar_transcripts(self, strucvar: StrucVar) -> TranscriptsStrucVar:
        """
        Get transcripts for a structural variant.

        :param strucvar: Structural variant
        :type strucvar: StrucVar
        :return: Transcripts
        :rtype: TranscriptsStrucVar
        :raises MehariException: if the request failed
        """
        url = (
            f"{self.api_base_url}/strucvars/csq?"
            f"genome_release={strucvar.genome_release.name.lower()}"
            f"&chromosome={strucvar.chrom}"
            f"&start={strucvar.start}"
            f"&stop={strucvar.stop}"
            f"&sv_type={strucvar.sv_type.name.upper()}"
        )
        logger.debug("GET request to: {}", url)

        cached_response = self.cache.get(url)
        if cached_response:
            try:
                return TranscriptsStrucVar.model_validate(cached_response)
            except ValidationError as e:
                logger.exception("Validation failed for cached data: {}", e)
                raise MehariException("Cached data is invalid") from e

        response = self.client.get(url)
        if response.status_code != 200:
            logger.error("Request failed: {}", response.text)
            raise MehariException(
                f"Request failed. Status code: {response.status_code}, Text: {response.text}"
            )
        try:
            response.raise_for_status()
            response_data = response.json()
            self.cache.add(url, response_data)
            return TranscriptsStrucVar.model_validate(response_data)
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
        :rtype: GeneTranscripts
        :raises MehariException: if the request failed
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

        response = httpx.get(url)
        if response.status_code != 200:
            logger.error("Request failed: {}", response.text)
            raise MehariException(
                f"Request failed. Status code: {response.status_code}, Text: {response.text}"
            )
        try:
            response.raise_for_status()
            response_data = response.json()
            self.cache.add(url, response_data)
            return GeneTranscripts.model_validate(response_data)
        except ValidationError as e:
            logger.exception("Validation failed: {}", e)
            raise MehariException("Mehari API returned invalid data") from e
