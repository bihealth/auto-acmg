"""Annonars API client."""

from typing import Optional, Union

import httpx
from loguru import logger
from pydantic import ValidationError

from src.core.cache import Cache
from src.core.config import settings
from src.defs.annonars_gene import AnnonarsGeneResponse
from src.defs.annonars_range import AnnonarsCustomRangeResult, AnnonarsRangeResponse
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.exceptions import AnnonarsException
from src.defs.seqvar import SeqVar
from src.defs.strucvar import StrucVar

#: Annonars API base URL
ANNONARS_API_BASE_URL = settings.AUTO_ACMG_API_ANNONARS_URL or f"{settings.API_REEV_URL}/annonars"


class AnnonarsClient:
    def __init__(self, *, api_base_url: Optional[str] = None):
        #: Annonars API base URL
        self.api_base_url = api_base_url or ANNONARS_API_BASE_URL
        #: HTTPX client
        self.client = httpx.Client()
        #: Persistent cache for API responses
        self.cache = Cache()

    def _get_variant_from_range(
        self, variant: Union[SeqVar, StrucVar], start: int, stop: int
    ) -> AnnonarsRangeResponse:
        """Pull all variants within a range.

        Args:
            variant (Union[SeqVar, StrucVar]): Sequence or structural variant.
            start (int): Start position.
            stop (int): Stop position.

        Returns:
            AnnonarsRangeResponse: Annonars response.
        """
        if abs(stop - start) > 5000:
            raise AnnonarsException("Range is too large for a single request.")

        gr = variant.genome_release.name.lower()
        chromosome = variant.chrom
        url = (
            f"{self.api_base_url}/annos/range?"
            f"genome_release={gr}"
            f"&chromosome={chromosome}"
            f"&start={start}"
            f"&stop={stop}"
        )
        logger.debug("GET request to: {}", url)

        cached_response = self.cache.get(url)
        if cached_response:
            try:
                return AnnonarsRangeResponse.model_validate(cached_response)
            except ValidationError as e:
                logger.exception("Validation failed for cached data: {}", e)
                raise AnnonarsException("Cached data is invalid") from e

        response = self.client.get(url)
        if response.status_code != 200:
            logger.error("Request failed: {}", response.text)
            raise AnnonarsException(
                f"Request failed. Status code: {response.status_code}, Text: {response.text}"
            )
        try:
            response.raise_for_status()
            response_data = response.json()
            self.cache.add(url, response_data)
            return AnnonarsRangeResponse.model_validate(response_data)
        except ValidationError as e:
            logger.exception("Validation failed: {}", e)
            raise AnnonarsException("Annonars returned non-validating data.") from e

    def get_variant_from_range(
        self, variant: Union[SeqVar, StrucVar], start: int, stop: int
    ) -> AnnonarsCustomRangeResult:
        """A wrapper for _get_variant_from_range that handles ranges larger than 5000.

        Args:
            variant (Union[SeqVar, StrucVar]): Sequence or structural variant.
            start (int): Start position.
            stop (int): Stop position.

        Returns:
            AnnonarsRangeResponse: Annonars response.
        """
        MAX_RANGE_SIZE = 5000
        res = AnnonarsCustomRangeResult()
        current_start = start
        while current_start < stop:
            current_stop = min(current_start + MAX_RANGE_SIZE - 1, stop)
            response = self._get_variant_from_range(variant, current_start, current_stop)
            if res.gnomad_genomes is None:
                res.gnomad_genomes = response.result.gnomad_genomes
            else:
                res.gnomad_genomes.extend(response.result.gnomad_genomes or [])
            if res.clinvar is None:
                res.clinvar = response.result.clinvar
            else:
                res.clinvar.extend(response.result.clinvar or [])
            current_start = current_stop + 1
        return res

    def get_variant_info(self, seqvar: SeqVar) -> AnnonarsVariantResponse:
        """Get variant information from Annonars.

        Args:
            seqvar (SeqVar): Sequence variant.

        Returns:
            Any: Annonars response.
        """
        url = (
            f"{self.api_base_url}/annos/variant?"
            f"genome_release={seqvar.genome_release.name.lower()}"
            f"&chromosome={seqvar.chrom}"
            f"&pos={seqvar.pos}"
            f"&reference={seqvar.delete}"
            f"&alternative={seqvar.insert}"
        )
        logger.debug("GET request to: {}", url)

        cached_response = self.cache.get(url)
        if cached_response:
            try:
                return AnnonarsVariantResponse.model_validate(cached_response)
            except ValidationError as e:
                logger.exception("Validation failed for cached data: {}", e)
                raise AnnonarsException("Cached data is invalid") from e

        response = self.client.get(url)
        if response.status_code != 200:
            logger.error("Request failed: {}", response.text)
            raise AnnonarsException(
                f"Request failed. Status code: {response.status_code}, Text: {response.text}"
            )
        try:
            response.raise_for_status()
            response_data = response.json()
            self.cache.add(url, response_data)
            return AnnonarsVariantResponse.model_validate(response_data)
        except ValidationError as e:
            logger.exception("Validation failed: {}", e)
            raise AnnonarsException("Annonars returned non-validating data.") from e

    def get_gene_info(self, hgnc_id: str) -> AnnonarsGeneResponse:
        """Get gene information from Annonars.

        Args:
            seqvar (SeqVar): Sequence variant.

        Returns:
            Any: Annonars response.
        """
        url = f"{self.api_base_url}/genes/info?hgnc_id={hgnc_id}"
        logger.debug("GET request to: {}", url)

        cached_response = self.cache.get(url)
        if cached_response:
            try:
                return AnnonarsGeneResponse.model_validate(cached_response)
            except ValidationError as e:
                logger.exception("Validation failed for cached data: {}", e)
                raise AnnonarsException("Cached data is invalid") from e

        response = self.client.get(url)
        if response.status_code != 200:
            logger.error("Request failed: {}", response.text)
            raise AnnonarsException(
                f"Request failed. Status code: {response.status_code}, Text: {response.text}"
            )
        try:
            response_data = response.json()
            self.cache.add(url, response_data)
            return AnnonarsGeneResponse.model_validate(response_data)
        except ValidationError as e:
            logger.exception("Validation failed: {}", e)
            raise AnnonarsException("Annonars returned non-validating data.") from e
