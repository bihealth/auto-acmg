"""Dotty API client."""

from typing import Optional

import httpx
from loguru import logger
from pydantic import ValidationError

from src.core.cache import Cache
from src.core.config import settings
from src.defs.dotty import DottySpdiResponse
from src.defs.genome_builds import GenomeRelease

#: Dotty API base URL
DOTTI_API_BASE_URL = settings.AUTO_ACMG_API_DOTTY_URL or f"{settings.API_REEV_URL}/dotty"


class DottyClient:
    def __init__(self, *, api_base_url: Optional[str] = None):
        #: Dotty API base URL
        self.api_base_url = api_base_url or DOTTI_API_BASE_URL
        #: HTTPX client
        self.client = httpx.Client()
        #: Persistent cache for API responses
        self.cache = Cache()

    def to_spdi(
        self, query: str, assembly: GenomeRelease = GenomeRelease.GRCh38
    ) -> DottySpdiResponse | None:
        """
        Converts a variant to SPDI format.

        :param query: Variant query
        :type query: str
        :param assembly: Genome assembly
        :type assembly: GRChAssemblyType
        :return: SPDI format
        :rtype: dict | None
        """
        url = f"{self.api_base_url}/api/v1/to-spdi?q={query}&assembly={assembly.name}"
        logger.debug("GET request to: {}", url)

        cached_response = self.cache.get(url)
        if cached_response:
            try:
                return DottySpdiResponse.model_validate(cached_response)
            except ValidationError as e:
                logger.exception("Validation failed for cached data: {}", e)
                return None

        response = self.client.get(url)
        if response.status_code != 200:
            logger.error("Request failed: {}", response.text)
            return None
        try:
            response.raise_for_status()
            response_data = response.json()
            self.cache.add(url, response_data)
            return DottySpdiResponse.model_validate(response_data)
        except httpx.RequestError as e:
            logger.exception("Request failed: {}", e)
            return None
        except ValidationError as e:
            logger.exception("Validation failed: {}", e)
            return None
