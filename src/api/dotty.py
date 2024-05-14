"""Dotty API client."""

from typing import Optional

import requests
from loguru import logger
from pydantic import ValidationError

from src.core.config import settings
from src.defs.dotty import DottySpdiResponse
from src.defs.genome_builds import GenomeRelease

#: Dotty API base URL
DOTTI_API_BASE_URL = f"{settings.API_REEV_URL}/dotty"


class DottyClient:
    def __init__(self, *, api_base_url: Optional[str] = None):
        self.api_base_url = api_base_url or DOTTI_API_BASE_URL

    def to_spdi(self, query: str, assembly: GenomeRelease = GenomeRelease.GRCh38) -> DottySpdiResponse | None:
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
        response = requests.get(url)
        try:
            response.raise_for_status()
            return DottySpdiResponse.model_validate(response.json())
        except requests.RequestException as e:
            logger.exception("Request failed: {}", e)
            return None
        except ValidationError as e:
            logger.exception("Validation failed: {}", e)
            return None
