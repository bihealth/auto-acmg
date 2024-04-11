"""Dotty API client."""

import requests
from pydantic import ValidationError

from src.core.config import settings
from src.defs.dotty import DottySpdiResponse
from src.defs.genome_builds import GenomeRelease

#: Dotty API base URL
DOTTI_API_BASE_URL = f"{settings.API_REEV_URL}/dotty"


class DottyClient:
    def __init__(self, api_base_url: str = DOTTI_API_BASE_URL):
        self.api_base_url = api_base_url

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
        response = requests.get(url)
        try:
            response.raise_for_status()
            return DottySpdiResponse.model_validate(response.json())
        except requests.RequestException:
            return None
        except ValidationError as e:
            return None
