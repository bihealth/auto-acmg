"""Dotty API client."""

import aiohttp

from src.core.config import settings
from src.genome_builds import GenomeRelease

#: Dotty API base URL
DOTTI_API_BASE_URL = f"{settings.API_REEV_URL}/dotty"


class DottyClient:
    def __init__(self, api_base_url: str = DOTTI_API_BASE_URL):
        self.api_base_url = api_base_url

    async def to_spdi(
        self, query: str, assembly: GenomeRelease = GenomeRelease.GRCh38
    ) -> dict | None:
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
        async with aiohttp.ClientSession() as session:
            async with session.get(url) as response:
                if response.status == 200:
                    return await response.json()
                else:
                    return None
