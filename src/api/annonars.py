"""Annonars API client."""

import aiohttp

from src.core.config import settings
from src.seqvar import SeqVar

#: Annonars API base URL
ANNONARS_API_BASE_URL = f"{settings.API_REEV_URL}/annonars"


class AnnonarsClient:
    def __init__(self, api_base_url: str = ANNONARS_API_BASE_URL):
        self.api_base_url = api_base_url

    async def get_variant_from_range(self, seqvar: SeqVar, start: int, stop: int) -> dict | None:
        """
        Pull all variants within a range.

        :param seqvar: Sequence variant
        :type seqvar: SeqVar
        :param start: Start position
        :type start: int
        :param stop: End position
        :type stop: int
        :return: Variants
        :rtype: dict | None
        """
        url = (
            f"{self.api_base_url}/annos/range?"
            f"genome_release={seqvar.genome_release.name.lower()}"
            f"&chromosome={seqvar.chrom}"
            f"&start={start}"
            f"&stop={stop}"
        )
        async with aiohttp.ClientSession() as session:
            async with session.get(url) as response:
                if response.status == 200:
                    return await response.json()
                else:
                    return None
