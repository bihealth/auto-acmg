"""Annonars API client."""

from typing import Any, Optional

import requests
from loguru import logger
from pydantic import ValidationError

from src.core.config import settings
from src.defs.annonars_range import AnnonarsRangeResponse
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.exceptions import AnnonarsException
from src.defs.seqvar import SeqVar

#: Annonars API base URL
ANNONARS_API_BASE_URL = f"{settings.API_REEV_URL}/annonars"


class AnnonarsClient:
    def __init__(self, api_base_url: str = ANNONARS_API_BASE_URL):
        self.api_base_url = api_base_url

    def get_variant_from_range(self, seqvar: SeqVar, start: int, stop: int) -> AnnonarsRangeResponse:
        """Pull all variants within a range.

        Args:
            seqvar (SeqVar): Sequence variant.
            start (int): Start position.
            stop (int): Stop position.

        Returns:
            AnnonarsRangeResponse: Annonars response.
        """
        url = (
            f"{self.api_base_url}/annos/range?"
            f"genome_release={seqvar.genome_release.name.lower()}"
            f"&chromosome={seqvar.chrom}"
            f"&start={start}"
            f"&stop={stop}"
        )
        logger.debug("GET request to: {}", url)
        response = requests.get(url)
        try:
            response.raise_for_status()
            return AnnonarsRangeResponse.model_validate(response.json())
        except requests.RequestException as e:
            logger.exception("Request failed: {}", e)
            raise AnnonarsException("Failed to get variant information.") from e
        except ValidationError as e:
            logger.exception("Validation failed: {}", e)
            raise AnnonarsException("Annonars returned non-validating data.") from e

    def get_variant_info(self, seqvar: SeqVar) -> Optional[AnnonarsVariantResponse]:
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
        response = requests.get(url)
        try:
            response.raise_for_status()
            return AnnonarsVariantResponse.model_validate(response.json())
        except requests.RequestException as e:
            logger.exception("Request failed: {}", e)
            raise AnnonarsException("Failed to get variant information.") from e
        except ValidationError as e:
            logger.exception("Validation failed: {}", e)
            raise AnnonarsException("Annonars returned non-validating data.") from e
