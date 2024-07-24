import hashlib
import json
import os
from typing import Optional

from loguru import logger

from src.core.config import settings


class Cache:
    """Cache class to store the results of the Annonars API calls."""

    #: Cache directory
    CACHE_DIR = os.path.join(settings.PATH_TO_ROOT, "src", "cache")

    def _get_cache_filename(self, url: str) -> str:
        """Generate a cache filename based on the MD5 hash of the URL."""
        url_hash = hashlib.md5(url.encode()).hexdigest()
        return os.path.join(self.CACHE_DIR, f"{url_hash}.json")

    def get(self, url: str) -> Optional[dict]:
        """Check if a cached response exists and return it."""
        cache_filename = self._get_cache_filename(url)
        if os.path.exists(cache_filename):
            logger.debug("Loading cached response from: {}", cache_filename)
            with open(cache_filename, "r") as cache_file:
                return json.load(cache_file)
        return None

    def add(self, url: str, response_data: dict) -> None:
        """Cache the response data."""
        cache_filename = self._get_cache_filename(url)
        logger.debug("Caching response to: {}", cache_filename)
        with open(cache_filename, "w") as cache_file:
            json.dump(response_data, cache_file)
