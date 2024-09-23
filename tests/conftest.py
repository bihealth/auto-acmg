"""pytest configuration file"""

from typing import Iterator

import pytest
from fastapi.testclient import TestClient

from src.core.config import settings
from src.main import app

# Override settings for testing
settings.API_REEV_URL = "https://reev.cubi.bihealth.org/internal/proxy"
settings.AUTO_ACMG_API_ANNONARS_URL = settings.API_REEV_URL + "/annonars"
settings.AUTO_ACMG_API_MEHARI_URL = settings.API_REEV_URL + "/mehari"
settings.AUTO_ACMG_API_DOTTY_URL = settings.API_REEV_URL + "/dotty"


@pytest.fixture(scope="module")
def client() -> Iterator[TestClient]:
    """Fixture with a test client for the FastAPI app."""
    with TestClient(app) as c:
        yield c
