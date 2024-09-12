"""pytest configuration file"""

from typing import Iterator

import pytest
from fastapi.testclient import TestClient

from src.core.config import Config
from src.main import app


@pytest.fixture
def api_base_url() -> str:
    return "https://reev.cubi.bihealth.org/internal/proxy"


@pytest.fixture
def config(api_base_url: str) -> Config:
    return Config(api_base_url=api_base_url)


@pytest.fixture(scope="module")
def client() -> Iterator[TestClient]:
    """Fixture with a test client for the FastAPI app."""
    with TestClient(app) as c:
        yield c
