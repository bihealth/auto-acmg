"""pytest configuration file"""

import pytest

from src.core.config import Config


@pytest.fixture
def api_base_url() -> str:
    return "https://reev.cubi.bihealth.org/internal/proxy"


@pytest.fixture
def config(api_base_url: str) -> Config:
    return Config(api_base_url=api_base_url)
