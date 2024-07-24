"""pytest configuration file"""

import os

import pytest

from src.core.config import Config


@pytest.fixture
def api_base_url() -> str:
    return "https://reev.cubi.bihealth.org/internal/proxy"


@pytest.fixture
def config(api_base_url: str) -> Config:
    return Config(api_base_url=api_base_url)


@pytest.fixture(autouse=True)
def set_env_vars():
    """Set and unset environment variables for tests."""
    # Disable cache for tests
    os.environ["USE_CACHE"] = "0"
    yield
    # Unset environment variables after tests
    os.environ.pop("USE_CACHE", None)
