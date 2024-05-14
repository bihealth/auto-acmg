"""pytest configuration file"""

import pytest

from src.core.config import HelperConfig


@pytest.fixture
def api_base_url() -> str:
    return "https://reev.cubi.bihealth.org/internal/proxy"


@pytest.fixture
def helper_config(api_base_url: str) -> HelperConfig:
    return HelperConfig(api_base_url=api_base_url)
