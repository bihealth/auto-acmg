import pytest
from pytest_httpx import HTTPXMock

from src.api.reev.dotty import DottyClient
from src.defs.dotty import DottySpdiResponse
from src.defs.genome_builds import GenomeRelease
from tests.utils import get_json_object

# -------- to_spdi ---------


@pytest.mark.asyncio
async def test_to_spdi_success(httpx_mock: HTTPXMock):
    """Test to_spdi method with a successful response."""
    mock_response = get_json_object("dotty/dotty_spdi_success.json")
    httpx_mock.add_response(
        method="GET",
        url="https://example.com/dotty/api/v1/to-spdi?q=test_query&assembly=GRCh38",
        json=mock_response,
        status_code=200,
    )

    client = DottyClient(api_base_url="https://example.com/dotty")
    response = client.to_spdi("test_query", GenomeRelease.GRCh38)
    assert response == DottySpdiResponse.model_validate(mock_response)


@pytest.mark.asyncio
async def test_to_spdi_failure(httpx_mock: HTTPXMock):
    """Test to_spdi method with a failed response."""
    mock_response = get_json_object("dotty/dotty_spdi_failure.json")
    httpx_mock.add_response(
        method="GET",
        url="https://example.com/dotty/api/v1/to-spdi?q=test_query&assembly=GRCh38",
        json=mock_response,
        status_code=200,
    )

    client = DottyClient(api_base_url="https://example.com/dotty")
    response = client.to_spdi("test_query", GenomeRelease.GRCh38)
    assert response == DottySpdiResponse.model_validate(mock_response)


@pytest.mark.asyncio
async def test_to_spdi_500(httpx_mock: HTTPXMock):
    """Test to_spdi method with a 500 response."""
    httpx_mock.add_response(
        method="GET",
        url="https://example.com/dotty/api/v1/to-spdi?q=test_query&assembly=GRCh38",
        status_code=500,
    )

    client = DottyClient(api_base_url="https://example.com/dotty")
    response = client.to_spdi("test_query", GenomeRelease.GRCh38)
    assert response is None
