import pytest
from aioresponses import aioresponses

from src.api.dotty import DottyClient
from src.genome_builds import GenomeRelease


@pytest.mark.asyncio
async def test_to_spdi_success():
    """Test to_spdi method with a successful response."""
    with aioresponses() as m:
        mock_response = {"success": True, "value": "mocked response data"}
        m.get(
            "https://example.com/dotty/api/v1/to-spdi?q=test_query&assembly=GRCh38",
            payload=mock_response,
            status=200,
        )

        client = DottyClient(api_base_url="https://example.com/dotty")
        response = await client.to_spdi("test_query", GenomeRelease.GRCh38)

        assert response == mock_response


@pytest.mark.asyncio
async def test_to_spdi_failure():
    """Test to_spdi method with a failed response."""
    with aioresponses() as m:
        m.get("https://example.com/dotty/api/v1/to-spdi?q=test_query&assembly=GRCh38", status=404)

        client = DottyClient(api_base_url="https://example.com/dotty")
        response = await client.to_spdi("test_query", GenomeRelease.GRCh38)

        assert response is None


@pytest.mark.asyncio
async def test_to_spdi_grch37_assembly():
    """Test to_spdi method with default assembly."""
    with aioresponses() as m:
        mock_response = {"success": True, "value": "mocked response data"}
        m.get(
            "https://example.com/dotty/api/v1/to-spdi?q=test_query&assembly=GRCh37",
            payload=mock_response,
            status=200,
        )

        client = DottyClient(api_base_url="https://example.com/dotty")
        response = await client.to_spdi("test_query", GenomeRelease.GRCh37)

        assert response == mock_response
