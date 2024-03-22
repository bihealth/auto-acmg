import responses

from src.api.dotty import DottyClient
from src.genome_builds import GenomeRelease


@responses.activate
def test_to_spdi_success():
    """Test to_spdi method with a successful response."""
    mock_response = {"success": True, "value": "mocked response data"}
    responses.add(
        responses.GET,
        "https://example.com/dotty/api/v1/to-spdi?q=test_query&assembly=GRCh38",
        json=mock_response,
        status=200,
    )

    client = DottyClient(api_base_url="https://example.com/dotty")
    response = client.to_spdi("test_query", GenomeRelease.GRCh38)
    # assert response == mock_response


@responses.activate
def test_to_spdi_failure():
    """Test to_spdi method with a failed response."""
    responses.add(
        responses.GET,
        "https://example.com/dotty/api/v1/to-spdi?q=test_query&assembly=GRCh38",
        status=404,
    )

    client = DottyClient(api_base_url="https://example.com/dotty")
    response = client.to_spdi("test_query", GenomeRelease.GRCh38)

    # assert response is None


@responses.activate
def test_to_spdi_grch37_assembly():
    """Test to_spdi method with GRCh37 assembly."""
    mock_response = {"success": True, "value": "mocked response data"}
    responses.add(
        responses.GET,
        "https://example.com/dotty/api/v1/to-spdi?q=test_query&assembly=GRCh37",
        json=mock_response,
        status=200,
    )

    client = DottyClient(api_base_url="https://example.com/dotty")
    response = client.to_spdi("test_query", GenomeRelease.GRCh37)

    # assert response == mock_response
