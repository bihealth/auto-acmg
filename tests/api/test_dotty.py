import responses

from src.api.dotty import DottyClient
from src.defs.dotty import DottySpdiResponse
from src.defs.genome_builds import GenomeRelease
from tests.utils import get_json_object


@responses.activate
def test_to_spdi_success():
    """Test to_spdi method with a successful response."""
    mock_response = get_json_object("dotty/dotty_spdi_success.json")
    responses.add(
        responses.GET,
        "https://example.com/dotty/api/v1/to-spdi?q=test_query&assembly=GRCh38",
        json=mock_response,
        status=200,
    )

    client = DottyClient(api_base_url="https://example.com/dotty")
    response = client.to_spdi("test_query", GenomeRelease.GRCh38)
    assert response == DottySpdiResponse.model_validate(mock_response)


@responses.activate
def test_to_spdi_failure():
    """Test to_spdi method with a failed response."""
    mock_response = get_json_object("dotty/dotty_spdi_failure.json")
    responses.add(
        responses.GET,
        "https://example.com/dotty/api/v1/to-spdi?q=test_query&assembly=GRCh38",
        json=mock_response,
        status=200,
    )

    client = DottyClient(api_base_url="https://example.com/dotty")
    response = client.to_spdi("test_query", GenomeRelease.GRCh38)
    assert response == DottySpdiResponse.model_validate(mock_response)


@responses.activate
def test_to_spdi_500():
    """Test to_spdi method with a 500 response."""
    responses.add(
        responses.GET,
        "https://example.com/dotty/api/v1/to-spdi?q=test_query&assembly=GRCh38",
        status=500,
    )

    client = DottyClient(api_base_url="https://example.com/dotty")
    response = client.to_spdi("test_query", GenomeRelease.GRCh38)
    assert response == None
