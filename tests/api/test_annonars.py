import responses

from src.api.annonars import AnnonarsClient
from src.genome_builds import GenomeRelease
from src.seqvar import SeqVar
from src.types.annonars import AnnonarsRangeResponse
from tests.utils import get_json_object

#: Example sequence variant
example_seqvar = SeqVar(
    chrom="1",
    pos=1000,
    delete="A",
    insert="T",
    genome_release=GenomeRelease.GRCh38,
    user_representation="1:1000A>T",
)


@responses.activate
def test_get_variant_from_range_success():
    """Test to_annonar method with a successful response."""
    mock_response = get_json_object("annonars_range_success.json")
    start = 1000
    stop = 2000
    responses.add(
        responses.GET,
        f"https://example.com/annonars/annos/range?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&start={start}&stop={stop}",
        json=mock_response,
        status=200,
    )

    client = AnnonarsClient(api_base_url="https://example.com/annonars")
    response = client.get_variant_from_range(example_seqvar, start, stop)
    assert response == AnnonarsRangeResponse(**mock_response)


@responses.activate
def test_get_variant_from_range_failure():
    """Test to_annonar method with a failed response."""
    mock_response = get_json_object("annonars_range_failure.json")
    start = 1000
    stop = 999
    responses.add(
        responses.GET,
        f"https://example.com/annonars/annos/range?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&start={start}&stop={stop}",
        json=mock_response,
        status=200,
    )

    client = AnnonarsClient(api_base_url="https://example.com/annonars")
    response = client.get_variant_from_range(example_seqvar, start, stop)
    assert response == AnnonarsRangeResponse(**mock_response)


@responses.activate
def test_get_variant_from_range_500():
    """Test to_annonar method with a 500 response."""
    start = 1000
    stop = 2000
    responses.add(
        responses.GET,
        f"https://example.com/annonars/annos/range?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&start={start}&stop={stop}",
        status=500,
    )

    client = AnnonarsClient(api_base_url="https://example.com/annonars")
    response = client.get_variant_from_range(example_seqvar, start, stop)
    assert response == None
