import responses

from src.api.annonars import AnnonarsClient
from src.defs.annonars_range import AnnonarsRangeResponse
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from tests.utils import get_json_object

#: Example sequence variant
example_seqvar = SeqVar(
    chrom="1",
    pos=1000,
    delete="A",
    insert="T",
    genome_release=GenomeRelease.GRCh38,
    user_repr="1:1000A>T",
)


@responses.activate
def test_get_variant_from_range_success():
    """Test to_annonar method with a successful response."""
    mock_response = get_json_object("annonars/annonars_range_success.json")
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
    assert response == AnnonarsRangeResponse.model_validate(mock_response)


@responses.activate
def test_get_variant_from_range_failure():
    """Test to_annonar method with a failed response."""
    mock_response = get_json_object("annonars/annonars_range_failure.json")
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
    assert response == AnnonarsRangeResponse.model_validate(mock_response)


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


@responses.activate
def test_get_variant_info_success():
    """Test get_variant_info method with a successful response."""
    mock_response = get_json_object("annonars/annonars_variant_success.json")
    responses.add(
        responses.GET,
        f"https://example.com/annonars/annos/variant?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&pos={example_seqvar.pos}&reference={example_seqvar.delete}&alternative={example_seqvar.insert}",
        json=mock_response,
        status=200,
    )

    client = AnnonarsClient(api_base_url="https://example.com/annonars")
    response = client.get_variant_info(example_seqvar)
    assert response == AnnonarsVariantResponse.model_validate(mock_response)


@responses.activate
def test_get_variant_info_failure():
    """Test get_variant_info method with a failed response."""
    mock_response = get_json_object("annonars/annonars_variant_failure.json")
    responses.add(
        responses.GET,
        f"https://example.com/annonars/annos/variant?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&pos={example_seqvar.pos}&reference={example_seqvar.delete}&alternative={example_seqvar.insert}",
        json=mock_response,
        status=200,
    )

    client = AnnonarsClient(api_base_url="https://example.com/annonars")
    response = client.get_variant_info(example_seqvar)
    assert response == AnnonarsVariantResponse.model_validate(mock_response)


@responses.activate
def test_get_variant_info_500():
    """Test get_variant_info method with a 500 response."""
    responses.add(
        responses.GET,
        f"https://example.com/annonars/annos/variant?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&pos={example_seqvar.pos}&reference={example_seqvar.delete}&alternative={example_seqvar.insert}",
        status=500,
    )

    client = AnnonarsClient(api_base_url="https://example.com/annonars")
    response = client.get_variant_info(example_seqvar)
    assert response == None
