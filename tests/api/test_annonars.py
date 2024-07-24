import pytest
from pytest_httpx import HTTPXMock

from src.api.annonars import AnnonarsClient
from src.defs.annonars_gene import AnnonarsGeneResponse
from src.defs.annonars_range import AnnonarsRangeResponse
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.exceptions import AnnonarsException
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


@pytest.mark.asyncio
async def test_get_variant_from_range_success(httpx_mock: HTTPXMock):
    """Test to_annonar method with a successful response."""
    mock_response = get_json_object("annonars/PCSK9_range.json")
    start = 1000
    stop = 2000
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/annonars/annos/range?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&start={start}&stop={stop}",
        json=mock_response,
        status_code=200,
    )

    client = AnnonarsClient(api_base_url="https://example.com/annonars")
    response = client.get_variant_from_range(example_seqvar, start, stop)
    assert response == AnnonarsRangeResponse.model_validate(mock_response)


@pytest.mark.asyncio
async def test_get_variant_from_range_failure(httpx_mock: HTTPXMock):
    """Test to_annonar method with a failed response."""
    mock_response = get_json_object("annonars/range_failure.json")
    start = 1000
    stop = 999
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/annonars/annos/range?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&start={start}&stop={stop}",
        json=mock_response,
        status_code=200,
    )

    client = AnnonarsClient(api_base_url="https://example.com/annonars")
    response = client.get_variant_from_range(example_seqvar, start, stop)
    assert response == AnnonarsRangeResponse.model_validate(mock_response)


@pytest.mark.asyncio
async def test_get_variant_from_range_500(httpx_mock: HTTPXMock):
    """Test to_annonar method with a 500 response."""
    start = 1000
    stop = 2000
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/annonars/annos/range?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&start={start}&stop={stop}",
        status_code=500,
    )

    client = AnnonarsClient(api_base_url="https://example.com/annonars")
    with pytest.raises(AnnonarsException):
        client.get_variant_from_range(example_seqvar, start, stop)


@pytest.mark.asyncio
async def test_get_variant_info_success(httpx_mock: HTTPXMock):
    """Test get_variant_info method with a successful response."""
    mock_response = get_json_object("annonars/RP11_variant.json")
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/annonars/annos/variant?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&pos={example_seqvar.pos}&reference={example_seqvar.delete}&alternative={example_seqvar.insert}",
        json=mock_response,
        status_code=200,
    )

    client = AnnonarsClient(api_base_url="https://example.com/annonars")
    response = client.get_variant_info(example_seqvar)
    assert response == AnnonarsVariantResponse.model_validate(mock_response)


@pytest.mark.asyncio
async def test_get_variant_info_failure(httpx_mock: HTTPXMock):
    """Test get_variant_info method with a failed response."""
    mock_response = get_json_object("annonars/variant_failure.json")
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/annonars/annos/variant?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&pos={example_seqvar.pos}&reference={example_seqvar.delete}&alternative={example_seqvar.insert}",
        json=mock_response,
        status_code=200,
    )

    client = AnnonarsClient(api_base_url="https://example.com/annonars")
    response = client.get_variant_info(example_seqvar)
    assert response == AnnonarsVariantResponse.model_validate(mock_response)


@pytest.mark.asyncio
async def test_get_variant_info_500(httpx_mock: HTTPXMock):
    """Test get_variant_info method with a 500 response."""
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/annonars/annos/variant?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&pos={example_seqvar.pos}&reference={example_seqvar.delete}&alternative={example_seqvar.insert}",
        status_code=500,
    )

    client = AnnonarsClient(api_base_url="https://example.com/annonars")
    with pytest.raises(AnnonarsException):
        client.get_variant_info(example_seqvar)


@pytest.mark.asyncio
async def test_get_gene_info_success(httpx_mock: HTTPXMock):
    """Test get_gene_info method with a successful response."""
    mock_response = get_json_object("annonars/BRCA1_gene.json")
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/annonars/genes/info?hgnc_id=HGNC:1100",
        json=mock_response,
        status_code=200,
    )

    client = AnnonarsClient(api_base_url="https://example.com/annonars")
    response = client.get_gene_info("HGNC:1100")
    assert response == AnnonarsGeneResponse.model_validate(mock_response)


@pytest.mark.asyncio
async def test_get_gene_info_failure(httpx_mock: HTTPXMock):
    """Test get_gene_info method with a failed response."""
    mock_response = get_json_object("annonars/gene_failure.json")
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/annonars/genes/info?hgnc_id=HGNC:1100",
        json=mock_response,
        status_code=200,
    )

    client = AnnonarsClient(api_base_url="https://example.com/annonars")
    # with pytest.raises(AnnonarsException):
    response = client.get_gene_info("HGNC:1100")
    assert response == AnnonarsGeneResponse.model_validate(mock_response)


@pytest.mark.asyncio
async def test_get_gene_info_500(httpx_mock: HTTPXMock):
    """Test get_gene_info method with a 500 response."""
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/annonars/genes/info?hgnc_id=HGNC:1100",
        status_code=500,
    )

    client = AnnonarsClient(api_base_url="https://example.com/annonars")
    with pytest.raises(AnnonarsException):
        client.get_gene_info("HGNC:1100")
