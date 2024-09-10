import pytest
from pytest_httpx import HTTPXMock

from src.api.mehari import MehariClient
from src.defs.exceptions import MehariException
from src.defs.genome_builds import GenomeRelease
from src.defs.mehari import GeneTranscripts, TranscriptsSeqVar, TranscriptsStrucVar
from src.defs.seqvar import SeqVar
from src.defs.strucvar import StrucVar, StrucVarType
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

#: Example structural variant
example_strucvar = StrucVar(
    sv_type=StrucVarType.DEL,
    chrom="1",
    start=1000,
    stop=2000,
    genome_release=GenomeRelease.GRCh38,
    user_repr="1:1000-2000DEL",
)

#: Example HGNC gene ID
example_hgnc_id = "HGNC:1234"


# ---------- get_seqvar_transcripts ----------------


@pytest.mark.asyncio
async def test_get_seqvar_transcripts_success(httpx_mock: HTTPXMock):
    """Test get_transcripts method with a successful response."""
    mock_response = get_json_object("mehari/DCDC2_seqvar.json")
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/mehari/seqvars/csq?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&position={example_seqvar.pos}&reference={example_seqvar.delete}&alternative={example_seqvar.insert}",
        json=mock_response,
        status_code=200,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    response = client.get_seqvar_transcripts(example_seqvar)
    assert response == TranscriptsSeqVar.model_validate(mock_response)


@pytest.mark.asyncio
async def test_get_seqvar_transcripts_failure(httpx_mock: HTTPXMock):
    """Test get_transcripts method with a failed response."""
    mock_response = get_json_object("mehari/seqvar_failure.json")
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/mehari/seqvars/csq?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&position={example_seqvar.pos}&reference={example_seqvar.delete}&alternative={example_seqvar.insert}",
        json=mock_response,
        status_code=200,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    response = client.get_seqvar_transcripts(example_seqvar)
    assert response == TranscriptsSeqVar.model_validate(mock_response)


@pytest.mark.asyncio
async def test_get_seqvar_transcripts_500(httpx_mock: HTTPXMock):
    """Test get_transcripts method with a 500 response."""
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/mehari/seqvars/csq?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&position={example_seqvar.pos}&reference={example_seqvar.delete}&alternative={example_seqvar.insert}",
        status_code=500,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    with pytest.raises(MehariException):
        client.get_seqvar_transcripts(example_seqvar)


# ---------- get_strucvar_transcripts ----------------


@pytest.mark.asyncio
async def test_get_strucvar_transcripts_success(httpx_mock: HTTPXMock):
    """Test get_strucvar_transcripts method with a successful response."""
    mock_response = get_json_object("mehari/KRTAP9-2_strucvar.json")
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/mehari/strucvars/csq?genome_release={example_strucvar.genome_release.name.lower()}&chromosome={example_strucvar.chrom}&start={example_strucvar.start}&stop={example_strucvar.stop}&sv_type={example_strucvar.sv_type.name.upper()}",
        json=mock_response,
        status_code=200,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    response = client.get_strucvar_transcripts(example_strucvar)
    assert response == TranscriptsStrucVar.model_validate(mock_response)


@pytest.mark.asyncio
async def test_get_strucvar_transcripts_failure(httpx_mock: HTTPXMock):
    """Test get_strucvar_transcripts method with a failed response."""
    mock_response = get_json_object("mehari/strucvar_failure.json")
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/mehari/strucvars/csq?genome_release={example_strucvar.genome_release.name.lower()}&chromosome={example_strucvar.chrom}&start={example_strucvar.start}&stop={example_strucvar.stop}&sv_type={example_strucvar.sv_type.name.upper()}",
        json=mock_response,
        status_code=200,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    response = client.get_strucvar_transcripts(example_strucvar)
    assert response == TranscriptsStrucVar.model_validate(mock_response)


@pytest.mark.asyncio
async def test_get_strucvar_transcripts_500(httpx_mock: HTTPXMock):
    """Test get_strucvar_transcripts method with a 500 response."""
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/mehari/strucvars/csq?genome_release={example_strucvar.genome_release.name.lower()}&chromosome={example_strucvar.chrom}&start={example_strucvar.start}&stop={example_strucvar.stop}&sv_type={example_strucvar.sv_type.name.upper()}",
        status_code=500,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    with pytest.raises(MehariException):
        client.get_strucvar_transcripts(example_strucvar)


# ---------- get_gene_transcripts ----------------


@pytest.mark.asyncio
async def test_get_gene_transcripts_success(httpx_mock: HTTPXMock):
    """Test get_gene_transcripts method with a successful response."""
    mock_response = get_json_object("mehari/HAL_gene.json")
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/mehari/genes/txs?hgncId={example_hgnc_id}&genomeBuild=GENOME_BUILD_GRCH38",
        json=mock_response,
        status_code=200,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    response = client.get_gene_transcripts(example_hgnc_id, GenomeRelease.GRCh38)
    assert response == GeneTranscripts.model_validate(mock_response)


@pytest.mark.asyncio
async def test_get_gene_transcripts_failure(httpx_mock: HTTPXMock):
    """Test get_gene_transcripts method with a failed response."""
    mock_response = get_json_object("mehari/gene_failure.json")
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/mehari/genes/txs?hgncId={example_hgnc_id}&genomeBuild=GENOME_BUILD_GRCH38",
        json=mock_response,
        status_code=200,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    with pytest.raises(MehariException):
        client.get_gene_transcripts(example_hgnc_id, GenomeRelease.GRCh38)


@pytest.mark.asyncio
async def test_get_gene_transcripts_500(httpx_mock: HTTPXMock):
    """Test get_gene_transcripts method with a 500 response."""
    httpx_mock.add_response(
        method="GET",
        url=f"https://example.com/mehari/genes/txs?hgncId={example_hgnc_id}&genomeBuild=GENOME_BUILD_GRCH38",
        status_code=500,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    with pytest.raises(MehariException):
        client.get_gene_transcripts(example_hgnc_id, GenomeRelease.GRCh38)
