import pytest
import responses

from src.api.mehari import MehariClient
from src.defs.exceptions import MehariException
from src.defs.genome_builds import GenomeRelease
from src.defs.mehari import GeneTranscripts, TranscriptsSeqVar
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

#: Example HGNC gene ID
example_hgnc_id = "HGNC:1234"


@responses.activate
def test_get_seqvar_transcripts_success():
    """Test get_transcripts method with a successful response."""
    mock_response = get_json_object("mehari/mehari_seqvar_success.json")
    responses.add(
        responses.GET,
        f"https://example.com/mehari/seqvars/csq?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&position={example_seqvar.pos}&reference={example_seqvar.delete}&alternative={example_seqvar.insert}",
        json=mock_response,
        status=200,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    response = client.get_seqvar_transcripts(example_seqvar)
    assert response == TranscriptsSeqVar.model_validate(mock_response)


@responses.activate
def test_get_seqvar_transcripts_failure():
    """Test get_transcripts method with a failed response."""
    mock_response = get_json_object("mehari/mehari_seqvar_failure.json")
    responses.add(
        responses.GET,
        f"https://example.com/mehari/seqvars/csq?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&position={example_seqvar.pos}&reference={example_seqvar.delete}&alternative={example_seqvar.insert}",
        json=mock_response,
        status=200,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    response = client.get_seqvar_transcripts(example_seqvar)
    assert response == TranscriptsSeqVar.model_validate(mock_response)


@responses.activate
def test_get_seqvar_transcripts_500():
    """Test get_transcripts method with a 500 response."""
    responses.add(
        responses.GET,
        f"https://example.com/mehari/seqvars/csq?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&position={example_seqvar.pos}&reference={example_seqvar.delete}&alternative={example_seqvar.insert}",
        status=500,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    with pytest.raises(MehariException):
        client.get_seqvar_transcripts(example_seqvar)


@responses.activate
def test_get_gene_transcripts_success():
    """Test get_gene_transcripts method with a successful response."""
    mock_response = get_json_object("mehari/mehari_genes_success.json")
    responses.add(
        responses.GET,
        f"https://example.com/mehari/genes/txs?hgncId={example_hgnc_id}&genomeBuild=GENOME_BUILD_GRCH38",
        json=mock_response,
        status=200,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    response = client.get_gene_transcripts(example_hgnc_id, GenomeRelease.GRCh38)
    assert response == GeneTranscripts.model_validate(mock_response)


@responses.activate
def test_get_gene_transcripts_failure():
    """Test get_gene_transcripts method with a failed response."""
    mock_response = get_json_object("mehari/mehari_genes_failure.json")
    responses.add(
        responses.GET,
        f"https://example.com/mehari/genes/txs?hgncId={example_hgnc_id}&genomeBuild=GENOME_BUILD_GRCH38",
        json=mock_response,
        status=200,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    with pytest.raises(MehariException):
        client.get_gene_transcripts(example_hgnc_id, GenomeRelease.GRCh38)


@responses.activate
def test_get_gene_transcripts_500():
    """Test get_gene_transcripts method with a 500 response."""
    responses.add(
        responses.GET,
        f"https://example.com/mehari/genes/txs?hgncId={example_hgnc_id}&genomeBuild=GENOME_BUILD_GRCH38",
        status=500,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    with pytest.raises(MehariException):
        client.get_gene_transcripts(example_hgnc_id, GenomeRelease.GRCh38)
