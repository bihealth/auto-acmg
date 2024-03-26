import responses

from src.api.mehari import MehariClient
from src.genome_builds import GenomeRelease
from src.models.mehari import GeneTranscripts, TranscriptsSeqVar
from src.seqvar import SeqVar
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
def test_get_seqvar_transcripts_success():
    """Test get_transcripts method with a successful response."""
    mock_response = get_json_object("mehari_seqvar_success.json")
    responses.add(
        responses.GET,
        f"https://example.com/mehari/seqvars/csq?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&position={example_seqvar.pos}&reference={example_seqvar.insert}&alternative={example_seqvar.delete}",
        json=mock_response,
        status=200,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    response = client.get_seqvar_transcripts(example_seqvar)
    assert response == TranscriptsSeqVar(**mock_response)


@responses.activate
def test_get_seqvar_transcripts_failure():
    """Test get_transcripts method with a failed response."""
    mock_response = get_json_object("mehari_seqvar_failure.json")
    responses.add(
        responses.GET,
        f"https://example.com/mehari/seqvars/csq?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&position={example_seqvar.pos}&reference={example_seqvar.insert}&alternative={example_seqvar.delete}",
        json=mock_response,
        status=200,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    response = client.get_seqvar_transcripts(example_seqvar)
    assert response == TranscriptsSeqVar(**mock_response)


@responses.activate
def test_get_seqvar_transcripts_500():
    """Test get_transcripts method with a 500 response."""
    responses.add(
        responses.GET,
        f"https://example.com/mehari/seqvars/csq?genome_release={example_seqvar.genome_release.name.lower()}&chromosome={example_seqvar.chrom}&position={example_seqvar.pos}&reference={example_seqvar.insert}&alternative={example_seqvar.delete}",
        status=500,
    )

    client = MehariClient(api_base_url="https://example.com/mehari")
    response = client.get_seqvar_transcripts(example_seqvar)
    assert response == None
