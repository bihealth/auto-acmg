from unittest.mock import MagicMock, patch

import pytest

from src.api.mehari import MehariClient
from src.defs.autopvs1 import SeqVarConsequence
from src.defs.genome_builds import GenomeRelease
from src.defs.mehari import GeneTranscripts, TranscriptsSeqVar
from src.defs.seqvar import SeqVar
from src.seqvar_pvs1 import SeqVarPVS1, SeqVarPVS1Helpers, SeqVarTranscriptsHelper
from tests.utils import get_json_object


@pytest.fixture
def seqvar():
    return SeqVar(GenomeRelease.GRCh38, "1", 1000, "A", "T", "1:1000A>T")


@pytest.fixture
def ts_helper(seqvar):
    return SeqVarTranscriptsHelper(seqvar)


@pytest.fixture
def seqvar_transcripts():
    return TranscriptsSeqVar.model_validate(get_json_object("mehari_seqvar_success.json")).result


@pytest.fixture
def gene_transcripts():
    return GeneTranscripts.model_validate(get_json_object("mehari_genes_success.json")).transcripts


# === SeqVarPVS1Helpers ===


# === SeqVarTranscriptsHelper ===


def test_get_ts_info_success(ts_helper):
    """Test get_ts_info method with a successful response."""
    # Mock the actual data that would be returned from the Mehari API
    ts_helper.seqvar_transcript = TranscriptsSeqVar.model_validate(
        get_json_object("mehari_seqvar_success.json")
    ).result
    ts_helper.gene_transcript = GeneTranscripts.model_validate(
        get_json_object("mehari_genes_success.json")
    ).transcripts
    ts_helper.consequence = SeqVarConsequence.InitiationCodon

    seqvar_transcript, gene_transcript, consequence = ts_helper.get_ts_info()

    assert seqvar_transcript is not None
    assert gene_transcript is not None
    assert consequence == SeqVarConsequence.InitiationCodon


def test_get_ts_info_failure(ts_helper):
    """Test get_ts_info method with a failed response."""
    seqvar_transcript, gene_transcript, consequence = ts_helper.get_ts_info()

    assert seqvar_transcript is None
    assert gene_transcript is None
    assert consequence == SeqVarConsequence.NotSet


@patch.object(MehariClient, "get_seqvar_transcripts")
@patch.object(MehariClient, "get_gene_transcripts")
def test_initialize_success(
    mock_get_gene_transcripts,
    mock_get_seqvar_transcripts,
    seqvar,
    ts_helper,
    seqvar_transcripts,
    gene_transcripts,
):
    # Mock successful responses
    mock_get_seqvar_transcripts.return_value = MagicMock(result=seqvar_transcripts)
    mock_get_gene_transcripts.return_value = MagicMock(transcripts=gene_transcripts)

    ts_helper.seqvar = seqvar
    ts_helper.initialize()

    assert ts_helper.seqvar_ts_info is seqvar_transcripts
    assert ts_helper.gene_ts_info is gene_transcripts
    assert ts_helper.HGNC_id is seqvar_transcripts[0].gene_id
    assert ts_helper.HGVSs is not None


@patch.object(MehariClient, "get_seqvar_transcripts")
@patch.object(MehariClient, "get_gene_transcripts")
def test_initialize_failure(
    mock_get_gene_transcripts, mock_get_seqvar_transcripts, seqvar, ts_helper
):
    # Mock failed responses
    mock_get_seqvar_transcripts.return_value = None
    mock_get_gene_transcripts.return_value = None

    ts_helper.seqvar = seqvar
    ts_helper.initialize()

    assert ts_helper.seqvar_ts_info is None
    assert ts_helper.gene_ts_info is None
    assert ts_helper.HGNC_id is ""
    assert len(ts_helper.HGVSs) == 0


def test_initialize_no_seqvar(ts_helper):
    ts_helper.initialize()
    assert len(ts_helper.seqvar_ts_info) == 0
    assert ts_helper.gene_ts_info is None
    assert ts_helper.HGNC_id is ""
    assert len(ts_helper.HGVSs) == 0


@pytest.mark.parametrize(
    "consequence_input, expected_consequence",
    [
        (["splice_region_variant"], SeqVarConsequence.SpliceSites),
        (["splice_donor_variant"], SeqVarConsequence.SpliceSites),
        (["frameshift_variant"], SeqVarConsequence.NonsenseFrameshift),
        (["initiator_codon_variant"], SeqVarConsequence.InitiationCodon),
        (["unknown_consequence"], SeqVarConsequence.NotSet),
        (["regulatory_region_amplification"], SeqVarConsequence.NotSet),
        ([""], SeqVarConsequence.NotSet),
        ([], SeqVarConsequence.NotSet),
    ],
)
def test_get_consequence_various_cases(consequence_input, expected_consequence, seqvar_transcripts):
    """Test get_consequence method with various cases."""
    mock_transcript = seqvar_transcripts[0]
    mock_transcript.consequences = consequence_input
    consequence = SeqVarTranscriptsHelper._get_consequence(mock_transcript)
    assert consequence == expected_consequence


def test_get_consequence_none_input():
    """Test get_consequence method with None input."""
    consequence = SeqVarTranscriptsHelper._get_consequence(None)
    assert consequence == SeqVarConsequence.NotSet


# TODO: Add more use cases for the choose_transcript method
# E.g. - multiple transcripts, multiple genes
#      - one transcript, multiple genes (or vice versa)
#      - no transcripts, no genes
#      - no transcripts, one gene (or vice versa)
@pytest.mark.parametrize(
    "hgvss, gene_ts_file, seqvar_ts_file, expected_hgvs",
    [
        (
            ["NM_001267039.2"],
            "larp7_mehari_gene.json",
            "larp7_mehari_seqvar.json",
            "NM_001267039.2",
        ),
        (
            ["NM_001267039.2", "NM_001370974.1"],
            "larp7_mehari_gene.json",
            "larp7_mehari_seqvar.json",
            "NM_001267039.2",
        ),
    ],
)
def test_choose_transcript_success(hgvss, gene_ts_file, seqvar_ts_file, expected_hgvs, ts_helper):
    """Test choose_transcript method."""
    ts_helper.seqvar_ts_info = TranscriptsSeqVar.model_validate(
        get_json_object(seqvar_ts_file)
    ).result
    ts_helper.gene_ts_info = GeneTranscripts.model_validate(
        get_json_object(gene_ts_file)
    ).transcripts

    seqvar_ts, gene_ts = ts_helper._choose_transcript(
        hgvss, ts_helper.seqvar_ts_info, ts_helper.gene_ts_info
    )
    assert seqvar_ts.feature_id == expected_hgvs


@pytest.mark.parametrize(
    "hgvss, gene_ts_file, seqvar_ts_file",
    [
        (["invalid"], "larp7_mehari_gene.json", "larp7_mehari_seqvar.json"),
        ([], "larp7_mehari_gene.json", "larp7_mehari_seqvar.json"),
    ],
)
def test_choose_transcript_invalid(hgvss, gene_ts_file, seqvar_ts_file, ts_helper):
    """Test choose_transcript method."""
    ts_helper.seqvar_ts_info = TranscriptsSeqVar.model_validate(
        get_json_object(seqvar_ts_file)
    ).result
    ts_helper.gene_ts_info = GeneTranscripts.model_validate(
        get_json_object(gene_ts_file)
    ).transcripts

    seqvar_ts, gene_ts = ts_helper._choose_transcript(
        hgvss, ts_helper.seqvar_ts_info, ts_helper.gene_ts_info
    )
    assert seqvar_ts == None


# === SeqVarPVS1 ===
