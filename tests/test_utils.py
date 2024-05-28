from unittest.mock import MagicMock, patch

import pytest

from src.api.mehari import MehariClient
from src.defs.auto_pvs1 import GenomicStrand, SeqVarConsequence
from src.defs.genome_builds import GenomeRelease
from src.defs.mehari import GeneTranscripts, TranscriptsSeqVar
from src.defs.seqvar import SeqVar
from src.pvs1.seqvar_pvs1 import SeqVarTranscriptsHelper
from tests.utils import get_json_object


@pytest.fixture
def seqvar():
    return SeqVar(GenomeRelease.GRCh38, "1", 1000, "A", "T", "1:1000A>T")


@pytest.fixture
def ts_helper(seqvar):
    return SeqVarTranscriptsHelper(seqvar)


@pytest.fixture
def seqvar_transcripts(file_name: str = "mehari/mehari_seqvar_success.json"):
    return TranscriptsSeqVar.model_validate(get_json_object(file_name)).result


@pytest.fixture
def gene_transcripts(file_name: str = "mehari/mehari_genes_success.json"):
    return GeneTranscripts.model_validate(get_json_object(file_name)).transcripts


#: Mock the Exon class
class MockExon:
    def __init__(self, altStartI, altEndI, altCdsStartI=None, altCdsEndI=None, cigar="", ord=None):
        self.altStartI = altStartI
        self.altEndI = altEndI
        self.altCdsStartI = altCdsStartI if altCdsStartI is not None else altStartI
        self.altCdsEndI = altCdsEndI if altCdsEndI is not None else altEndI
        self.cigar = cigar
        self.ord = ord


#: Mock the CdsInfo class
class MockCdsInfo:
    def __init__(
        self, start_codon, stop_codon, cds_start, cds_end, exons, cds_strand=GenomicStrand.Plus
    ):
        self.start_codon = start_codon
        self.stop_codon = stop_codon
        self.cds_start = cds_start
        self.cds_end = cds_end
        self.cds_strand = cds_strand
        self.exons = exons


# === SeqVarTranscriptsHelper ===


def test_get_ts_info_success(ts_helper):
    """Test get_ts_info method with a successful response."""
    # Mock the actual data that would be returned from the Mehari API
    ts_helper.seqvar_ts_info = TranscriptsSeqVar.model_validate(
        get_json_object("mehari/mehari_seqvar_success.json")
    )
    ts_helper.seqvar_transcript = TranscriptsSeqVar.model_validate(
        get_json_object("mehari/mehari_seqvar_success.json")
    ).result
    ts_helper.gene_ts_info = GeneTranscripts.model_validate(
        get_json_object("mehari/mehari_genes_success.json")
    )
    ts_helper.gene_transcript = GeneTranscripts.model_validate(
        get_json_object("mehari/mehari_genes_success.json")
    ).transcripts
    ts_helper.consequence = SeqVarConsequence.InitiationCodon

    seqvar_transcript, gene_transcript, seqvar_ts_info, gene_ts_info, consequence = (
        ts_helper.get_ts_info()
    )

    assert seqvar_transcript is not None
    assert gene_transcript is not None
    assert seqvar_ts_info is not None
    assert gene_ts_info is not None
    assert consequence == SeqVarConsequence.InitiationCodon


def test_get_ts_info_failure(ts_helper):
    """Test get_ts_info method with a failed response."""
    seqvar_transcript, gene_transcript, seqvar_ts_info, gene_ts_info, consequence = (
        ts_helper.get_ts_info()
    )

    assert seqvar_transcript is None
    assert gene_transcript is None
    assert seqvar_ts_info == []
    assert gene_ts_info == []
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

    assert ts_helper.seqvar_ts_info == []
    assert ts_helper.gene_ts_info == []
    assert ts_helper.seqvar_ts_info == []
    assert ts_helper.gene_ts_info == []
    assert ts_helper.HGNC_id is ""
    assert len(ts_helper.HGVSs) == 0


@patch.object(MehariClient, "get_seqvar_transcripts")
@patch.object(MehariClient, "get_gene_transcripts")
def test_initialize_no_seqvar(mock_get_gene_transcripts, mock_get_seqvar_transcripts, ts_helper):
    # Mock the MehariClient methods for CI
    mock_get_gene_transcripts, mock_get_seqvar_transcripts = MagicMock(), MagicMock()
    ts_helper.initialize()
    assert ts_helper.gene_ts_info == []
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
            "mehari/larp7_mehari_gene.json",
            "mehari/larp7_mehari_seqvar.json",
            "NM_001267039.2",
        ),
        (
            ["NM_001267039.2", "NM_001370974.1"],
            "mehari/larp7_mehari_gene.json",
            "mehari/larp7_mehari_seqvar.json",
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
        (["invalid"], "mehari/larp7_mehari_gene.json", "mehari/larp7_mehari_seqvar.json"),
        ([], "mehari/larp7_mehari_gene.json", "mehari/larp7_mehari_seqvar.json"),
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
