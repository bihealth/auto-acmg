from unittest.mock import MagicMock, patch

import pytest

from src.api.mehari import MehariClient
from src.defs.autopvs1 import (
    PVS1Prediction,
    PVS1PredictionSeqVarPath,
    SeqVarConsequence,
)
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


#: Mock the Exon class
class MockExon:
    def __init__(self, altStartI, altEndI, altCdsStartI=None, altCdsEndI=None, cigar="", ord=None):
        self.altStartI = altStartI
        self.altEndI = altEndI
        self.altCdsStartI = altCdsStartI if altCdsStartI is not None else altStartI
        self.altCdsEndI = altCdsEndI if altCdsEndI is not None else altEndI
        self.cigar = cigar
        self.ord = ord


# === SeqVarPVS1Helpers ===


# TODO: Check if the termination number is correct
@pytest.mark.parametrize(
    "pHGVS,expected_termination",
    [
        ("NM_031475.2:p.Gln98*", 98),
        ("NM_031475.2:p.Ala586Glyfs*73", 586 + 73),
        ("NP_000305.3:p.Arg378SerfsTer5", 378 + 5),
        ("p.Arg97Glyfs*26", 97 + 26),
        ("p.Arg97GlyfsX26", 97 + 26),
        ("p.Arg97GlyfsTer26", 97 + 26),
        # ("p.Arg97fs", -1),  # No termination number provided
        ("p.Gly100Ter", 100),
        ("p.Cys24*", 24),
        ("p.Ala2X", 2),
        ("p.Tyr10Ter", 10),
        ("p.Gln98*", 98),
        ("p.Ala586Gly", -1),  # No frameshift or termination codon
    ],
)
def test_get_pHGVS_termination(pHGVS, expected_termination):
    """Test the _get_pHGVS_termination method."""
    termination = SeqVarPVS1Helpers._get_pHGVS_termination(pHGVS)
    assert termination == expected_termination, f"Failed for pHGVS: {pHGVS}"


# TODO: Check if the exon number is correct
@pytest.mark.parametrize(
    "exons, pHGVS, hgnc_id, expected_result",
    [
        # Specific gene with HGNC ID "HGNC:4284" should always return True
        ([MockExon(0, 100, 0, 100)], "p.Gln98*", "HGNC:4284", True),
        # Variant in the last exon should not undergo NMD
        ([MockExon(0, 200, 0, 150), MockExon(200, 400, 150, 350)], "p.Gln300*", "HGNC:1234", False),
        # Variant in the penultimate exon, more than 50 nt from the end, should undergo NMD
        # ([MockExon(0, 100, 0, 100), MockExon(100, 300, 100, 300)], "p.Gln50*", "HGNC:1234", True),
        # Single exon variants should not undergo NMD
        ([MockExon(0, 100, 0, 100)], "p.Gln50*", "HGNC:1234", False),
        # Variant in the penultimate exon, within the last 50 nt, should not undergo NMD
        ([MockExon(0, 100, 0, 100), MockExon(100, 300, 100, 300)], "p.Gln95*", "HGNC:1234", False),
    ],
)
def test_undergo_nmd(exons, pHGVS, hgnc_id, expected_result):
    """Test the _undergo_nmd method."""
    result = SeqVarPVS1Helpers()._undergo_nmd(exons, pHGVS, hgnc_id)
    assert result == expected_result, f"Failed for hgnc_id: {hgnc_id}, pHGVS: {pHGVS}"


@pytest.mark.parametrize(
    "transcript_tags,expected_result",
    [
        ([], False),
        (["NonRelevantTag"], False),
        (["ManeSelect"], True),
        (["SomeOtherTag", "ManeSelect"], True),
        (["maneselect"], False),  # Case-sensitive check
        (["MANESELECT"], False),  # Case-sensitive check
        (["SomeTag", "AnotherTag"], False),
    ],
)
def test_in_biologically_relevant_transcript(transcript_tags, expected_result):
    """Test the _in_biologically_relevant_transcript method."""
    result = SeqVarPVS1Helpers._in_biologically_relevant_transcript(transcript_tags)
    assert result == expected_result, f"Failed for transcript_tags: {transcript_tags}"


def test_critical4protein_function():
    """Test the _critical4protein_function method."""
    pass


def test_lof_is_frequent_in_population():
    """Test the _lof_is_frequent_in_population method."""
    pass


# TODO: Check if the exon number is correct
@pytest.mark.parametrize(
    "exons, pHGVS, expected_result",
    [
        # ([MockExon(0, 300)], "p.Gln100*", True),  # Simple case where LoF removes more than 10% of a 100-codon protein
        (
            [MockExon(0, 300)],
            "p.Gln90*",
            False,
        ),  # LoF variant at codon 90, removes exactly 10% of a 100-codon protein
        # ([MockExon(0, 300), MockExon(300, 900)], "p.Gln200fs*1", True),  # Frameshift early in the protein
        ([MockExon(0, 900)], "p.Gln850*", False),  # Truncation removes less than 10% of the protein
        ([MockExon(0, 100), MockExon(100, 500)], "p.Arg50X", True),  # Early nonsense mutation
        ([MockExon(0, 500)], "p.Arg490Ter", False),  # Truncation very close to the end
    ],
)
def test_lof_removes_more_then_10_percent_of_protein(exons, pHGVS, expected_result):
    """Test the _lof_removes_more_then_10_percent_of_protein method."""
    result = SeqVarPVS1Helpers._lof_removes_more_then_10_percent_of_protein(pHGVS, exons)
    assert (
        result == expected_result
    ), f"Expected {expected_result} for pHGVS: {pHGVS} with exon lengths: {[exon.altEndI - exon.altStartI for exon in exons]}"


def test_exon_skipping_or_cryptic_ss_disruption():
    """Test the _exon_skipping_or_cryptic_ss_disruption method."""
    pass


def test_alternative_start_codon():
    """Test the _alternative_start_codon method."""
    pass


def test_upstream_pathogenic_variant():
    """Test the _upstream_pathogenic_variant method."""
    pass


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


def test_init(seqvar):
    """Test the initialization of SeqVarPVS1."""
    pvs1 = SeqVarPVS1(seqvar)
    assert pvs1.seqvar == seqvar
    assert pvs1._seqvar_transcript is None
    assert pvs1._gene_transcript is None
    assert pvs1._consequence == SeqVarConsequence.NotSet
    assert pvs1.HGVS == ""
    assert pvs1.pHGVS == ""
    assert pvs1.tHGVS == ""
    assert pvs1.HGNC_id == ""
    assert len(pvs1.transcript_tags) == 0
    assert len(pvs1.exons) == 0
    assert pvs1.cds_pos is None
    assert pvs1.prediction == PVS1Prediction.NotPVS1
    assert pvs1.prediction_path == PVS1PredictionSeqVarPath.NotSet


@patch.object(SeqVarTranscriptsHelper, "get_ts_info", autospec=True)
@patch.object(SeqVarTranscriptsHelper, "initialize", autospec=True)
def test_get_pvs1_prediction_success(
    mock_initialize, mock_get_ts_info, seqvar, seqvar_transcripts, gene_transcripts
):
    mock_get_ts_info.return_value = (
        seqvar_transcripts[0],
        gene_transcripts[0],
        SeqVarConsequence.InitiationCodon,
    )

    pvs1 = SeqVarPVS1(seqvar)
    pvs1.initialize()

    assert pvs1._seqvar_transcript is seqvar_transcripts[0]
    assert pvs1._gene_transcript is gene_transcripts[0]
    assert pvs1._consequence is SeqVarConsequence.InitiationCodon


@patch.object(SeqVarTranscriptsHelper, "get_ts_info", autospec=True)
@patch.object(SeqVarTranscriptsHelper, "initialize", autospec=True)
def test_get_pvs1_prediction_failure(mock_initialize, mock_get_ts_info, seqvar):
    mock_get_ts_info.return_value = (None, None, SeqVarConsequence.NotSet)

    pvs1 = SeqVarPVS1(seqvar)
    pvs1.initialize()

    assert pvs1._seqvar_transcript is None
    assert pvs1._gene_transcript is None
    assert pvs1._consequence is SeqVarConsequence.NotSet


# TODO: Add integration tests for the verify_PVS1 method


def test_get_prediction_default(seqvar):
    pvs1 = SeqVarPVS1(seqvar)
    prediction, path = pvs1.get_prediction()
    assert prediction == PVS1Prediction.NotPVS1
    assert path == PVS1PredictionSeqVarPath.NotSet
