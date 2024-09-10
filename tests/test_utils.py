from typing import List
from unittest.mock import MagicMock, patch

import pytest

from src.api.mehari import MehariClient
from src.defs.auto_acmg import GenomicStrand, SpliceType
from src.defs.auto_pvs1 import SeqVarPVS1Consequence
from src.defs.exceptions import AlgorithmError
from src.defs.genome_builds import GenomeRelease
from src.defs.mehari import GeneTranscripts, TranscriptsSeqVar, TranscriptsStrucVar
from src.defs.seqvar import SeqVar
from src.defs.strucvar import StrucVar, StrucVarType
from src.utils import SeqVarTranscriptsHelper, SplicingPrediction, StrucVarTranscriptsHelper
from tests.utils import get_json_object


@pytest.fixture
def seqvar():
    return SeqVar(GenomeRelease.GRCh38, "1", 1000, "A", "T", "1:1000A>T")


@pytest.fixture
def ts_helper(seqvar):
    return SeqVarTranscriptsHelper(seqvar)


@pytest.fixture
def seqvar_transcripts(file_name: str = "mehari/DCDC2_seqvar.json"):
    return TranscriptsSeqVar.model_validate(get_json_object(file_name)).result


@pytest.fixture
def gene_transcripts(file_name: str = "mehari/HAL_gene.json"):
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
        self,
        start_codon,
        stop_codon,
        cds_start,
        cds_end,
        exons,
        cds_strand=GenomicStrand.Plus,
    ):
        self.start_codon = start_codon
        self.stop_codon = stop_codon
        self.cds_start = cds_start
        self.cds_end = cds_end
        self.cds_strand = cds_strand
        self.exons = exons


# === SplicingPrediction ===


@pytest.fixture
def seqvar_ss():
    return SeqVar(genome_release=GenomeRelease.GRCh38, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def exons_ss():
    return [MagicMock(altEndI=110, altStartI=90)]


@pytest.fixture
def config_ss():
    mock_config = MagicMock()
    mock_config.api_base_url_annonars = "http://mock_api_base_url"
    mock_config.seqrepo_data_dir = "/mock/seqrepo/data/dir"
    return mock_config


# Patching the external dependencies (SeqRepo, AnnonarsClient, load_matrix5, load_matrix3)
@patch("src.utils.SeqRepo")
@patch("src.utils.AnnonarsClient")
@patch("src.utils.load_matrix5", return_value="mock_matrix5")
@patch("src.utils.load_matrix3", return_value="mock_matrix3")
@patch.object(SplicingPrediction, "_initialize_maxentscore")
def test_splicing_prediction_init(
    mock_initialize_maxentscore,
    mock_load_matrix3,
    mock_load_matrix5,
    mock_annonars_client,
    mock_seqrepo,
    seqvar_ss,
    exons_ss,
    config_ss,
):
    """Test the initialization of SplicingPrediction class."""
    mock_annonars_client.return_value = MagicMock()
    mock_seqrepo.return_value = MagicMock()

    splicing_prediction = SplicingPrediction(
        seqvar=seqvar_ss,
        strand=GenomicStrand.Plus,
        consequences=["splice_donor_variant"],
        exons=exons_ss,
        config=config_ss,
    )

    # Asserts
    assert splicing_prediction.seqvar == seqvar_ss, "SeqVar should be initialized correctly."
    assert (
        splicing_prediction.strand == GenomicStrand.Plus
    ), "Strand should be initialized correctly."
    assert splicing_prediction.exons == exons_ss, "Exons should be initialized correctly."
    assert (
        splicing_prediction.splice_type == SpliceType.Donor
    ), "Splice type should be correctly determined."
    assert splicing_prediction.config == config_ss, "Config should be initialized correctly."
    assert (
        splicing_prediction.annonars_client == mock_annonars_client.return_value
    ), "AnnonarsClient should be initialized correctly."
    assert (
        splicing_prediction.sr == mock_seqrepo.return_value
    ), "SeqRepo should be initialized correctly."
    assert splicing_prediction.matrix5 == "mock_matrix5", "Matrix5 should be loaded correctly."
    assert splicing_prediction.matrix3 == "mock_matrix3", "Matrix3 should be loaded correctly."
    assert (
        splicing_prediction.maxentscore_ref == -1.00
    ), "Initial maxentscore_ref should be set correctly."
    assert (
        splicing_prediction.maxentscore_alt == -1.00
    ), "Initial maxentscore_alt should be set correctly."
    assert (
        splicing_prediction.maxent_foldchange == 1.00
    ), "Initial maxent_foldchange should be set correctly."
    (
        mock_initialize_maxentscore.assert_called_once(),
        "MaxEntScan initialization should be called.",
    )


@pytest.mark.parametrize(
    "strand,expected_splice_type",
    [
        (GenomicStrand.Plus, SpliceType.Donor),
        (GenomicStrand.Minus, SpliceType.Donor),
    ],
)
@patch("src.utils.SeqRepo")
@patch("src.utils.AnnonarsClient")
@patch("src.utils.load_matrix5", return_value="mock_matrix5")
@patch("src.utils.load_matrix3", return_value="mock_matrix3")
@patch.object(SplicingPrediction, "_initialize_maxentscore")
def test_splicing_prediction_splice_type(
    mock_initialize_maxentscore,
    mock_load_matrix3,
    mock_load_matrix5,
    mock_annonars_client,
    mock_seqrepo,
    strand,
    expected_splice_type,
    seqvar_ss,
    exons_ss,
    config_ss,
):
    """Test splice type determination during initialization."""
    splicing_prediction = SplicingPrediction(
        seqvar=seqvar_ss,
        strand=strand,
        consequences=["splice_donor_variant"],
        exons=exons_ss,
        config=config_ss,
    )

    assert (
        splicing_prediction.splice_type == expected_splice_type
    ), "Splice type should be correctly determined based on consequences."


# Patching the internal methods of SplicingPrediction that are called within _initialize_maxentscore
@patch("src.utils.SeqRepo")
@patch("src.utils.AnnonarsClient")
@patch.object(SplicingPrediction, "_find_reference_sequence", return_value=(90, 111, "ATGCGTACG"))
@patch.object(SplicingPrediction, "_generate_alt_sequence", return_value="ATGCGTACG")
@patch.object(SplicingPrediction, "_calculate_maxentscore", return_value=(5.0, 3.0, 0.6))
def test_initialize_maxentscore(
    mock_calculate_maxentscore,
    mock_generate_alt_sequence,
    mock_find_reference_sequence,
    mock_annonars_client,
    mock_seqrepo,
    seqvar_ss,
    exons_ss,
    config_ss,
):
    """Test the _initialize_maxentscore method in SplicingPrediction."""

    # Initialize the mocked SeqRepo and AnnonarsClient
    mock_annonars_client.return_value = MagicMock()
    mock_seqrepo.return_value = MagicMock()

    splicing_prediction = SplicingPrediction(
        seqvar=seqvar_ss,
        strand=GenomicStrand.Plus,
        consequences=["splice_donor_variant"],
        exons=exons_ss,
        config=config_ss,
    )

    # Call the _initialize_maxentscore method
    splicing_prediction._initialize_maxentscore()

    # Asserts
    assert mock_find_reference_sequence.call_count == 2
    assert mock_generate_alt_sequence.call_count == 2
    (
        mock_generate_alt_sequence.assert_called_with("ATGCGTACG", 90),
        "The method _generate_alt_sequence should be called with the correct parameters.",
    )
    assert mock_calculate_maxentscore.call_count == 2
    (
        mock_calculate_maxentscore.assert_called_with(
            "ATGCGTACG", "ATGCGTACG", splicing_prediction.splice_type
        ),
        "The method _calculate_maxentscore should be called with the correct parameters.",
    )

    assert (
        splicing_prediction.maxentscore_ref == 5.0
    ), "maxentscore_ref should be set to the correct value."
    assert (
        splicing_prediction.maxentscore_alt == 3.0
    ), "maxentscore_alt should be set to the correct value."
    assert (
        splicing_prediction.maxent_foldchange == 0.6
    ), "maxent_foldchange should be set to the correct value."


def test_determine_splice_type_acceptor():
    """Test determine_splice_type with splice acceptor variant."""
    consequences = ["missense_variant", "splice_acceptor_variant"]
    result = SplicingPrediction.determine_splice_type(consequences)
    assert result == SpliceType.Acceptor, "Should return SpliceType.Acceptor."


def test_determine_splice_type_donor():
    """Test determine_splice_type with splice donor variant."""
    consequences = ["missense_variant", "splice_donor_variant"]
    result = SplicingPrediction.determine_splice_type(consequences)
    assert result == SpliceType.Donor, "Should return SpliceType.Donor."


def test_determine_splice_type_acceptor_and_donor():
    """Test determine_splice_type with both splice acceptor and donor variants."""
    consequences = ["splice_acceptor_variant", "splice_donor_variant"]
    result = SplicingPrediction.determine_splice_type(consequences)
    assert result == SpliceType.Acceptor, "Should return SpliceType.Acceptor when both are present."


def test_determine_splice_type_unknown():
    """Test determine_splice_type with no splice-related variants."""
    consequences = ["missense_variant", "synonymous_variant"]
    result = SplicingPrediction.determine_splice_type(consequences)
    assert (
        result == SpliceType.Unknown
    ), "Should return SpliceType.Unknown when no splice variants are present."


def test_determine_splice_type_empty():
    """Test determine_splice_type with an empty consequence list."""
    consequences: List[str] = []
    result = SplicingPrediction.determine_splice_type(consequences)
    assert (
        result == SpliceType.Unknown
    ), "Should return SpliceType.Unknown for an empty consequence list."


def test_reverse_complement_uppercase():
    """Test reverse complement with uppercase DNA sequence."""
    seq = "ATCG"
    result = SplicingPrediction.reverse_complement(seq)
    assert result == "CGAT", "Should return the reverse complement 'CGAT' for input 'ATCG'."


def test_reverse_complement_lowercase():
    """Test reverse complement with lowercase DNA sequence."""
    seq = "atcg"
    result = SplicingPrediction.reverse_complement(seq)
    assert result == "cgat", "Should return the reverse complement 'cgat' for input 'atcg'."


def test_reverse_complement_mixed_case():
    """Test reverse complement with mixed case DNA sequence."""
    seq = "AtCg"
    result = SplicingPrediction.reverse_complement(seq)
    assert result == "cGaT", "Should return the reverse complement 'cGaT' for input 'AtCg'."


def test_reverse_complement_empty():
    """Test reverse complement with an empty string."""
    seq = ""
    result = SplicingPrediction.reverse_complement(seq)
    assert result == "", "Should return an empty string for an empty input."


def test_reverse_complement_invalid_character():
    """Test reverse complement with a non-standard character."""
    seq = "ATXG"
    with pytest.raises(KeyError):
        SplicingPrediction.reverse_complement(seq)


# Mocking constants used in the method
CHROM_REFSEQ_37 = {"1": "NC_000001.10"}
CHROM_REFSEQ_38 = {"1": "NC_000001.11"}


@pytest.fixture
@patch("src.utils.SeqRepo")
def splicing_prediction(mock_seqrepo, seqvar_ss):
    """Fixture for initializing SplicingPrediction with a mocked SeqRepo."""
    mock_seqrepo.return_value = MagicMock()

    return SplicingPrediction(
        seqvar=seqvar_ss,
        strand=GenomicStrand.Plus,
        exons=[],
    )


@patch("src.utils.SplicingPrediction.reverse_complement")
def test_get_sequence_valid_grch37(mock_reverse_complement, splicing_prediction, seqvar_ss):
    """Test sequence retrieval with valid GRCh37 genome release."""
    seqvar_ss.genome_release = GenomeRelease.GRCh37
    splicing_prediction.sr = {CHROM_REFSEQ_37["1"]: "ACGTACGTACGT"}

    result = splicing_prediction.get_sequence(1, 5)
    assert result == "CGTA", "The sequence should be correctly retrieved for GRCh37."


@patch("src.utils.SplicingPrediction.reverse_complement")
def test_get_sequence_valid_grch38(mock_reverse_complement, splicing_prediction, seqvar_ss):
    """Test sequence retrieval with valid GRCh38 genome release."""
    splicing_prediction.sr = {CHROM_REFSEQ_38["1"]: "TGCATGCATGCA"}

    result = splicing_prediction.get_sequence(1, 5)
    assert result == "GCAT", "The sequence should be correctly retrieved for GRCh38."


def test_get_sequence_invalid_genome_release(splicing_prediction, seqvar_ss):
    """Test sequence retrieval with invalid genome release."""
    seqvar_ss.genome_release = "InvalidRelease"

    with pytest.raises(AlgorithmError) as exc_info:
        splicing_prediction.get_sequence(1, 5)
    assert "Invalid genome release" in str(
        exc_info.value
    ), "Should raise an error for invalid genome release."


@patch("src.utils.SplicingPrediction.reverse_complement")
def test_get_sequence_failure(mock_reverse_complement, splicing_prediction, seqvar_ss):
    """Test sequence retrieval failure."""
    splicing_prediction.sr = MagicMock()
    splicing_prediction.sr.__getitem__.side_effect = Exception("SeqRepo error")

    with pytest.raises(AlgorithmError) as exc_info:
        splicing_prediction.get_sequence(1, 5)
    assert "Failed to get sequence" in str(
        exc_info.value
    ), "Should raise an error if sequence retrieval fails."


@patch("src.utils.SplicingPrediction.reverse_complement", return_value="TACG")
def test_get_sequence_minus_strand(mock_reverse_complement, splicing_prediction, seqvar_ss):
    """Test sequence retrieval with reverse complement on minus strand."""
    splicing_prediction.strand = GenomicStrand.Minus
    splicing_prediction.sr = {CHROM_REFSEQ_38["1"]: "ACGTACGTACGT"}

    result = splicing_prediction.get_sequence(1, 5)
    mock_reverse_complement.assert_called_once_with("CGTA")
    assert (
        result == "TACG"
    ), "The sequence should be the reverse complement when the strand is minus."


def test_get_cryptic_ss_unknown_splice_type(splicing_prediction):
    """Test cryptic splice site retrieval when splice type is unknown."""
    result = splicing_prediction.get_cryptic_ss("ACGTACGTACGT", SpliceType.Unknown)
    assert result == [], "Should return an empty list for unknown splice type."


@patch.object(SplicingPrediction, "_find_cryptic_donor_sites", return_value=[(105, "GT", 4.5)])
def test_get_cryptic_ss_donor(mock_find_cryptic_donor_sites, splicing_prediction):
    """Test cryptic splice site retrieval for donor splice type."""
    splicing_prediction.seqvar.pos = 100
    splicing_prediction.maxentscore_ref = 5.0
    refseq = "ACGTACGTACGT"

    result = splicing_prediction.get_cryptic_ss(refseq, SpliceType.Donor)

    mock_find_cryptic_donor_sites.assert_called()
    # Since the _find_cryptic_donor_sites method is mocked, the result contains many repeated sites
    assert ((105, "GT", 4.5)) in result, "Should return the correct cryptic donor sites."


@patch.object(SplicingPrediction, "_find_cryptic_acceptor_sites", return_value=[(95, "AG", 5.2)])
def test_get_cryptic_ss_acceptor(mock_find_cryptic_acceptor_sites, splicing_prediction):
    """Test cryptic splice site retrieval for acceptor splice type."""
    splicing_prediction.seqvar.pos = 100
    splicing_prediction.maxentscore_ref = 5.0
    refseq = "ACGTACGTACGT"

    result = splicing_prediction.get_cryptic_ss(refseq, SpliceType.Acceptor)

    mock_find_cryptic_acceptor_sites.assert_called()
    # Since the _find_cryptic_acceptor_sites method is mocked, the result contains repeated sites
    assert ((95, "AG", 5.2)) in result, "Should return the correct cryptic acceptor sites."


@patch.object(SplicingPrediction, "_find_cryptic_donor_sites", return_value=[])
@patch.object(SplicingPrediction, "_find_cryptic_acceptor_sites", return_value=[])
def test_get_cryptic_ss_no_sites_found(
    mock_find_cryptic_acceptor_sites, mock_find_cryptic_donor_sites, splicing_prediction
):
    """Test cryptic splice site retrieval when no sites are found."""
    splicing_prediction.seqvar.pos = 100
    splicing_prediction.maxentscore_ref = 5.0
    refseq = "ACGTACGTACGT"

    result = splicing_prediction.get_cryptic_ss(refseq, SpliceType.Donor)

    mock_find_cryptic_donor_sites.assert_called()
    assert result == [], "Should return an empty list if no cryptic sites are found."

    result = splicing_prediction.get_cryptic_ss(refseq, SpliceType.Acceptor)

    mock_find_cryptic_acceptor_sites.assert_called()
    assert result == [], "Should return an empty list if no cryptic sites are found."


@patch.object(SplicingPrediction, "get_sequence", return_value="ACGTGTACG")
@patch("src.utils.maxent.score5", return_value=4.5)
def test_find_cryptic_donor_sites_plus_strand(mock_score5, mock_get_sequence, splicing_prediction):
    """Test finding cryptic donor sites on the Plus strand with a valid context."""
    refseq = "ACGTGTACG"
    pos = 100
    refscore = 3.0

    result = splicing_prediction._find_cryptic_donor_sites(pos, refseq, refscore)

    mock_get_sequence.assert_called_once_with(pos, pos + 9)
    mock_score5.assert_called_once_with("ACGTGTACG", matrix=splicing_prediction.matrix5)

    assert result == [(pos, "ACGTGTACG", 4.5)], "Should return the correct cryptic donor site."


@patch.object(SplicingPrediction, "get_sequence", return_value="ACGTGTACG")
@patch("src.utils.maxent.score5", return_value=4.5)
@patch.object(SplicingPrediction, "reverse_complement", return_value="CGTACGTAC")
def test_find_cryptic_donor_sites_minus_strand(
    mock_reverse_complement, mock_score5, mock_get_sequence, splicing_prediction
):
    """Test finding cryptic donor sites on the Minus strand with a valid context."""
    splicing_prediction.strand = GenomicStrand.Minus
    refseq = "CGTACGTAC"
    pos = 100
    refscore = 3.0

    result = splicing_prediction._find_cryptic_donor_sites(pos, refseq, refscore)

    mock_get_sequence.assert_called_once_with(pos, pos + 9)
    mock_reverse_complement.assert_called_once_with("ACGTGTACG")
    mock_score5.assert_called_once_with("CGTACGTAC", matrix=splicing_prediction.matrix5)

    assert result == [
        (pos, "CGTACGTAC", 4.5)
    ], "Should return the correct cryptic donor site for the Minus strand."


@patch.object(SplicingPrediction, "get_sequence", return_value="ACGTGTACG")
@patch("src.utils.maxent.score5", return_value=0.8)
def test_find_cryptic_donor_sites_below_threshold(
    mock_score5, mock_get_sequence, splicing_prediction
):
    """Test when maxent score is below the threshold."""
    refseq = "ACGTGTACG"
    pos = 100
    refscore = 3.0

    result = splicing_prediction._find_cryptic_donor_sites(pos, refseq, refscore)

    mock_get_sequence.assert_called_once_with(pos, pos + 9)
    mock_score5.assert_called_once_with("ACGTGTACG", matrix=splicing_prediction.matrix5)

    assert result == [], "Should return an empty list when the maxent score is below the threshold."


@patch.object(SplicingPrediction, "get_sequence", return_value="ACGTACTAC")
@patch("src.utils.maxent.score5", return_value=4.5)
def test_find_cryptic_donor_sites_no_gt_in_splice_context(
    mock_score5, mock_get_sequence, splicing_prediction
):
    """Test when the splice context does not contain 'GT'."""
    refseq = "ACGAAATAC"
    pos = 100
    refscore = 3.0

    result = splicing_prediction._find_cryptic_donor_sites(pos, refseq, refscore)

    mock_get_sequence.assert_called_once_with(pos, pos + 9)
    mock_score5.assert_called_once_with("ACGTACTAC", matrix=splicing_prediction.matrix5)

    assert (
        result == []
    ), "Should return an empty list when the splice context does not contain 'GT'."


@patch.object(SplicingPrediction, "get_sequence", return_value="ACGTACGTAGTACGTACGTAGTAC")
@patch("src.utils.maxent.score3", return_value=4.5)
def test_find_cryptic_acceptor_sites_plus_strand(
    mock_score3, mock_get_sequence, splicing_prediction
):
    """Test finding cryptic acceptor sites on the Plus strand with a valid context."""
    refseq = "ACGTACGTAGTACGTACGTAGTAC"
    pos = 100
    refscore = 3.0

    result = splicing_prediction._find_cryptic_acceptor_sites(pos, refseq, refscore)

    mock_get_sequence.assert_called_once_with(pos, pos + 23)
    mock_score3.assert_called_once_with(
        "ACGTACGTAGTACGTACGTAGTAC", matrix=splicing_prediction.matrix3
    )

    assert result == [
        (pos, "ACGTACGTAGTACGTACGTAGTAC", 4.5)
    ], "Should return the correct cryptic acceptor site."


@patch.object(SplicingPrediction, "get_sequence", return_value="ACGTACGTAGTACGTACGTAGTAC")
@patch("src.utils.maxent.score3", return_value=4.5)
@patch.object(SplicingPrediction, "reverse_complement", return_value="GTACGTACGTACGTACTACGTACG")
def test_find_cryptic_acceptor_sites_minus_strand(
    mock_reverse_complement, mock_score3, mock_get_sequence, splicing_prediction
):
    """Test finding cryptic acceptor sites on the Minus strand with a valid context."""
    splicing_prediction.strand = GenomicStrand.Minus
    refseq = "GTACGTACGTACGTACTACGTACG"
    pos = 100
    refscore = 3.0

    result = splicing_prediction._find_cryptic_acceptor_sites(pos, refseq, refscore)

    mock_get_sequence.assert_called_once_with(pos, pos + 23)
    mock_reverse_complement.assert_called_once_with("ACGTACGTAGTACGTACGTAGTAC")
    mock_score3.assert_called_once_with(
        "GTACGTACGTACGTACTACGTACG", matrix=splicing_prediction.matrix3
    )

    assert result == [
        (pos, "GTACGTACGTACGTACTACGTACG", 4.5)
    ], "Should return the correct cryptic acceptor site for the Minus strand."


@patch.object(SplicingPrediction, "get_sequence", return_value="ACGTACGTAGTACGTACGTAGTAC")
@patch("src.utils.maxent.score3", return_value=0.8)
def test_find_cryptic_acceptor_sites_below_threshold(
    mock_score3, mock_get_sequence, splicing_prediction
):
    """Test when maxent score is below the threshold."""
    refseq = "ACGTACGTAGTACGTACGTAGTAC"
    pos = 100
    refscore = 3.0

    result = splicing_prediction._find_cryptic_acceptor_sites(pos, refseq, refscore)

    mock_get_sequence.assert_called_once_with(pos, pos + 23)
    mock_score3.assert_called_once_with(
        "ACGTACGTAGTACGTACGTAGTAC", matrix=splicing_prediction.matrix3
    )

    assert result == [], "Should return an empty list when the maxent score is below the threshold."


@patch.object(SplicingPrediction, "get_sequence", return_value="ACGTACGTAGTACGTACGTACTAC")
@patch("src.utils.maxent.score3", return_value=4.5)
def test_find_cryptic_acceptor_sites_no_ag_in_splice_context(
    mock_score3, mock_get_sequence, splicing_prediction
):
    """Test when the splice context does not contain 'AG'."""
    refseq = "ACGTACGTAGTACGTACGAAAAAC"
    pos = 100
    refscore = 3.0

    result = splicing_prediction._find_cryptic_acceptor_sites(pos, refseq, refscore)

    mock_get_sequence.assert_called_once_with(pos, pos + 23)
    mock_score3.assert_called_once_with(
        "ACGTACGTAGTACGTACGTACTAC", matrix=splicing_prediction.matrix3
    )

    assert (
        result == []
    ), "Should return an empty list when the splice context does not contain 'AG'."


# === SeqVarTranscriptsHelper ===


def test_get_ts_info_seqvar_success(ts_helper):
    """Test get_ts_info method with a successful response."""
    # Mock the actual data that would be returned from the Mehari API
    ts_helper.seqvar_ts_info = TranscriptsSeqVar.model_validate(
        get_json_object("mehari/DCDC2_seqvar.json")
    )
    ts_helper.seqvar_transcript = TranscriptsSeqVar.model_validate(
        get_json_object("mehari/DCDC2_seqvar.json")
    ).result
    ts_helper.gene_ts_info = GeneTranscripts.model_validate(get_json_object("mehari/HAL_gene.json"))
    ts_helper.gene_transcript = GeneTranscripts.model_validate(
        get_json_object("mehari/HAL_gene.json")
    ).transcripts
    ts_helper.consequence = SeqVarPVS1Consequence.InitiationCodon

    seqvar_transcript, gene_transcript, seqvar_ts_info, gene_ts_info, consequence = (
        ts_helper.get_ts_info()
    )

    assert seqvar_transcript is not None
    assert gene_transcript is not None
    assert seqvar_ts_info is not None
    assert gene_ts_info is not None
    assert consequence == SeqVarPVS1Consequence.InitiationCodon


def test_get_ts_info_seqvar_failure(ts_helper):
    """Test get_ts_info method with a failed response."""
    seqvar_transcript, gene_transcript, seqvar_ts_info, gene_ts_info, consequence = (
        ts_helper.get_ts_info()
    )

    assert seqvar_transcript is None
    assert gene_transcript is None
    assert seqvar_ts_info == []
    assert gene_ts_info == []
    assert consequence == SeqVarPVS1Consequence.NotSet


@patch.object(MehariClient, "get_seqvar_transcripts")
@patch.object(MehariClient, "get_gene_transcripts")
def test_initialize_seqvar_success(
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
def test_initialize_seqvar_failure(
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
    assert ts_helper.HGNC_id == ""
    assert len(ts_helper.HGVSs) == 0


@patch.object(MehariClient, "get_seqvar_transcripts")
@patch.object(MehariClient, "get_gene_transcripts")
def test_initialize_seqvar_no_transcripts(
    mock_get_gene_transcripts, mock_get_seqvar_transcripts, ts_helper
):
    # Mock the MehariClient methods for CI
    _mock_get_gene_transcripts, _mock_get_seqvar_transcripts = MagicMock(), MagicMock()
    ts_helper.initialize()
    assert ts_helper.gene_ts_info == []
    assert ts_helper.HGNC_id == ""
    assert len(ts_helper.HGVSs) == 0


@pytest.mark.parametrize(
    "consequence_input, expected_consequence",
    [
        (["splice_region_variant"], SeqVarPVS1Consequence.SpliceSites),
        (["splice_donor_variant"], SeqVarPVS1Consequence.SpliceSites),
        (["frameshift_variant"], SeqVarPVS1Consequence.NonsenseFrameshift),
        (["initiator_codon_variant"], SeqVarPVS1Consequence.InitiationCodon),
        (["unknown_consequence"], SeqVarPVS1Consequence.NotSet),
        (["regulatory_region_amplification"], SeqVarPVS1Consequence.NotSet),
        ([""], SeqVarPVS1Consequence.NotSet),
        ([], SeqVarPVS1Consequence.NotSet),
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
    assert consequence == SeqVarPVS1Consequence.NotSet


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
            "mehari/LARP7_gene.json",
            "mehari/LARP7_seqvar.json",
            "NM_001267039.2",
        ),
        (
            ["NM_001267039.2", "NM_001370974.1"],
            "mehari/LARP7_gene.json",
            "mehari/LARP7_seqvar.json",
            "NM_001267039.2",
        ),
    ],
)
def test_choose_transcript_seqvar_success(
    hgvss, gene_ts_file, seqvar_ts_file, expected_hgvs, ts_helper
):
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
        (["invalid"], "mehari/LARP7_gene.json", "mehari/LARP7_seqvar.json"),
        ([], "mehari/LARP7_gene.json", "mehari/LARP7_seqvar.json"),
    ],
)
def test_choose_transcript_seqvar_invalid(hgvss, gene_ts_file, seqvar_ts_file, ts_helper):
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
    assert seqvar_ts is None


# === StrucVarTranscriptsHelper ===


@pytest.fixture
def strucvar():
    return StrucVar(
        sv_type=StrucVarType.DEL,
        genome_release=GenomeRelease.GRCh38,
        chrom="1",
        start=1000,
        stop=2000,
        user_repr="1:1000-2000DEL",
    )


@pytest.fixture
def strucvar_ts_helper(strucvar):
    return StrucVarTranscriptsHelper(strucvar)


@pytest.fixture
def strucvar_transcripts(file_name: str = "mehari/KRTAP9-2_strucvar.json"):
    return TranscriptsStrucVar.model_validate(get_json_object(file_name)).result


@pytest.fixture
def gene_transcripts_strucvar():
    # Simulate a list of TranscriptGene objects with varying lengths and MANE select tags
    return [
        MagicMock(
            id="T1",
            tags=["TRANSCRIPT_TAG_MANE_SELECT"],
            genomeAlignments=[MagicMock(exons=[MagicMock(altStartI=100, altEndI=200)])],
        ),
        MagicMock(
            id="T2",
            tags=[],
            genomeAlignments=[MagicMock(exons=[MagicMock(altStartI=100, altEndI=300)])],
        ),
        MagicMock(
            id="T3",
            tags=["TRANSCRIPT_TAG_MANE_SELECT"],
            genomeAlignments=[MagicMock(exons=[MagicMock(altStartI=100, altEndI=250)])],
        ),
        MagicMock(
            id="T4",
            tags=[],
            genomeAlignments=[MagicMock(exons=[MagicMock(altStartI=100, altEndI=500)])],
        ),
    ]


def test_get_ts_info_strucvar_success(strucvar_ts_helper):
    """Test get_ts_info method with a successful response."""
    # Mock the actual data that would be returned from the Mehari API
    strucvar_ts_helper.strucvar_ts_info = TranscriptsStrucVar.model_validate(
        get_json_object("mehari/KRTAP9-2_strucvar.json")
    )
    strucvar_ts_helper.strucvar_transcript = TranscriptsStrucVar.model_validate(
        get_json_object("mehari/KRTAP9-2_strucvar.json")
    ).result[0]  # assuming the first result is the most relevant
    strucvar_ts_helper.gene_ts_info = GeneTranscripts.model_validate(
        get_json_object("mehari/HAL_gene.json")
    )
    strucvar_ts_helper.gene_transcript = GeneTranscripts.model_validate(
        get_json_object("mehari/HAL_gene.json")
    ).transcripts[0]  # assuming the first transcript is the most relevant

    strucvar_transcript, gene_transcript, strucvar_ts_info, gene_ts_info = (
        strucvar_ts_helper.get_ts_info()
    )

    assert strucvar_transcript is not None
    assert gene_transcript is not None
    assert strucvar_ts_info is not None
    assert gene_ts_info is not None


def test_get_ts_info_strucvar_failure(strucvar_ts_helper):
    """Test get_ts_info method with a failed response."""
    strucvar_transcript, gene_transcript, strucvar_ts_info, gene_ts_info = (
        strucvar_ts_helper.get_ts_info()
    )

    assert strucvar_transcript is None
    assert gene_transcript is None
    assert strucvar_ts_info == []
    assert gene_ts_info == []


@patch.object(MehariClient, "get_strucvar_transcripts")
@patch.object(MehariClient, "get_gene_transcripts")
def test_initialize_strucvar_success(
    mock_get_gene_transcripts,
    mock_get_strucvar_transcripts,
    strucvar,
    strucvar_ts_helper,
    strucvar_transcripts,
    gene_transcripts,
):
    """Test initialize method with a successful response."""
    # Mock successful responses
    mock_get_strucvar_transcripts.return_value = MagicMock(result=strucvar_transcripts)
    mock_get_gene_transcripts.return_value = MagicMock(transcripts=gene_transcripts)

    strucvar_ts_helper.strucvar = strucvar
    strucvar_ts_helper.initialize()

    assert strucvar_ts_helper.strucvar_ts_info is strucvar_transcripts
    assert strucvar_ts_helper.gene_ts_info is gene_transcripts


@patch.object(MehariClient, "get_strucvar_transcripts")
@patch.object(MehariClient, "get_gene_transcripts")
def test_initialize_strucvar_failure(
    mock_get_gene_transcripts,
    mock_get_strucvar_transcripts,
    strucvar,
    strucvar_ts_helper,
):
    """Test initialize method with a failed response."""
    # Mock failed responses
    mock_get_strucvar_transcripts.return_value = None
    mock_get_gene_transcripts.return_value = None

    strucvar_ts_helper.strucvar = strucvar
    strucvar_ts_helper.initialize()

    assert strucvar_ts_helper.strucvar_ts_info == []
    assert strucvar_ts_helper.gene_ts_info == []


def test_choose_transcript_mane_select(gene_transcripts_strucvar):
    """Test that the MANE Select tagged transcript is chosen when available."""
    chosen_transcript = StrucVarTranscriptsHelper._choose_transcript(gene_transcripts_strucvar)
    assert chosen_transcript is not None
    assert chosen_transcript.id == "T3", "The longest MANE transcript should be chosen."


def test_choose_transcript_longest_when_no_mane(gene_transcripts_strucvar):
    """Test that the longest transcript is chosen when no MANE Select tag is present."""
    # Remove MANE select tags for this test
    for transcript in gene_transcripts_strucvar:
        transcript.tags = []
    chosen_transcript = StrucVarTranscriptsHelper._choose_transcript(gene_transcripts_strucvar)
    assert chosen_transcript is not None
    assert (
        chosen_transcript.id == "T4"
    ), "The longest transcript should be chosen when no MANE transcript is available."


def test_choose_transcript_none_available():
    """Test the behavior when no transcripts are available."""
    chosen_transcript = StrucVarTranscriptsHelper._choose_transcript([])
    assert chosen_transcript is None, "Should return None when no transcripts are available."


def test_choose_transcript_exon_length_priority(gene_transcripts_strucvar):
    """Test that exon length is considered for MANE transcripts with equal priority."""
    # Adjusting to create equal MANE priority but different lengths
    gene_transcripts_strucvar[0].genomeAlignments[0].exons[0] = MagicMock(
        altStartI=100, altEndI=300
    )  # Longer exon for T1
    gene_transcripts_strucvar[2].genomeAlignments[0].exons[0] = MagicMock(
        altStartI=100, altEndI=250
    )  # Shorter exon for T3
    chosen_transcript = StrucVarTranscriptsHelper._choose_transcript(gene_transcripts_strucvar)
    assert chosen_transcript is not None
    assert chosen_transcript.id == "T1", "The MANE transcript with longer exons should be chosen."
