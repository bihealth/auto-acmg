from unittest.mock import MagicMock, patch

import pytest

from src.api.annonars import AnnonarsClient
from src.api.mehari import MehariClient
from src.criteria.auto_pvs1 import AutoPVS1, SeqVarPVS1Helper
from src.defs.annonars_range import AnnonarsRangeResponse
from src.defs.auto_acmg import AutoACMGPrediction, AutoACMGStrength, GenomicStrand
from src.defs.auto_pvs1 import (
    PVS1Prediction,
    PVS1PredictionSeqVarPath,
    SeqvarConsequenceMapping,
    SeqVarPVS1Consequence,
)
from src.defs.exceptions import AlgorithmError, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.mehari import Exon, GeneTranscripts, TranscriptsSeqVar
from src.defs.seqvar import SeqVar
from src.utils import SeqVarTranscriptsHelper, SplicingPrediction
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


# === SeqVarPVS1Helper ===


@pytest.mark.parametrize(
    "var_pos,exons,strand,expected_result",
    [
        (
            100,
            [MockExon(0, 100, 0, 100)],
            GenomicStrand.Plus,
            (100, 100),
        ),
        (100, [MockExon(0, 100, 0, 100)], GenomicStrand.Plus, (100, 100)),
        (
            150,
            [MockExon(0, 100, 0, 100), MockExon(100, 200, 100, 200), MockExon(200, 300, 200, 300)],
            GenomicStrand.Plus,
            (150, 300),
        ),
        (
            150,
            [MockExon(0, 100, 0, 100), MockExon(100, 200, 100, 200), MockExon(200, 300, 200, 300)],
            GenomicStrand.Plus,
            (150, 300),
        ),
    ],
)
def test_calc_alt_reg(var_pos, exons, strand, expected_result):
    """Test the _calc_alt_reg method."""
    result = SeqVarPVS1Helper()._calc_alt_reg(var_pos, exons, strand)
    assert result == expected_result


@pytest.mark.parametrize(
    "gene_transcripts_file,transcript_id,var_pos,expected_result",
    [
        (
            "mehari/HAL_gene.json",
            "NM_002108.4",
            96370184,
            (96366439, 96370184),
        ),  # Strand minus
        (
            "mehari/F10_gene_NM_000504.4.json",
            "NM_000504.4",
            113139456,
            (113139456, 113149529),
        ),  # Strand plus
        (
            "mehari/PCID2_gene.json",
            "NM_001127202.4",
            113184385,
            (113177535, 113184385),
        ),  # Strand minus
    ],
)
def test_calc_alt_reg_real_data(gene_transcripts_file, transcript_id, var_pos, expected_result):
    """Test the _calc_alt_reg method."""
    gene_transcripts = GeneTranscripts.model_validate(
        get_json_object(gene_transcripts_file)
    ).transcripts
    tsx = None
    for transcript in gene_transcripts:
        if transcript.id == transcript_id:
            tsx = transcript
    if tsx is None:  # Should never happen
        raise ValueError(f"Transcript {transcript_id} not found in the gene transcripts")

    exons = tsx.genomeAlignments[0].exons
    strand = GenomicStrand.from_string(tsx.genomeAlignments[0].strand)
    start_pos, end_pos = SeqVarPVS1Helper()._calc_alt_reg(var_pos, exons, strand)
    assert start_pos == expected_result[0]
    assert end_pos == expected_result[1]


@pytest.mark.parametrize(
    "annonars_range_response, expected_result",
    [
        ("annonars/GAA_range.json", (535, 2606)),
        ("annonars/CDH1_range.json", (0, 2)),
    ],
)
def test_count_pathogenic_vars(annonars_range_response, expected_result, seqvar):
    """Test the _count_pathogenic_vars method."""
    with patch.object(AnnonarsClient, "_get_variant_from_range") as mock_get_variant_from_range:
        mock_get_variant_from_range.return_value = AnnonarsRangeResponse.model_validate(
            get_json_object(annonars_range_response)
        )
        result = SeqVarPVS1Helper()._count_pathogenic_vars(seqvar, 1, 1000)  # Real range is mocked
        assert result == expected_result


@pytest.mark.parametrize(
    "var_pos,exons,expected_result",
    [
        (100, [MockExon(0, 100, 0, 100)], (0, 100)),
        (150, [MockExon(0, 100, 0, 100), MockExon(100, 200, 100, 200)], (100, 200)),
        (150, [MockExon(0, 100, 0, 100), MockExon(100, 200, 100, 200)], (100, 200)),
    ],
)
def test_find_aff_exon_pos(var_pos, exons, expected_result):
    """Test the _find_aff_exon_pos method."""
    result = SeqVarPVS1Helper()._find_aff_exon_pos(var_pos, exons)
    assert result == expected_result


@pytest.mark.parametrize(
    "value,expected_result",
    [
        (
            SeqVarPVS1Consequence.NonsenseFrameshift,
            [
                "3_prime_utr_variant",
                "3_prime_UTR_variant",
                "frameshift_variant",
                "stop_gained",
            ],
        ),
        (
            SeqVarPVS1Consequence.InitiationCodon,
            [
                "downstream_gene_variant",
                "start_lost",
                "initiator_codon_variant",
                "start_retained_variant",
            ],
        ),
        (
            SeqVarPVS1Consequence.SpliceSites,
            [
                "splice_region_variant",
                "splice_donor_variant",
                "splice_donor_5th_base_variant",
                "splice_donor_region_variant",
                "splice_polypyrimidine_tract_variant",
                "splice_acceptor_variant",
            ],
        ),
        (SeqVarPVS1Consequence.Missense, ["missense_variant"]),
    ],
)
def test_get_conseq(value, expected_result):
    """Test the _get_conseq method."""
    result = SeqVarPVS1Helper()._get_conseq(value)
    assert result == expected_result


@pytest.mark.parametrize(
    "annonars_range_response, expected_result",
    [
        ("annonars/GAA_range.json", (56, 158)),
        ("annonars/CDH1_range.json", (0, 0)),
    ],
)
def test_count_lof_vars(annonars_range_response, expected_result, seqvar):
    """Test the _count_lof_vars method."""
    with patch.object(AnnonarsClient, "_get_variant_from_range") as mock_get_variant_from_range:
        mock_get_variant_from_range.return_value = AnnonarsRangeResponse.model_validate(
            get_json_object(annonars_range_response)
        )
        result = SeqVarPVS1Helper()._count_lof_vars(seqvar, 1, 1000)  # Real range is mocked
        assert result == expected_result


@pytest.mark.parametrize(
    "hgvs, cds_info, expected_result",
    [
        # Test no alternative start codon is found
        (
            "NM_000001",
            {
                "NM_000001": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000002": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
            },
            None,
        ),
        # Test an alternative start codon is found
        (
            "NM_000001",
            {
                "NM_000001": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000002": MockCdsInfo(
                    start_codon=150, stop_codon=1000, cds_start=150, cds_end=1000, exons=[]
                ),
            },
            150,
        ),
        # Test multiple transcripts, one with an alternative start
        (
            "NM_000001",
            {
                "NM_000001": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000002": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000003": MockCdsInfo(
                    start_codon=200, stop_codon=1000, cds_start=200, cds_end=1000, exons=[]
                ),
            },
            200,
        ),
        # Test multiple transcripts, none with an alternative start
        (
            "NM_000001",
            {
                "NM_000001": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000002": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000003": MockCdsInfo(
                    start_codon=100, stop_codon=2000, cds_start=100, cds_end=1000, exons=[]
                ),
            },
            None,
        ),
    ],
)
def test_closest_alt_start_cdn(hgvs, cds_info, expected_result):
    """Test the _closest_alt_start_cdn method."""
    result = SeqVarPVS1Helper()._closest_alt_start_cdn(cds_info, hgvs)
    assert result == expected_result


def test_closest_alt_start_cdn_invalid():
    """Test the _closest_alt_start_cdn method."""
    hgvs = "NM_000"
    cds_info = {
        "NM_000001": MockCdsInfo(
            start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
        ),
        "NM_000002": MockCdsInfo(
            start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
        ),
        "NM_000003": MockCdsInfo(
            start_codon=200, stop_codon=1000, cds_start=200, cds_end=1000, exons=[]
        ),
    }
    with pytest.raises(MissingDataError):
        SeqVarPVS1Helper()._closest_alt_start_cdn(cds_info, hgvs)  # type: ignore


@pytest.mark.parametrize(
    "seqvar, exons, expected_result",
    [
        (
            SeqVar(GenomeRelease.GRCh38, "1", 50, "A", "T"),
            [MockExon(2, 100, 0, 100)],
            (2, 100),
        ),  # One exon
        (
            SeqVar(GenomeRelease.GRCh38, "1", 150, "A", "T"),
            [MockExon(0, 100, 0, 100), MockExon(100, 200, 100, 200)],
            (100, 200),
        ),  # Multiple exons
        (
            SeqVar(GenomeRelease.GRCh38, "1", 98, "A", "T"),
            [MockExon(100, 200, 100, 200)],
            (100, 200),
        ),  # Upstream intron
        (
            SeqVar(GenomeRelease.GRCh38, "1", 220, "A", "T"),
            [MockExon(100, 200, 100, 200)],
            (100, 200),
        ),  # Downstream intron
    ],
)
def test_skipping_exon_pos(seqvar, exons, expected_result):
    """Test the _skipping_exon_pos method."""
    result = SeqVarPVS1Helper()._skipping_exon_pos(seqvar, exons)
    assert result == expected_result


def test_skipping_exon_pos_invalid():
    """Test the _skipping_exon_pos method."""
    seqvar = SeqVar(GenomeRelease.GRCh38, "1", 50, "A", "T")
    exons = [MockExon(100, 200, 100, 200)]
    with pytest.raises(AlgorithmError):
        SeqVarPVS1Helper()._skipping_exon_pos(seqvar, exons)  # type: ignore


@pytest.mark.parametrize(
    "gene_transcripts_file,transcript_id,hgnc_id,var_pos,expected_result",
    [
        (
            "mehari/F10_gene_NM_000504.4.json",
            "NM_000504.4",
            "HGNC:3528",
            500,
            True,
        ),  # Strand plus
        (
            "mehari/F10_gene_NM_000504.4.json",
            "NM_000504.4",
            "HGNC:3528",
            1000,
            False,
        ),  # Strand plus
        (
            "mehari/PCID2_gene.json",
            "NM_001127202.4",
            "HGNC:25653",
            900,
            True,
        ),  # Strand minus. Not a frameshift!
        (
            "mehari/PCID2_gene.json",
            "NM_001127202.4",
            "HGNC:25653",
            1100,
            False,
        ),  # Strand minus. Not a real variant
    ],
)
def test_undergo_nmd(gene_transcripts_file, transcript_id, hgnc_id, var_pos, expected_result):
    """
    Test the _undergo_nmd method. Note, that we don't mock the `_get_variant_position` and
    `_calculate_5_prime_UTR_length` methods.
    """
    gene_transcripts = GeneTranscripts.model_validate(
        get_json_object(gene_transcripts_file)
    ).transcripts
    tsx = None
    for transcript in gene_transcripts:
        if transcript.id == transcript_id:
            tsx = transcript
    if tsx is None:  # Should never happen
        raise ValueError(f"Transcript {transcript_id} not found in the gene transcripts")
    exons = tsx.genomeAlignments[0].exons
    strand = GenomicStrand.from_string(tsx.genomeAlignments[0].strand)
    result = SeqVarPVS1Helper().undergo_nmd(var_pos, hgnc_id, strand, exons)
    assert result == expected_result


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
def test_in_bio_relevant_tsx(transcript_tags, expected_result):
    """Test the _in_bio_relevant_tsx method."""
    result = SeqVarPVS1Helper().in_bio_relevant_tx(transcript_tags)
    assert result == expected_result, f"Failed for transcript_tags: {transcript_tags}"


@pytest.mark.parametrize(
    "pathogenic_variants, total_variants, strand, expected_result",
    [
        (6, 100, GenomicStrand.Plus, True),  # Test pathogenic variants exceed the threshold
        (3, 100, GenomicStrand.Plus, False),  # Test pathogenic variants do not exceed the threshold
        (0, 0, GenomicStrand.Plus, False),  # Test no variants are found
        (0, 100, GenomicStrand.Plus, False),  # Test no pathogenic variants are found
        (100, 0, GenomicStrand.Plus, False),  # Test more pathogenic variants than total variants
    ],
)
def test_crit4prot_func(
    seqvar, pathogenic_variants, total_variants, strand, expected_result, monkeypatch
):
    """Test the _crit4prot_func method."""
    # Mock exons, since they won't affect the outcome
    exons = [MagicMock(spec=Exon)]
    # Mocking _calc_alt_reg to return a mocked range
    mock_calculate = MagicMock(return_value=(1, 1000))
    monkeypatch.setattr(SeqVarPVS1Helper, "_calc_alt_reg", mock_calculate)
    # Mocking _count_pathogenic_vars to return controlled counts of pathogenic and total variants
    mock_count_pathogenic = MagicMock(return_value=(pathogenic_variants, total_variants))
    monkeypatch.setattr(SeqVarPVS1Helper, "_count_pathogenic_vars", mock_count_pathogenic)

    result = SeqVarPVS1Helper().crit4prot_func(seqvar, exons, strand)  # type: ignore
    assert result == expected_result


@pytest.mark.parametrize(
    "frequent_lof_variants, lof_variants, strand, expected_result",
    [
        (
            11,
            100,
            GenomicStrand.Plus,
            True,
        ),  # Test case where frequent LoF variants exceed the 10% threshold
        (
            5,
            100,
            GenomicStrand.Plus,
            False,
        ),  # Test case where frequent LoF variants do not exceed the 10% threshold
        (0, 0, GenomicStrand.Plus, False),  # Test case where no LoF variants are found
        (
            0,
            100,
            GenomicStrand.Plus,
            False,
        ),  # Test case where no frequent LoF variants are found
        (0, 0, GenomicStrand.Plus, False),  # Test case where no LoF variants are found
    ],
)
def test_lof_freq_in_pop(
    seqvar, frequent_lof_variants, lof_variants, strand, expected_result, monkeypatch
):
    """Test the _lof_freq_in_pop method."""
    # Mocking exons, since they won't affect the outcome
    exons = [MagicMock(spec=Exon)]
    # Mocking _calc_alt_reg to return a mocked range
    mock_calculate = MagicMock(return_value=(1, 1000))
    monkeypatch.setattr(SeqVarPVS1Helper, "_find_aff_exon_pos", mock_calculate)
    # Mocking _count_lof_vars to return controlled counts of frequent and total LoF variants
    mock_count_lof_vars = MagicMock(return_value=(frequent_lof_variants, lof_variants))
    monkeypatch.setattr(SeqVarPVS1Helper, "_count_lof_vars", mock_count_lof_vars)

    result = SeqVarPVS1Helper().lof_freq_in_pop(seqvar, exons, strand)  # type: ignore
    assert result == expected_result


@pytest.mark.parametrize(
    "prot_pos, prot_length, expected_result",
    [
        (
            4,
            100,
            False,
        ),  # Test case where the variant remove less than 10% of the protein
        (
            99,
            100,
            True,
        ),  # Test case where the variant removes more than 10% of the protein
    ],
)
def test_lof_rm_gt_10pct_of_prot(prot_pos, prot_length, expected_result):
    """Test the _lof_rm_gt_10pct_of_prot method."""
    result = SeqVarPVS1Helper().lof_rm_gt_10pct_of_prot(prot_pos, prot_length)
    assert result == expected_result


# ============== exon_skip_or_cryptic_ss_disrupt ==============


@pytest.fixture
def seqvar_ss():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def exons():
    return [MagicMock(start=90, end=110), MagicMock(start=120, end=140)]


@pytest.fixture
def consequences():
    return ["splice_acceptor_variant", "splice_donor_variant"]


@pytest.fixture
def helper():
    return SeqVarPVS1Helper()


@pytest.mark.skip(reason="Patching is not working properly")
@patch("src.criteria.auto_pvs1.SeqVarPVS1Helper._skipping_exon_pos", return_value=(90, 111))
@patch("src.criteria.auto_pvs1.SplicingPrediction.get_sequence", return_value="ATGC" * 10)
@patch("src.criteria.auto_pvs1.SplicingPrediction.determine_splice_type", return_value="donor")
@patch("src.criteria.auto_pvs1.SplicingPrediction.get_cryptic_ss", return_value=[(95, "GT", 5)])
def test_exon_skip_or_cryptic_ss_disrupt_exon_skipping(
    mock_skipping_exon_pos, helper, seqvar_ss, exons, consequences
):
    """Test when the exon length is not a multiple of 3, predicting exon skipping."""
    result = helper.exon_skip_or_cryptic_ss_disrupt(
        seqvar_ss, exons, consequences, GenomicStrand.Plus
    )
    assert result is True
    assert (
        "Exon length is not a multiple of 3. Predicted to cause exon skipping."
        in helper.comment_pvs1
    )


@pytest.mark.skip(reason="Patching is not working properly")
@patch("src.criteria.auto_pvs1.SeqVarPVS1Helper._skipping_exon_pos", return_value=(90, 110))
@patch("src.criteria.auto_pvs1.SplicingPrediction.get_sequence", return_value="ATGC" * 10)
@patch("src.criteria.auto_pvs1.SplicingPrediction.determine_splice_type", return_value="donor")
@patch("src.criteria.auto_pvs1.SplicingPrediction.get_cryptic_ss", return_value=[(95, "GT", 5)])
def test_exon_skip_or_cryptic_ss_disrupt_cryptic_splice_site_disruption(
    mock_get_sequence,
    mock_determine_splice_type,
    mock_get_cryptic_ss,
    mock_skipping_exon_pos,
    helper,
    seqvar_ss,
    exons,
    consequences,
):
    """Test when there is a cryptic splice site disruption."""
    result = helper.exon_skip_or_cryptic_ss_disrupt(
        seqvar_ss, exons, consequences, GenomicStrand.Plus
    )
    assert result is True
    assert "Cryptic splice site disruption predicted." in helper.comment_pvs1


@pytest.mark.skip(reason="Patching is not working properly")
@patch("src.criteria.auto_pvs1.SeqVarPVS1Helper._skipping_exon_pos", return_value=(90, 110))
@patch("src.criteria.auto_pvs1.SplicingPrediction.get_sequence", return_value="ATGC" * 10)
@patch("src.criteria.auto_pvs1.SplicingPrediction.determine_splice_type", return_value="donor")
@patch("src.criteria.auto_pvs1.SplicingPrediction.get_cryptic_ss", return_value=[])
def test_exon_skip_or_cryptic_ss_disrupt_preserve_reading_frame(
    mock_get_sequence,
    mock_determine_splice_type,
    mock_get_cryptic_ss,
    mock_skipping_exon_pos,
    helper,
    seqvar_ss,
    exons,
    consequences,
):
    """Test when there is no cryptic splice site and exon length is a multiple of 3, preserving the reading frame."""
    result = helper.exon_skip_or_cryptic_ss_disrupt(
        seqvar_ss, exons, consequences, GenomicStrand.Plus
    )
    assert result is False
    assert (
        "No cryptic splice site found. Predicted to preserve reading frame." in helper.comment_pvs1
    )


def test_exon_skip_or_cryptic_ss_disrupt_missing_strand(helper, seqvar_ss, exons, consequences):
    """Test when the strand is not set, raising MissingDataError."""
    with pytest.raises(MissingDataError):
        helper.exon_skip_or_cryptic_ss_disrupt(seqvar_ss, exons, consequences, GenomicStrand.NotSet)


@pytest.mark.parametrize(
    "hgvs, cds_info, expected_result",
    [
        # Test no alternative start codon is found
        (
            "NM_000001",
            {
                "NM_000001": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000002": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
            },
            False,
        ),
        # Test an alternative start codon is found
        (
            "NM_000001",
            {
                "NM_000001": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000002": MockCdsInfo(
                    start_codon=150, stop_codon=1000, cds_start=150, cds_end=1000, exons=[]
                ),
            },
            True,
        ),
        # Test multiple transcripts, one with an alternative start
        (
            "NM_000001",
            {
                "NM_000001": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000002": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000003": MockCdsInfo(
                    start_codon=200, stop_codon=1000, cds_start=200, cds_end=1000, exons=[]
                ),
            },
            True,
        ),
        # Test multiple transcripts, none with an alternative start
        (
            "NM_000001",
            {
                "NM_000001": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000002": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000003": MockCdsInfo(
                    start_codon=100, stop_codon=2000, cds_start=100, cds_end=1000, exons=[]
                ),
            },
            False,
        ),
    ],
)
def test_alt_start_cdn(hgvs, cds_info, expected_result):
    """Test the _alt_start_cdn method."""
    result = SeqVarPVS1Helper().alt_start_cdn(cds_info, hgvs)
    assert result == expected_result


@pytest.mark.parametrize(
    "exons, strand, pathogenic_variants, hgvs, cds_info, expected_result",
    [
        (
            [MagicMock(altStartI=1, altEndI=200, altCdsStartI=1, altCdsEndI=200)],
            GenomicStrand.Plus,
            1,
            "NM_000001",
            {
                "NM_000001": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000002": MockCdsInfo(
                    start_codon=150, stop_codon=1000, cds_start=150, cds_end=1000, exons=[]
                ),
            },
            True,
        ),  # Test pathogenic variant is found
        (
            [MagicMock(altStartI=1, altEndI=200, altCdsStartI=1, altCdsEndI=200)],
            GenomicStrand.Plus,
            0,
            "NM_000001",
            {
                "NM_000001": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000002": MockCdsInfo(
                    start_codon=150, stop_codon=1000, cds_start=150, cds_end=1000, exons=[]
                ),
            },
            False,
        ),  # Test no pathogenic variants found
        (
            [MagicMock(altStartI=1, altEndI=200, altCdsStartI=1, altCdsEndI=200)],
            GenomicStrand.Minus,
            1,
            "NM_000001",
            {
                "NM_000001": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000002": MockCdsInfo(
                    start_codon=150, stop_codon=1000, cds_start=150, cds_end=1000, exons=[]
                ),
            },
            True,
        ),  # Test pathogenic variants found on minus strand
        (
            [MagicMock(altStartI=1, altEndI=200, altCdsStartI=1, altCdsEndI=200)],
            GenomicStrand.Minus,
            0,
            "NM_000001",
            {
                "NM_000001": MockCdsInfo(
                    start_codon=100, stop_codon=1000, cds_start=100, cds_end=1000, exons=[]
                ),
                "NM_000002": MockCdsInfo(
                    start_codon=150, stop_codon=1000, cds_start=150, cds_end=1000, exons=[]
                ),
            },
            False,
        ),  # Test no pathogenic variants found on minus strand
    ],
)
def test_up_pathogenic_vars(
    seqvar, exons, strand, pathogenic_variants, hgvs, cds_info, expected_result, monkeypatch
):
    """Test the _up_pathogenic_vars method."""
    # Mocking _count_pathogenic_vars to return a controlled number of pathogenic variants
    mock_count_pathogenic = MagicMock(return_value=(pathogenic_variants, 10))
    monkeypatch.setattr(SeqVarPVS1Helper, "_count_pathogenic_vars", mock_count_pathogenic)

    result = SeqVarPVS1Helper().up_pathogenic_vars(seqvar, exons, strand)
    assert result == expected_result


# ============== AutoPVS1 ==============


@pytest.fixture
def auto_pvs1():
    return AutoPVS1()


@pytest.fixture
def seqvar_a():
    return SeqVar(genome_release=GenomeRelease.GRCh38, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def var_data():
    return MagicMock(
        hgnc_id="HGNC:12345",
        prot_pos=100,
        prot_length=500,
        transcript_tags=["important_tag"],
        consequence=MagicMock(mehari=["splice_acceptor_variant"], cadd=""),
        tx_pos_utr=0,
        exons=[],
        strand=GenomicStrand.Plus,
        transcript_id="transcript_1",
    )


def test_convert_consequence_mapped_mehari(auto_pvs1, var_data):
    """Test _convert_consequence when a Mehari consequence is mapped."""
    var_data.consequence.mehari = ["nonsense"]
    SeqvarConsequenceMapping["nonsense"] = SeqVarPVS1Consequence.NonsenseFrameshift

    result = auto_pvs1._convert_consequence(var_data)

    assert (
        result == SeqVarPVS1Consequence.NonsenseFrameshift
    ), "The consequence should be mapped correctly from Mehari."


def test_convert_consequence_mapped_cadd(auto_pvs1, var_data):
    """Test _convert_consequence when a CADD consequence is mapped."""
    var_data.consequence.cadd = "splice_donor_variant"
    SeqvarConsequenceMapping["splice_donor_variant"] = SeqVarPVS1Consequence.SpliceSites

    result = auto_pvs1._convert_consequence(var_data)

    assert (
        result == SeqVarPVS1Consequence.SpliceSites
    ), "The consequence should be mapped correctly from CADD."


def test_convert_consequence_no_mapping(auto_pvs1, var_data):
    """Test _convert_consequence when no consequence is mapped."""
    var_data.consequence.mehari = ["nonexistent_consequence"]
    var_data.consequence.cadd = "another_nonexistent_consequence"

    result = auto_pvs1._convert_consequence(var_data)

    assert (
        result == SeqVarPVS1Consequence.NotSet
    ), "The consequence should be set to NotSet when no mapping is found."


def test_convert_consequence_mapped_mehari_priority(auto_pvs1, var_data):
    """Test _convert_consequence gives priority to Mehari consequence."""
    var_data.consequence.mehari = ["frameshift_variant"]
    var_data.consequence.cadd = "splice_acceptor_variant"
    SeqvarConsequenceMapping["frameshift_variant"] = SeqVarPVS1Consequence.NonsenseFrameshift
    SeqvarConsequenceMapping["splice_acceptor_variant"] = SeqVarPVS1Consequence.SpliceSites

    result = auto_pvs1._convert_consequence(var_data)

    assert (
        result == SeqVarPVS1Consequence.NonsenseFrameshift
    ), "The Mehari consequence should take priority if present."


@patch.object(AutoPVS1, "undergo_nmd", return_value=True)
@patch.object(AutoPVS1, "in_bio_relevant_tx", return_value=True)
@patch.object(AutoPVS1, "crit4prot_func", return_value=False)
@patch.object(AutoPVS1, "lof_freq_in_pop", return_value=False)
@patch.object(AutoPVS1, "lof_rm_gt_10pct_of_prot", return_value=False)
@patch.object(
    AutoPVS1, "_convert_consequence", return_value=SeqVarPVS1Consequence.NonsenseFrameshift
)
def test_verify_pvs1_nonsense_frameshift(
    mock_convert_consequence,
    mock_lof_rm_gt_10pct_of_prot,
    mock_lof_freq_in_pop,
    mock_crit4prot_func,
    mock_in_bio_relevant_tx,
    mock_undergo_nmd,
    auto_pvs1,
    seqvar_a,
    var_data,
):
    """Test verify_pvs1 for nonsense/frameshift variant prediction."""
    var_data.hgnc_id = "HGNC:9588"
    var_data.prot_pos = 200

    prediction, path, comment = auto_pvs1.verify_pvs1(seqvar_a, var_data)

    assert prediction == PVS1Prediction.PVS1, "Expected PVS1 prediction for PTEN gene."
    assert path == PVS1PredictionSeqVarPath.PTEN, "Expected PTEN prediction path."
    assert "The alteration results in a premature stop codon" in comment


@patch.object(AutoPVS1, "exon_skip_or_cryptic_ss_disrupt", return_value=True)
@patch.object(AutoPVS1, "undergo_nmd", return_value=True)
@patch.object(AutoPVS1, "in_bio_relevant_tx", return_value=True)
@patch.object(AutoPVS1, "_convert_consequence", return_value=SeqVarPVS1Consequence.SpliceSites)
def test_verify_pvs1_splice_sites(
    mock_convert_consequence,
    mock_in_bio_relevant_tx,
    mock_undergo_nmd,
    mock_exon_skip_or_cryptic_ss_disrupt,
    auto_pvs1,
    seqvar_a,
    var_data,
):
    """Test verify_pvs1 for splice site variant prediction."""
    prediction, path, comment = auto_pvs1.verify_pvs1(seqvar_a, var_data)

    assert prediction == PVS1Prediction.PVS1, "Expected PVS1 prediction for splice site variant."
    assert path == PVS1PredictionSeqVarPath.SS1, "Expected SS1 prediction path."
    assert "Analysing as splice site variant" in comment


@patch.object(AutoPVS1, "alt_start_cdn", return_value=True)
@patch.object(AutoPVS1, "_convert_consequence", return_value=SeqVarPVS1Consequence.InitiationCodon)
def test_verify_pvs1_initiation_codon(
    mock_convert_consequence, mock_alt_start_cdn, auto_pvs1, seqvar_a, var_data
):
    """Test verify_pvs1 for initiation codon variant prediction."""
    prediction, path, comment = auto_pvs1.verify_pvs1(seqvar_a, var_data)

    assert prediction == PVS1Prediction.NotPVS1, "Expected NotPVS1 prediction for alt start codon."
    assert path == PVS1PredictionSeqVarPath.IC3, "Expected IC3 prediction path."
    assert "Analysing as initiation codon variant" in comment


@patch.object(AutoPVS1, "_convert_consequence", return_value=SeqVarPVS1Consequence.Missense)
def test_verify_pvs1_unsupported_consequence(
    mock_convert_consequence, auto_pvs1, seqvar_a, var_data
):
    """Test verify_pvs1 for unsupported consequence."""
    prediction, path, comment = auto_pvs1.verify_pvs1(seqvar_a, var_data)

    assert prediction == PVS1Prediction.UnsupportedConsequence, "Expected UnsupportedConsequence."
    assert path == PVS1PredictionSeqVarPath.NotSet, "Expected NotSet prediction path."
    assert "PVS1 criteria cannot be applied" in comment


@patch.object(AutoPVS1, "verify_pvs1")
def test_predict_pvs1_met(mock_verify_pvs1, auto_pvs1, seqvar, var_data):
    """Test predict_pvs1 when PVS1 criteria is met."""
    mock_verify_pvs1.return_value = (
        PVS1Prediction.PVS1,
        PVS1PredictionSeqVarPath.NF1,
        "PVS1 criteria met.",
    )

    result = auto_pvs1.predict_pvs1(seqvar, var_data)

    assert result.name == "PVS1", "Criteria name should be PVS1."
    assert result.prediction == AutoACMGPrediction.Met, "Prediction should be Met."
    assert (
        result.strength == AutoACMGStrength.PathogenicVeryStrong
    ), "Strength should be Very Strong."
    assert "PVS1 criteria met." in result.summary, "Summary should contain the correct comment."
    assert "Prediction path: " in result.description, "Description should contain prediction path."


@patch.object(AutoPVS1, "verify_pvs1")
def test_predict_pvs1_not_met(mock_verify_pvs1, auto_pvs1, seqvar, var_data):
    """Test predict_pvs1 when PVS1 criteria is not met."""
    mock_verify_pvs1.return_value = (
        PVS1Prediction.NotPVS1,
        PVS1PredictionSeqVarPath.NF4,
        "PVS1 criteria not met.",
    )

    result = auto_pvs1.predict_pvs1(seqvar, var_data)

    assert result.name == "PVS1", "Criteria name should be PVS1."
    assert result.prediction == AutoACMGPrediction.NotMet, "Prediction should be Not Met."
    assert (
        result.strength == AutoACMGStrength.PathogenicVeryStrong
    ), "Strength should be Very Strong."
    assert "PVS1 criteria not met." in result.summary, "Summary should contain the correct comment."
    assert "Prediction path: " in result.description, "Description should contain prediction path."


@patch.object(AutoPVS1, "verify_pvs1")
def test_predict_pvs1_strong(mock_verify_pvs1, auto_pvs1, seqvar, var_data):
    """Test predict_pvs1 when PVS1 criteria is met with Strong strength."""
    mock_verify_pvs1.return_value = (
        PVS1Prediction.PVS1_Strong,
        PVS1PredictionSeqVarPath.NF3,
        "PVS1 strong criteria met.",
    )

    result = auto_pvs1.predict_pvs1(seqvar, var_data)

    assert result.name == "PVS1", "Criteria name should be PVS1."
    assert result.prediction == AutoACMGPrediction.Met, "Prediction should be Met."
    assert result.strength == AutoACMGStrength.PathogenicStrong, "Strength should be Strong."
    assert (
        "PVS1 strong criteria met." in result.summary
    ), "Summary should contain the correct comment."
    assert "Prediction path: " in result.description, "Description should contain prediction path."


@patch.object(AutoPVS1, "verify_pvs1")
def test_predict_pvs1_moderate(mock_verify_pvs1, auto_pvs1, seqvar, var_data):
    """Test predict_pvs1 when PVS1 criteria is met with Moderate strength."""
    mock_verify_pvs1.return_value = (
        PVS1Prediction.PVS1_Moderate,
        PVS1PredictionSeqVarPath.SS6,
        "PVS1 moderate criteria met.",
    )

    result = auto_pvs1.predict_pvs1(seqvar, var_data)

    assert result.name == "PVS1", "Criteria name should be PVS1."
    assert result.prediction == AutoACMGPrediction.Met, "Prediction should be Met."
    assert result.strength == AutoACMGStrength.PathogenicModerate, "Strength should be Moderate."
    assert (
        "PVS1 moderate criteria met." in result.summary
    ), "Summary should contain the correct comment."
    assert "Prediction path: " in result.description, "Description should contain prediction path."


@patch.object(AutoPVS1, "verify_pvs1")
def test_predict_pvs1_supporting(mock_verify_pvs1, auto_pvs1, seqvar, var_data):
    """Test predict_pvs1 when PVS1 criteria is met with Supporting strength."""
    mock_verify_pvs1.return_value = (
        PVS1Prediction.PVS1_Supporting,
        PVS1PredictionSeqVarPath.IC1,
        "PVS1 supporting criteria met.",
    )

    result = auto_pvs1.predict_pvs1(seqvar, var_data)

    assert result.name == "PVS1", "Criteria name should be PVS1."
    assert result.prediction == AutoACMGPrediction.Met, "Prediction should be Met."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "Strength should be Supporting."
    assert (
        "PVS1 supporting criteria met." in result.summary
    ), "Summary should contain the correct comment."
    assert "Prediction path: " in result.description, "Description should contain prediction path."
