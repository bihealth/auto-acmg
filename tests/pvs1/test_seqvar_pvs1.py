from unittest.mock import MagicMock, patch

import pytest

from src.api.annonars import AnnonarsClient
from src.api.mehari import MehariClient
from src.defs.annonars_range import AnnonarsRangeResponse
from src.defs.auto_pvs1 import (
    GenomicStrand,
    PVS1Prediction,
    PVS1PredictionSeqVarPath,
    SeqVarConsequence,
)
from src.defs.exceptions import AlgorithmError, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.mehari import Exon, GeneTranscripts, TranscriptsSeqVar
from src.defs.seqvar import SeqVar
from src.pvs1.seqvar_pvs1 import SeqVarPVS1, SeqVarPVS1Helper, SeqVarTranscriptsHelper
from src.utils import SplicingPrediction
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
    "main_hgvs, main_hgvs_p, transcripts_data, expected_result",
    [
        # Case where main transcript has valid protein HGVS
        ("NM_000001.1", "p.Gly100Ser", [], "NM_000001.1:p.Gly100Ser"),
        # Case where main transcript HGVS protein is not set, but another transcript has it
        ("NM_000001.1", "", [("NM_000002.1", "p.Arg200Gln")], "NM_000001.1:p.Arg200Gln"),
        # Case where main transcript and others do not have valid protein HGVS
        ("NM_000001.1", "", [("NM_000002.1", ""), ("NM_000003.1", "p.?")], "NM_000001.1:p.?"),
        # Case with no valid protein HGVS notation in any transcript
        ("NM_000001.1", "p.?", [("NM_000002.1", "p.?"), ("NM_000003.1", "")], "NM_000001.1:p.?"),
        # Case where multiple transcripts have valid HGVS, but the first valid one is chosen
        (
            "NM_000001.1",
            "",
            [("NM_000002.1", ""), ("NM_000003.1", "p.Lys300Thr")],
            "NM_000001.1:p.Lys300Thr",
        ),
    ],
)
def test_choose_hgvs_p(main_hgvs, main_hgvs_p, transcripts_data, expected_result):
    # Mocking the main and other transcripts
    main_transcript = MagicMock(hgvs_p=main_hgvs_p)
    transcripts = [MagicMock(feature_id=id, hgvs_p=hgvs_p) for id, hgvs_p in transcripts_data]

    # Invoke the method under test
    result = SeqVarPVS1Helper._choose_hgvs_p(main_hgvs, main_transcript, transcripts)  # type: ignore

    # Verify the result
    assert result == expected_result


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
        ("'NM_000218.3:p.Y662S'", -1),  # Missense mutation
    ],
)
def test_get_pHGVS_termination(pHGVS, expected_termination):
    """Test the _get_pHGVS_termination method."""
    termination = SeqVarPVS1Helper._get_pHGVS_termination(pHGVS)
    assert termination == expected_termination, f"Failed for pHGVS: {pHGVS}"


@pytest.mark.parametrize(
    "tHGVS,expected_result",
    [
        ("c.*1102=", 1102),  # Stop codon at position 1102
        ("c.2506G>T", 2506),  # Missense mutation at position 2506
        ("c.2506_2507insT", 2506),  # Insertion at position 2506
        ("c.2506_2507del", 2506),  # Deletion at position 2506
        ("c.2506_2507delinsT", 2506),  # Deletion-insertion at position 2506
        ("c.1234+5G>T", 1234),  # Intron mutation (splice site)
        ("c.1234-1G>A", 1234),  # Intron mutation (splice site)
        ("c.234_235del", 234),  # Deletion with range
        ("c.234_235insA", 234),  # Insertion with range
        ("c.234_235delinsA", 234),  # Deletion-insertion with range
        ("c.234+5_234+7del", 234),  # Complex range with intron positions
        ("c.234-5_234-7del", 234),  # Complex range with intron positions
        ("invalid", -1),  # Invalid HGVS
    ],
)
def test_get_cds_position(tHGVS, expected_result):
    """Test the _get_cds_position method."""
    result = SeqVarPVS1Helper()._get_cds_position(tHGVS)
    assert result == expected_result


@pytest.mark.parametrize(
    "gene_transcripts_file,transcript_id,expected_result",
    [
        ("mehari/mehari_genes_success.json", "NM_002108.4", 348),
        ("mehari/f10_mehari_genes_NM_000504.4.json", "NM_000504.4", 57),
        ("mehari/pcid2_mehari_genes.json", "NM_001127202.4", 35),
    ],
)
def test_calculate_5_prime_UTR_length(gene_transcripts_file, transcript_id, expected_result):
    """Test the _calculate_5_prime_UTR_length method."""
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
    cds_start = tsx.genomeAlignments[0].cdsStart
    cds_end = tsx.genomeAlignments[0].cdsEnd
    strand = GenomicStrand.from_string(tsx.genomeAlignments[0].strand)
    result = SeqVarPVS1Helper()._calculate_5_prime_UTR_length(exons, cds_start, cds_end, strand)
    assert result == expected_result


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
def test_calculate_altered_region(var_pos, exons, strand, expected_result):
    """Test the _calculate_altered_region method."""
    result = SeqVarPVS1Helper._calculate_altered_region(var_pos, exons, strand)
    assert result == expected_result


@pytest.mark.parametrize(
    "gene_transcripts_file,transcript_id,var_pos,expected_result",
    [
        (
            "mehari/mehari_genes_success.json",
            "NM_002108.4",
            96370184,
            (96366439, 96370184),
        ),  # Strand minus
        (
            "mehari/f10_mehari_genes_NM_000504.4.json",
            "NM_000504.4",
            113139456,
            (113139456, 113149529),
        ),  # Strand plus
        (
            "mehari/pcid2_mehari_genes.json",
            "NM_001127202.4",
            113184385,
            (113177535, 113184385),
        ),  # Strand minus
    ],
)
def test_calculate_altered_region_real_data(
    gene_transcripts_file, transcript_id, var_pos, expected_result
):
    """Test the _calculate_altered_region method."""
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
    start_pos, end_pos = SeqVarPVS1Helper()._calculate_altered_region(var_pos, exons, strand)
    assert start_pos == expected_result[0]
    assert end_pos == expected_result[1]


@pytest.mark.parametrize(
    "annonars_range_response, expected_result",
    [
        ("annonars/NM_000152.4:c.1A>G_annonars_range.json", (501, 2205)),
        ("annonars/CDH1_range.json", (0, 2)),
    ],
)
def test_count_pathogenic_variants(annonars_range_response, expected_result, seqvar):
    """Test the _count_pathogenic_variants method."""
    with patch.object(AnnonarsClient, "get_variant_from_range") as mock_get_variant_from_range:
        mock_get_variant_from_range.return_value = AnnonarsRangeResponse.model_validate(
            get_json_object(annonars_range_response)
        )
        result = SeqVarPVS1Helper()._count_pathogenic_variants(
            seqvar, 1, 1000
        )  # Real range is mocked
        assert result == expected_result


@pytest.mark.parametrize(
    "value,expected_result",
    [
        (
            SeqVarConsequence.NonsenseFrameshift,
            [
                "3_prime_utr_variant",
                "3_prime_UTR_variant",
                "frameshift_variant",
                "stop_gained",
            ],
        ),
        (
            SeqVarConsequence.InitiationCodon,
            [
                "upstream_gene_variant",
                "downstream_gene_variant",
                "start_lost",
                "initiator_codon_variant",
                "start_retained_variant",
            ],
        ),
        (
            SeqVarConsequence.SpliceSites,
            [
                "splice_region_variant",
                "splice_donor_variant",
                "splice_donor_5th_base_variant",
                "splice_donor_region_variant",
                "splice_polypyrimidine_tract_variant",
                "splice_acceptor_variant",
            ],
        ),
        (SeqVarConsequence.Missense, ["missense_variant"]),
    ],
)
def test_get_consequence(value, expected_result):
    """Test the _get_consequence method."""
    result = SeqVarPVS1Helper._get_consequence(value)
    assert result == expected_result


@pytest.mark.parametrize(
    "annonars_range_response, expected_result",
    [
        ("annonars/NM_000152.4:c.1A>G_annonars_range.json", (56, 158)),
        ("annonars/CDH1_range.json", (0, 0)),
    ],
)
def test_count_lof_variants(annonars_range_response, expected_result, seqvar):
    """Test the _count_pathogenic_variants method."""
    with patch.object(AnnonarsClient, "get_variant_from_range") as mock_get_variant_from_range:
        mock_get_variant_from_range.return_value = AnnonarsRangeResponse.model_validate(
            get_json_object(annonars_range_response)
        )
        result = SeqVarPVS1Helper()._count_lof_variants(seqvar, 1, 1000)  # Real range is mocked
        assert result == expected_result


@pytest.mark.parametrize(
    "gene_transcripts_file,transcript_id,hgnc_id,var_pos,expected_result",
    [
        (
            "mehari/f10_mehari_genes_NM_000504.4.json",
            "NM_000504.4",
            "HGNC:3528",
            500,
            True,
        ),  # Strand plus
        (
            "mehari/f10_mehari_genes_NM_000504.4.json",
            "NM_000504.4",
            "HGNC:3528",
            1000,
            False,
        ),  # Strand plus
        (
            "mehari/pcid2_mehari_genes.json",
            "NM_001127202.4",
            "HGNC:25653",
            900,
            True,
        ),  # Strand minus. Not a frameshift!
        (
            "mehari/pcid2_mehari_genes.json",
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
    result = SeqVarPVS1Helper()._undergo_nmd(var_pos, hgnc_id, strand, exons)
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
def test_in_biologically_relevant_transcript(transcript_tags, expected_result):
    """Test the _in_biologically_relevant_transcript method."""
    result = SeqVarPVS1Helper._in_biologically_relevant_transcript(transcript_tags)
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
def test_critical4protein_function(
    seqvar, pathogenic_variants, total_variants, strand, expected_result, monkeypatch
):
    """Test the _critical4protein_function method."""
    # Mock exons, since they won't affect the outcome
    exons = [MagicMock(spec=Exon)]
    # Mocking _calculate_altered_region to return a mocked range
    mock_calculate = MagicMock(return_value=(1, 1000))
    monkeypatch.setattr(SeqVarPVS1Helper, "_calculate_altered_region", mock_calculate)
    # Mocking _count_pathogenic_variants to return controlled counts of pathogenic and total variants
    mock_count_pathogenic = MagicMock(return_value=(pathogenic_variants, total_variants))
    monkeypatch.setattr(SeqVarPVS1Helper, "_count_pathogenic_variants", mock_count_pathogenic)

    result = SeqVarPVS1Helper()._critical4protein_function(seqvar, exons, strand)  # type: ignore
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
def test_lof_is_frequent_in_population(
    seqvar, frequent_lof_variants, lof_variants, strand, expected_result, monkeypatch
):
    # Mocking exons, since they won't affect the outcome
    exons = [MagicMock(spec=Exon)]
    # Mocking _calculate_altered_region to return a mocked range
    mock_calculate = MagicMock(return_value=(1, 1000))
    monkeypatch.setattr(SeqVarPVS1Helper, "_find_affected_exon_position", mock_calculate)
    # Mocking _count_lof_variants to return controlled counts of frequent and total LoF variants
    mock_count_lof_variants = MagicMock(return_value=(frequent_lof_variants, lof_variants))
    monkeypatch.setattr(SeqVarPVS1Helper, "_count_lof_variants", mock_count_lof_variants)

    result = SeqVarPVS1Helper()._lof_is_frequent_in_population(seqvar, exons, strand)  # type: ignore
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
def test_lof_removes_more_then_10_percent_of_protein(prot_pos, prot_length, expected_result):
    """Test the _lof_removes_more_then_10_percent_of_protein method."""
    result = SeqVarPVS1Helper()._lof_removes_more_then_10_percent_of_protein(prot_pos, prot_length)
    assert result == expected_result


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
    "skipping_exon_pos_output, consequences, cryptic_ss_output, expected",
    [
        (
            (90, 120),  # _skipping_exon_pos output
            ["splice_donor_variant"],  # consequences
            [(95, "some_seq", 5.0)],  # get_cryptic_ss output
            True,
        ),
        (
            (90, 123),
            ["splice_acceptor_variant"],
            [],
            False,
        ),
        (
            (90, 120),
            ["splice_donor_variant"],
            [(101, "some_seq", 5.0)],
            True,
        ),
        (
            (90, 120),
            ["splice_donor_variant"],
            [(103, "some_seq", 5.0)],
            False,
        ),
    ],
)
def test_exon_skipping_or_cryptic_ss_disruption(
    seqvar, skipping_exon_pos_output, consequences, cryptic_ss_output, expected
):
    """Test the _exon_skipping_or_cryptic_ss_disruption method."""
    exons = [MockExon(90, 120, 90, 120)]
    with patch.object(
        SeqVarPVS1Helper, "_skipping_exon_pos", return_value=skipping_exon_pos_output
    ):
        with patch.object(SplicingPrediction, "get_sequence", return_value="some_sequence"):
            with patch.object(SplicingPrediction, "get_cryptic_ss", return_value=cryptic_ss_output):
                result = SeqVarPVS1Helper()._exon_skipping_or_cryptic_ss_disruption(
                    seqvar, exons, consequences  # type: ignore
                )
                assert result == expected, f"Expected {expected}, but got {result}"


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
def test_closest_alternative_start_codon(hgvs, cds_info, expected_result):
    """Test the _closest_alternative_start_codon method."""
    result = SeqVarPVS1Helper()._closest_alternative_start_codon(cds_info, hgvs)
    assert result == expected_result


def test_closest_alternative_start_codon_invalid():
    """Test the _closest_alternative_start_codon method."""
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
        SeqVarPVS1Helper()._closest_alternative_start_codon(cds_info, hgvs)  # type: ignore


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
def test_alternative_start_codon(hgvs, cds_info, expected_result):
    """Test the _alternative_start_codon method."""
    result = SeqVarPVS1Helper()._alternative_start_codon(cds_info, hgvs)
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
def test_upstream_pathogenic_variants(
    seqvar, exons, strand, pathogenic_variants, hgvs, cds_info, expected_result, monkeypatch
):
    """Test the _upstream_pathogenic_variants method."""
    # Mocking _count_pathogenic_variants to return a controlled number of pathogenic variants
    mock_count_pathogenic = MagicMock(return_value=(pathogenic_variants, 10))
    monkeypatch.setattr(SeqVarPVS1Helper, "_count_pathogenic_variants", mock_count_pathogenic)

    result = SeqVarPVS1Helper()._upstream_pathogenic_variants(seqvar, exons, strand, cds_info, hgvs)
    assert result == expected_result


# === SeqVarPVS1 ===


def test_init(seqvar):
    """Test the initialization of SeqVarPVS1."""
    pvs1 = SeqVarPVS1(seqvar)
    assert pvs1.config is not None
    assert pvs1.annonars_client is not None
    assert pvs1.seqvar == seqvar
    assert pvs1._seqvar_transcript is None
    assert pvs1._gene_transcript is None
    assert pvs1._consequence == SeqVarConsequence.NotSet
    assert pvs1.HGVS == ""
    assert pvs1.HGNC_id == ""
    assert len(pvs1.transcript_tags) == 0
    assert len(pvs1.exons) == 0
    assert pvs1.tx_pos_utr == -1
    assert pvs1.prot_pos == -1
    assert pvs1.prot_length == -1
    assert pvs1.cds_info == {}
    assert pvs1.strand == None
    assert pvs1.prediction == PVS1Prediction.NotPVS1
    assert pvs1.prediction_path == PVS1PredictionSeqVarPath.NotSet
