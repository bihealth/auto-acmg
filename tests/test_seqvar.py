from unittest.mock import Mock

import pytest

from src.api.dotty import DottyClient
from src.core.exceptions import InvalidPos, ParseError
from src.genome_builds import GenomeRelease
from src.seqvar import SeqVar, SeqVarResolver


@pytest.fixture
def seqvar_resolver():
    return SeqVarResolver()


# ===== SeqVar tests =====


def test_seqvar_initialization():
    """Test SeqVar initialization."""
    variant = SeqVar(GenomeRelease.GRCh37, "chr1", 100, "A", "T")
    assert variant.genome_release == GenomeRelease.GRCh37
    assert variant.chrom == "1"
    assert variant.pos == 100
    assert variant.delete == "A"
    assert variant.insert == "T"
    assert variant.user_representation == "GRCh37-1-100-A-T"


def test_seqvar_initialization_user_representation():
    """Test SeqVar initialization with custom user representation."""
    variant = SeqVar(GenomeRelease.GRCh37, "chr1", 100, "A", "T", user_representation="1:100:A:T")
    assert variant.user_representation == "1:100:A:T"


@pytest.mark.parametrize(
    "input_chrom, expected_normalized_chrom",
    [
        ("chr1", "1"),
        ("CHRm", "MT"),
        ("chrX", "X"),
        ("23", "23"),
        ("X", "X"),
    ],
)
def test_seqvar_normalize_chromosome(input_chrom, expected_normalized_chrom):
    """Test SeqVar._normalize_chromosome method."""
    resolver = SeqVarResolver()
    assert resolver._normalize_chrom(input_chrom) == expected_normalized_chrom


# ===== SeqVarResolver tests =====


def test_sr_validate_seqvar_valid(seqvar_resolver):
    """Test that valid SeqVar is returned unchanged."""
    variant = SeqVar(GenomeRelease.GRCh37, "1", 100, "A", "T")
    assert seqvar_resolver._validate_seqvar(variant) is variant


@pytest.mark.parametrize(
    "genome_release, chrom, pos, delete, insert",
    [
        (GenomeRelease.GRCh37, "1", 0, "A", "T"),  # Position 0 is invalid
        (GenomeRelease.GRCh37, "1", -1, "C", "G"),  # Negative position is invalid
        (GenomeRelease.GRCh37, "1", 10e8, "A", "T"),  # Out of range position is invalid
    ],
)
def test_sr_validate_seqvar_invalid_pos(
    seqvar_resolver, genome_release, chrom, pos, delete, insert
):
    """Test that InvalidPos is raised for invalid position."""
    with pytest.raises(InvalidPos):
        variant = SeqVar(genome_release, chrom, pos, delete, insert)
        seqvar_resolver._validate_seqvar(variant)


@pytest.mark.parametrize(
    "input_chrom, expected_normalized_chrom",
    [
        ("chr1", "1"),
        ("CHRm", "MT"),
        ("chrX", "X"),
        ("23", "23"),
        ("X", "X"),
    ],
)
def test_sr_normalize_chrom(seqvar_resolver, input_chrom, expected_normalized_chrom):
    """Test SeqVarResolver._normalize_chrom method."""
    assert seqvar_resolver._normalize_chrom(input_chrom) == expected_normalized_chrom


# TODO: Add more test cases for gnomad and spdi representations
@pytest.mark.parametrize(
    "representation, expected",
    [
        ("GRCh38-1-100-A-T", SeqVar(GenomeRelease.GRCh38, "1", 100, "A", "T")),
        ("GRCh37-1-100-A-T", SeqVar(GenomeRelease.GRCh37, "1", 100, "A", "T")),
    ],
)
def test_sr_parse_separated_seqvar(seqvar_resolver, representation, expected):
    variant = seqvar_resolver._parse_separated_seqvar(representation)
    assert variant.genome_release == expected.genome_release
    assert variant.chrom == expected.chrom
    assert variant.pos == expected.pos
    assert variant.delete == expected.delete
    assert variant.insert == expected.insert
    assert variant.user_representation == expected.user_representation


@pytest.mark.parametrize(
    "representation",
    [
        "GRCh38:1-100-A-T",  # Mixed separator
        "1:100:A:T:T",  # Extra field
        "A-T",  # Missing fields
    ],
)
def test_sr_parse_separated_seqvar_fail(seqvar_resolver, representation):
    with pytest.raises(ParseError):
        seqvar_resolver._parse_separated_seqvar(representation)


# TODO: Add more test cases for canonical spdi representation
@pytest.mark.parametrize(
    "value, expected",
    [
        (
            "NC_000001.10:100:A:T",
            SeqVar(GenomeRelease.GRCh37, "1", 100, "A", "T", "NC_000001.10:100:A:T"),
        ),
        (
            "NC_000001.11:200:G:C",
            SeqVar(GenomeRelease.GRCh38, "1", 200, "G", "C", "NC_000001.11:200:G:C"),
        ),
    ],
)
def test_sr_parse_canonical_spdi_seqvar_success(seqvar_resolver, value, expected):
    variant = seqvar_resolver._parse_canonical_spdi_seqvar(value)
    assert variant.genome_release == expected.genome_release
    assert variant.chrom == expected.chrom
    assert variant.pos == expected.pos
    assert variant.delete == expected.delete
    assert variant.insert == expected.insert
    assert variant.user_representation == expected.user_representation


@pytest.mark.parametrize(
    "value",
    [
        "NC_999999.99:100:A:T",  # Non-existent NC sequence
        "NC_000001.10:100:A",  # Missing field
        "XYZ_000001.10:100:A:T",  # Incorrect NC sequence format
    ],
)
def test_sr_parse_canonical_spdi_seqvar_fail(seqvar_resolver, value):
    with pytest.raises(ParseError):
        seqvar_resolver._parse_canonical_spdi_seqvar(value)


def test_sr_resolve_parse_error(seqvar_resolver):
    with pytest.raises(ParseError):
        seqvar_resolver.resolve_seqvar("invalid_variant_representation", GenomeRelease.GRCh38)


# TODO: Fix the following tests, where DottyClient is mocked
# # Mocking DottyClient's response
# @pytest.fixture
# def mock_dotty_client(monkeypatch):
#     mock_client = Mock(spec=DottyClient)
#     mock_spdi = Mock()
#     mock_spdi.success = True
#     mock_spdi.value = Mock(
#         assembly="GRCh38",
#         contig="1",
#         pos=100,
#         reference_deleted="A",
#         alternate_inserted="T"
#     )
#     mock_client.to_spdi.return_value = mock_spdi
#     monkeypatch.setattr("src.sequence_variant.DottyClient", lambda: mock_client)
#     return mock_client

# @pytest.mark.parametrize(
#     "value, genome_release, expected",
#     [
#         ("GRCh38-1-100-A-T", GenomeRelease.GRCh38, SeqVar(GenomeRelease.GRCh38, "1", 100, "A", "T")),
#         ("NC_000001.11:200:G:C", GenomeRelease.GRCh38, SeqVar(GenomeRelease.GRCh38, "1", 200, "G", "C")),
#         # Add dbSNP ID that DottyClient can resolve
#         ("rs123456", GenomeRelease.GRCh38, SeqVar(GenomeRelease.GRCh38, "1", 100, "A", "T", "rs123456")),
#     ],
# )
# def test_resolve_seqvar_success(value, genome_release, expected, mock_dotty_client, seqvar_resolver):
#     variant = seqvar_resolver.resolve_seqvar(value, genome_release)
#     assert variant.genome_release == expected.genome_release
#     assert variant.chrom == expected.chrom
#     assert variant.pos == expected.pos
#     assert variant.delete == expected.delete
#     assert variant.insert == expected.insert
#     assert variant.user_representation == expected.user_representation

# @pytest.mark.parametrize(
#     "value, genome_release",
#     [
#         ("invalid_format", GenomeRelease.GRCh38),  # Unsupported format
#         ("NC_000999.1:100:A:T", GenomeRelease.GRCh38),  # Non-existent NC sequence
#         # Add any other cases that should fail
#     ],
# )
# def test_resolve_seqvar_failure(value, genome_release, seqvar_resolver):
#     with pytest.raises(ParseError):
#         seqvar_resolver.resolve_seqvar(value, genome_release)
