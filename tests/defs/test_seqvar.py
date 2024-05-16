from unittest.mock import MagicMock, patch

import pytest

from src.api.dotty import DottyClient
from src.defs.dotty import DottySpdiResponse
from src.defs.exceptions import InvalidPos, ParseError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar, SeqVarResolver
from tests.utils import get_json_object


@pytest.fixture
def seqvar_resolver():
    return SeqVarResolver(config=None)


@pytest.fixture
def dotty_response_success():
    return DottySpdiResponse.model_validate(get_json_object("dotty/dotty_spdi_success.json"))


# ===== SeqVar tests =====


def test_seqvar_initialization():
    """Test SeqVar initialization."""
    variant = SeqVar(GenomeRelease.GRCh37, "chr1", 100, "A", "T")
    assert variant.genome_release == GenomeRelease.GRCh37
    assert variant.chrom == "1"
    assert variant.pos == 100
    assert variant.delete == "A"
    assert variant.insert == "T"
    assert variant.user_repr == "GRCh37-1-100-A-T"


def test_seqvar_initialization_user_representation():
    """Test SeqVar initialization with custom user representation."""
    variant = SeqVar(GenomeRelease.GRCh37, "chr1", 100, "A", "T", user_repr="1:100:A:T")
    assert variant.user_repr == "1:100:A:T"


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
def test_seqvar_normalize_chromosome(input_chrom, expected_normalized_chrom, seqvar_resolver):
    """Test SeqVar._normalize_chromosome method."""
    assert seqvar_resolver._normalize_chrom(input_chrom) == expected_normalized_chrom


# ===== SeqVarResolver tests =====


def test_validate_seqvar_valid(seqvar_resolver):
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
def test_validate_seqvar_invalid_pos(seqvar_resolver, genome_release, chrom, pos, delete, insert):
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
def test_normalize_chrom(seqvar_resolver, input_chrom, expected_normalized_chrom):
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
def test_parse_separated_seqvar(seqvar_resolver, representation, expected):
    variant = seqvar_resolver._parse_separated_seqvar(representation)
    assert variant.genome_release == expected.genome_release
    assert variant.chrom == expected.chrom
    assert variant.pos == expected.pos
    assert variant.delete == expected.delete
    assert variant.insert == expected.insert
    assert variant.user_repr == expected.user_repr


@pytest.mark.parametrize(
    "representation",
    [
        "GRCh38:1-100-A-T",  # Mixed separator
        "1:100:A:T:T",  # Extra field
        "A-T",  # Missing fields
    ],
)
def test_parse_separated_seqvar_fail(seqvar_resolver, representation):
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
def test_parse_canonical_spdi_seqvar_success(seqvar_resolver, value, expected):
    variant = seqvar_resolver._parse_canonical_spdi_seqvar(value)
    assert variant.genome_release == expected.genome_release
    assert variant.chrom == expected.chrom
    assert variant.pos == expected.pos
    assert variant.delete == expected.delete
    assert variant.insert == expected.insert
    assert variant.user_repr == expected.user_repr


@pytest.mark.parametrize(
    "value",
    [
        "NC_999999.99:100:A:T",  # Non-existent NC sequence
        "NC_000001.10:100:A",  # Missing field
        "XYZ_000001.10:100:A:T",  # Incorrect NC sequence format
    ],
)
def test_parse_canonical_spdi_seqvar_fail(seqvar_resolver, value):
    with pytest.raises(ParseError):
        seqvar_resolver._parse_canonical_spdi_seqvar(value)


@patch.object(DottyClient, "to_spdi")
def test_resolve_seqvar_success(mock_to_spdi, seqvar_resolver, dotty_response_success):
    mock_to_spdi.return_value = dotty_response_success
    variant = seqvar_resolver.resolve_seqvar("Example.3:c.1085delT", GenomeRelease.GRCh38)
    assert variant.genome_release == GenomeRelease.GRCh37
    assert variant.chrom == "4"
    assert variant.pos == 113568536
    assert variant.delete == "G"
    assert variant.insert == "GA"


@patch.object(DottyClient, "to_spdi")
def test_resolve_seqvar_failure(mock_to_spdi, seqvar_resolver):
    mock_to_spdi.side_effect = MagicMock(success=False)
    with pytest.raises(ParseError):
        seqvar_resolver.resolve_seqvar("Example.3:c.1085delT", GenomeRelease.GRCh38)
