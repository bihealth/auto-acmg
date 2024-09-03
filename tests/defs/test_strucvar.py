import pytest

from src.defs.exceptions import InvalidPos, ParseError
from src.defs.genome_builds import GenomeRelease
from src.defs.strucvar import StrucVar, StrucVarResolver, StrucVarType


@pytest.fixture
def strucvar_resolver():
    return StrucVarResolver(config=None)


# ===== StrucVar tests =====


def test_strucvar_initialization():
    """Test StrucVar initialization."""
    variant = StrucVar(StrucVarType.DEL, GenomeRelease.GRCh37, "chr1", 100, 200)
    assert variant.sv_type == StrucVarType.DEL
    assert variant.genome_release == GenomeRelease.GRCh37
    assert variant.chrom == "1"
    assert variant.start == 100
    assert variant.stop == 200


def test_strucvar_initialization_user_representation():
    """Test StrucVar initialization with custom user representation."""
    variant = StrucVar(
        StrucVarType.DEL, GenomeRelease.GRCh37, "chr1", 100, 200, user_repr="1:100-200"
    )


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
def test_strucvar_normalize_chromosome(input_chrom, expected_normalized_chrom):
    """Test StrucVar._normalize_chromosome method."""
    resolver = StrucVarResolver()
    assert resolver._normalize_chromosome(input_chrom) == expected_normalized_chrom


def test_strucvar_eq():
    """Test StrucVar.__eq__ method."""
    variant1 = StrucVar(StrucVarType.DEL, GenomeRelease.GRCh37, "1", 100, 200)
    variant2 = StrucVar(StrucVarType.DEL, GenomeRelease.GRCh37, "1", 100, 200)
    assert variant1 == variant2


# ===== StrucVarResolver tests =====


def test_validate_strucvar_valid(strucvar_resolver):
    """Test that valid StrucVar is returned unchanged."""
    variant = StrucVar(StrucVarType.DUP, GenomeRelease.GRCh38, "1", 100, 500)
    validated_variant = strucvar_resolver._validate_strucvar(variant)
    assert validated_variant is variant


@pytest.mark.parametrize(
    "sv_type, genome_release, chrom, start, stop",
    [
        (StrucVarType.DEL, GenomeRelease.GRCh37, "1", 0, 100),  # Position 0 is invalid
        (StrucVarType.DEL, GenomeRelease.GRCh37, "1", -1, 100),  # Negative position is invalid
        (
            StrucVarType.DEL,
            GenomeRelease.GRCh37,
            "1",
            10e8,
            100,
        ),  # Out of range position is invalid
    ],
)
def test_validate_strucvar_invalid(strucvar_resolver, sv_type, genome_release, chrom, start, stop):
    """Test that invalid StrucVar raises an exception."""
    with pytest.raises(InvalidPos):
        variant = StrucVar(sv_type, genome_release, chrom, start, stop)
        strucvar_resolver._validate_strucvar(variant)


@pytest.mark.parametrize(
    "input_str,expected_output",
    [
        ("DEL-GRCh38-1-100-200", StrucVar(StrucVarType.DEL, GenomeRelease.GRCh38, "1", 100, 200)),
        (
            "DUP:GRCh37:chrX:1000:1500",
            StrucVar(StrucVarType.DUP, GenomeRelease.GRCh37, "X", 1000, 1500),
        ),
    ],
)
def test_parse_separated_strucvar(strucvar_resolver, input_str, expected_output):
    """Test the _parse_separated_strucvar method."""
    result = strucvar_resolver._parse_separated_strucvar(input_str)
    assert result.sv_type == expected_output.sv_type
    assert result.genome_release == expected_output.genome_release
    assert result.chrom == expected_output.chrom
    assert result.start == expected_output.start
    assert result.stop == expected_output.stop


# TODO: Add more test cases for canonical spdi representation
@pytest.mark.parametrize(
    "input_str",
    [
        "DEL-GRCh38-1-100",  # Missing one positional value
        "DCT-GRCh38-1-200-300",  # Invalid SV type
        "DUP-GRCh38-1-100-200-300",  # Extra positional value
    ],
)
def test_parse_separated_strucvar_failure(strucvar_resolver, input_str):
    """Test that invalid input raises an exception."""
    with pytest.raises(ParseError):
        strucvar_resolver._parse_separated_strucvar(input_str)


@pytest.mark.parametrize(
    "input_str,expected_output",
    [
        ("DEL:GRCh38:1:100:200", StrucVar(StrucVarType.DEL, GenomeRelease.GRCh38, "1", 100, 200)),
        (
            "DUP:GRCh37:chrX:1000:1500",
            StrucVar(StrucVarType.DUP, GenomeRelease.GRCh37, "X", 1000, 1500),
        ),
    ],
)
def test_resolve_strucvar_success(strucvar_resolver, input_str, expected_output):
    """Test the resolve_strucvar method."""
    result = strucvar_resolver.resolve_strucvar(input_str, GenomeRelease.GRCh38)
    assert result.sv_type == expected_output.sv_type
    assert result.genome_release == expected_output.genome_release
    assert result.start == expected_output.start
    assert result.stop == expected_output.stop


@pytest.mark.parametrize(
    "input_str",
    [
        "DEL-GRCh38-1-100",  # Missing one positional value
        "DCT-GRCh38-1-200-300",  # Invalid SV type
        "DUP-GRCh38-1-100-200-300",  # Extra positional value
    ],
)
def test_resolve_strucvar_failure(strucvar_resolver, input_str):
    """Test that invalid input raises an exception."""
    with pytest.raises(ParseError):
        strucvar_resolver.resolve_strucvar(input_str, GenomeRelease.GRCh38)
