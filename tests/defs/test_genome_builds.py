import pytest

from src.defs.genome_builds import GenomeRelease, MappingError, refseq_to_genome_build


def test_genome_release_from_string_valid():
    """Test that GenomeRelease.from_string returns the correct GenomeRelease."""
    assert GenomeRelease.from_string("hg19") == GenomeRelease.GRCh37
    assert GenomeRelease.from_string("hg38") == GenomeRelease.GRCh38
    assert GenomeRelease.from_string("grch37") == GenomeRelease.GRCh37
    assert GenomeRelease.from_string("grch38") == GenomeRelease.GRCh38


def test_genome_release_from_string_invalid():
    """Test that GenomeRelease.from_string returns None for invalid values."""
    assert GenomeRelease.from_string("hg20") is None


def test_genome_release_list():
    """Test that GenomeRelease.list returns the correct list of genome releases."""
    expected = ["GRCh37", "GRCh38"]
    assert GenomeRelease.list() == expected


@pytest.mark.parametrize(
    "refseq_acc, expected_build",
    [
        ("NC_000001.10", GenomeRelease.GRCh37),
        ("NC_000001.11", GenomeRelease.GRCh38),
    ],
)
def test_refseq_to_genome_build_valid(refseq_acc, expected_build):
    """Test that refseq_to_genome_build returns the correct GenomeRelease."""
    assert refseq_to_genome_build(refseq_acc) == expected_build


def test_refseq_to_genome_build_invalid():
    """Test that refseq_to_genome_build raises MappingError for invalid RefSeq accessions."""
    with pytest.raises(MappingError):
        refseq_to_genome_build("NC_UNKNOWN")
