"""End-to-end system tests for `auto_acmg` using upstream server."""

from typing import List, Tuple

import pytest

from src.auto_acmg import AutoACMG, AutoACMGCriteria
from src.core.config import Config
from src.defs.auto_acmg import ACMGCriteria
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from tests.utils import load_test_data

#: Type for ACMG criteria test data.
AcmgTestData = List[Tuple[str, GenomeRelease, ACMGCriteria, str]]
#: Test data.
ACMG_TEST_DATA: AcmgTestData = load_test_data("tests/assets/e2e_variants/other_criteria.csv")


@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction, comment",
    ACMG_TEST_DATA,
)
def acmg_criteria_test_helper(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: ACMGCriteria,
    comment: str,
    config: Config,
):
    # First, resolve variant
    auto_acmg = AutoACMG(variant_name, genome_release, config=config)
    seqvar = auto_acmg.resolve_variant()
    assert isinstance(seqvar, SeqVar)
    # Then, predict ACMG criteria
    auto_acmg_criteria = AutoACMGCriteria(seqvar, genome_release, config=config)
    result = auto_acmg_criteria.predict()
    assert result == expected_prediction


@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction, comment",
    ACMG_TEST_DATA,
)
def test_acmg_criteria_seqvar_csv(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: ACMGCriteria,
    comment: str,
    config: Config,
):
    """Test ACMG criteria predictions, variants read from CSV file."""
    acmg_criteria_test_helper(variant_name, genome_release, expected_prediction, comment, config)


@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction, comment",
    [
        (
            "3-178936091-C-T",
            GenomeRelease.GRCh37,
            ACMGCriteria(),
            "Hotspot with 4 pathogenic variants",
        ),
        (
            "7-140453136-C-T",
            GenomeRelease.GRCh37,
            ACMGCriteria(),
            "UniProt domain with 2 pathogenic variants",
        ),
        (
            "1-55516888-G-A",
            GenomeRelease.GRCh37,
            ACMGCriteria(),
            "No hotspot or critical domain",
        ),
    ],
)
def test_acmg_criteria_seqvar_inline(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: ACMGCriteria,
    comment: str,
    config: Config,
):
    """Test ACMG criteria predictions, variants defined inline."""
    acmg_criteria_test_helper(variant_name, genome_release, expected_prediction, comment, config)
