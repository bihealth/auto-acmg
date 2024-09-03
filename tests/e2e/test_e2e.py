"""End-to-end system tests for `auto_acmg` using upstream server."""

from typing import List, Tuple

import pytest

from src.auto_acmg import AutoACMG
from src.core.config import Config
from src.defs.auto_acmg import AutoACMGSeqVarResult
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from tests.utils import load_test_data

#: Type for ACMG criteria test data.
AcmgTestData = List[Tuple[str, GenomeRelease, AutoACMGSeqVarResult, str]]
#: Test data.
ACMG_TEST_DATA: AcmgTestData = load_test_data("tests/assets/e2e/variants.csv")


@pytest.mark.default_cassette("e2e.yaml")
@pytest.mark.vcr
@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction, comment",
    ACMG_TEST_DATA,
)
def acmg_criteria_test_helper(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: List[str],
    comment: str,
    config: Config,
):
    # First, resolve variant
    auto_acmg = AutoACMG(variant_name, genome_release, config=config)
    seqvar = auto_acmg.resolve_variant()
    assert isinstance(seqvar, SeqVar)
    # Then, predict ACMG criteria
    prediction = auto_acmg.predict()
    assert isinstance(prediction, AutoACMGSeqVarResult)


@pytest.mark.default_cassette("e2e.yaml")
@pytest.mark.vcr
@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction, comment",
    ACMG_TEST_DATA,
)
def test_acmg_criteria_seqvar_csv(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: List[str],
    comment: str,
    config: Config,
):
    """Test ACMG criteria predictions, variants read from CSV file."""
    acmg_criteria_test_helper(variant_name, genome_release, expected_prediction, comment, config)


@pytest.mark.default_cassette("e2e.yaml")
@pytest.mark.vcr
@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction, comment",
    [
        # (
        #     "3-178936091-C-T",
        #     GenomeRelease.GRCh37,
        #     ACMGCriteria(),
        #     "Hotspot with 4 pathogenic variants",
        # ),
    ],
)
def test_acmg_criteria_seqvar_inline(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: str,
    comment: str,
    config: Config,
):
    """Test ACMG criteria predictions, variants defined inline."""
    acmg_criteria_test_helper(variant_name, genome_release, expected_prediction, comment, config)
