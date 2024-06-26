"""End-to-end system tests for ``auto_acmg`` using upstream server."""

import csv
from typing import Any, List, Tuple

import pytest

from src.auto_acmg import AutoACMG
from src.core.config import Config
from src.defs.auto_pvs1 import PVS1Prediction, PVS1PredictionSeqVarPath
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.pvs1.auto_pvs1 import AutoPVS1
from tests.utils import load_test_data_pvs1

#: Type for PVS1 test data.
Pvs1TestData = List[Tuple[str, GenomeRelease, PVS1Prediction, PVS1PredictionSeqVarPath]]
#: Test data.
PVS1_TEST_DATA: Pvs1TestData = load_test_data_pvs1("tests/assets/e2e_variants/pvs1.csv")


@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction, expected_path",
    PVS1_TEST_DATA,
)
def pvs1_seqvar_test_helper(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: PVS1Prediction,
    expected_path: PVS1PredictionSeqVarPath,
    config: Config,
):
    # first, resolve variant
    auto_acmg = AutoACMG(variant_name, genome_release, config=config)
    variant = auto_acmg.resolve_variant()
    assert isinstance(variant, SeqVar)
    # then, predict PVS1
    pvs1 = AutoPVS1(variant, config=config)
    prediction, path, _ = pvs1.predict()
    assert prediction == expected_prediction
    assert path == expected_path


@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction, expected_path",
    PVS1_TEST_DATA,
)
def test_pvs1_seqvar_csv(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: PVS1Prediction,
    expected_path: PVS1PredictionSeqVarPath,
    config: Config,
):
    """Test PVS1 predictions, variants read from CSV file."""
    pvs1_seqvar_test_helper(
        variant_name, genome_release, expected_prediction, expected_path, config
    )


@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction, expected_path",
    [
        (
            "4-113568536-G-GA",
            GenomeRelease.GRCh37,
            PVS1Prediction.PVS1,
            PVS1PredictionSeqVarPath.NF1,
        )
    ],
)
def test_pvs1_seqvar_inline(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: PVS1Prediction,
    expected_path: PVS1PredictionSeqVarPath,
    config: Config,
):
    """Test PVS1 predictions, variants defined inline."""
    pvs1_seqvar_test_helper(
        variant_name, genome_release, expected_prediction, expected_path, config
    )
