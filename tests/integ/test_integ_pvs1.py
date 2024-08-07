"""Integration tests for `PVS1` criteria using upstream server."""

from typing import List, Tuple

import pytest

from src.auto_acmg import AutoACMG
from src.core.config import Config
from src.criteria.auto_pvs1 import AutoPVS1
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGResult
from src.defs.auto_pvs1 import PVS1Prediction, PVS1PredictionSeqVarPath
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from tests.utils import load_test_data_pvs1

#: Type for PVS1 test data.
Pvs1TestData = List[Tuple[str, GenomeRelease, PVS1Prediction, PVS1PredictionSeqVarPath]]
#: Test data.
PVS1_TEST_DATA: Pvs1TestData = load_test_data_pvs1("tests/assets/integ/pvs1.csv")


@pytest.mark.default_cassette("integ_pvs1.yaml")
@pytest.mark.vcr
@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction, expected_path",
    PVS1_TEST_DATA,
)
def test_pvs1(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: PVS1Prediction,
    expected_path: PVS1PredictionSeqVarPath,
    config: Config,
):
    # First, resolve variant
    auto_acmg = AutoACMG(variant_name, genome_release, config=config)
    seqvar = auto_acmg.resolve_variant()
    assert isinstance(seqvar, SeqVar)
    # Then, setup the data
    auto_acmg_result = auto_acmg.parse_data(seqvar)
    assert isinstance(auto_acmg_result, AutoACMGResult)
    # Then, predict PVS1
    auto_pvs1 = AutoPVS1()
    pvs1 = auto_pvs1.predict_pvs1(seqvar, auto_acmg_result.data)
    assert pvs1.name == "PVS1"
