"""Integration tests for `PM1` criteria using upstream server."""

from typing import List, Tuple

import pytest

from src.auto_acmg import VCEP_MAPPING, AutoACMG
from src.core.config import Config
from src.defs.auto_acmg import AutoACMGPrediction, AutoACMGSeqVarResult
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.auto_pm1 import AutoPM1
from tests.utils import load_test_data_pm1

#: Type for PM1 test data.
Pm1TestData = List[Tuple[str, GenomeRelease, AutoACMGPrediction]]
#: Test data.
PM1_TEST_DATA: Pm1TestData = load_test_data_pm1("tests/assets/integ/pm1.csv")


@pytest.mark.default_cassette("integ_pm1.yaml")
@pytest.mark.vcr
@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction",
    PM1_TEST_DATA,
)
def test_pm1(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: bool,
    config: Config,
):
    # First, resolve variant
    auto_acmg = AutoACMG(variant_name, genome_release, config=config)
    seqvar = auto_acmg.resolve_variant()
    assert isinstance(seqvar, SeqVar)
    # Then, setup the data
    auto_acmg_result = auto_acmg.parse_seqvar_data(seqvar)
    assert isinstance(auto_acmg_result, AutoACMGSeqVarResult)
    # Then, predict PM1
    if auto_acmg_result.data.hgnc_id in VCEP_MAPPING:
        predictor_class = VCEP_MAPPING[auto_acmg_result.data.hgnc_id]
        predictor = predictor_class(seqvar, auto_acmg_result, config)
        pm1 = predictor.predict_pm1(seqvar, auto_acmg_result.data)
    else:
        auto_pm1 = AutoPM1()
        pm1 = auto_pm1.predict_pm1(seqvar, auto_acmg_result.data)
    assert pm1.name == "PM1"
    assert pm1.prediction == expected_prediction
