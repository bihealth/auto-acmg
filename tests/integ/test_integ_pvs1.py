"""Integration tests for `PVS1` criteria using upstream server."""

from typing import Dict, List, Tuple

import pytest

from src.auto_acmg import VCEP_MAPPING, AutoACMG
from src.core.config import Config
from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarResult,
    AutoACMGStrength,
)
from src.defs.auto_pvs1 import PVS1Prediction, PVS1PredictionSeqVarPath
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.auto_pvs1 import AutoPVS1
from tests.utils import load_test_data_pvs1

#: Type for PVS1 test data.
Pvs1TestData = List[Tuple[str, GenomeRelease, PVS1Prediction, PVS1PredictionSeqVarPath]]
#: Test data.
PVS1_TEST_DATA: Pvs1TestData = load_test_data_pvs1("tests/assets/integ/pvs1.csv")

evidence_strength_mapping: Dict[PVS1Prediction, AutoACMGStrength] = {
    PVS1Prediction.PVS1: AutoACMGStrength.PathogenicVeryStrong,
    PVS1Prediction.PVS1_Strong: AutoACMGStrength.PathogenicStrong,
    PVS1Prediction.PVS1_Moderate: AutoACMGStrength.PathogenicModerate,
    PVS1Prediction.PVS1_Supporting: AutoACMGStrength.PathogenicSupporting,
    PVS1Prediction.NotPVS1: AutoACMGStrength.PathogenicVeryStrong,
    PVS1Prediction.UnsupportedConsequence: AutoACMGStrength.PathogenicVeryStrong,
    PVS1Prediction.NotSet: AutoACMGStrength.PathogenicVeryStrong,
}

met_pred: List[PVS1Prediction] = [
    PVS1Prediction.PVS1,
    PVS1Prediction.PVS1_Strong,
    PVS1Prediction.PVS1_Moderate,
    PVS1Prediction.PVS1_Supporting,
]


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
    auto_acmg_result = auto_acmg.parse_seqvar_data(seqvar)
    assert isinstance(auto_acmg_result, AutoACMGSeqVarResult)
    # Then, predict PVS1
    if auto_acmg_result.data.hgnc_id in VCEP_MAPPING:
        predictor_class = VCEP_MAPPING[auto_acmg_result.data.hgnc_id]
        predictor = predictor_class(seqvar, auto_acmg_result, config)
        pvs1 = predictor.predict_pvs1(seqvar, auto_acmg_result.data)
    else:
        auto_pvs1 = AutoPVS1()
        pvs1 = auto_pvs1.predict_pvs1(seqvar, auto_acmg_result.data)
    expected_pred = (
        AutoACMGPrediction.Met if expected_prediction in met_pred else AutoACMGPrediction.NotMet
    )
    assert pvs1.name == "PVS1"
    assert pvs1.prediction == expected_pred
    assert pvs1.strength == evidence_strength_mapping[expected_prediction]
