from unittest.mock import patch

import pytest

from src.core.config import Config
from src.defs.auto_acmg import AutoACMGStrucVarResult
from src.defs.genome_builds import GenomeRelease
from src.defs.strucvar import StrucVar, StrucVarType
from src.strucvar.default_predictor import DefaultStrucVarPredictor


@pytest.fixture
def strucvar():
    return StrucVar(
        sv_type=StrucVarType.DEL,
        genome_release=GenomeRelease.GRCh38,
        chrom="1",
        start=100,
        stop=200,
    )


@pytest.fixture
def auto_acmg_result():
    return AutoACMGStrucVarResult()


@pytest.fixture
def default_predictor(strucvar, auto_acmg_result):
    config = Config(api_base_url="https://api.example.com")
    return DefaultStrucVarPredictor(strucvar, auto_acmg_result, config)


@patch.object(DefaultStrucVarPredictor, "predict_pvs1", return_value="PVS1 Not Applicable")
def test_predict(mock_predict_pvs1, default_predictor):
    """Test the predict method to ensure it calls predict_pvs1 and logs the result correctly."""
    result = default_predictor.predict()
    mock_predict_pvs1.assert_called_once_with(
        default_predictor.strucvar, default_predictor.result.data
    )
    assert (
        result.criteria.pvs1 == "PVS1 Not Applicable"
    ), "PVS1 prediction should be captured in the result."
    assert (
        result is default_predictor.result
    ), "Predict should return the AutoACMGStrucVarResult object."


def test_predict_no_pvs1_implemented(default_predictor):
    """Test predict when PVS1 is not implemented to check if it handles absence of prediction."""
    # Assuming predict_pvs1 returns None or does not make changes when not implemented
    with patch.object(
        DefaultStrucVarPredictor, "predict_pvs1", return_value=None
    ) as mock_predict_pvs1:
        result = default_predictor.predict()
        mock_predict_pvs1.assert_called_once()
        assert "pvs1" not in result.criteria, "PVS1 should not be set if no implementation exists."
