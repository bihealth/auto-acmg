from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGSeqVarData, AutoACMGSeqVarResult
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh38, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def auto_acmg_result():
    return AutoACMGSeqVarResult()


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


@pytest.fixture
def default_predictor(seqvar, auto_acmg_result):
    return DefaultSeqVarPredictor(seqvar=seqvar, result=auto_acmg_result, config=MagicMock())


@patch("src.seqvar.default_predictor.DefaultSeqVarPredictor.predict_pvs1")
@patch("src.seqvar.default_predictor.DefaultSeqVarPredictor.predict_ps1pm5")
@patch("src.seqvar.default_predictor.DefaultSeqVarPredictor.predict_pm1")
@patch("src.seqvar.default_predictor.DefaultSeqVarPredictor.predict_pm2ba1bs1bs2")
@patch("src.seqvar.default_predictor.DefaultSeqVarPredictor.predict_pm4bp3")
@patch("src.seqvar.default_predictor.DefaultSeqVarPredictor.predict_pp2bp1")
@patch("src.seqvar.default_predictor.DefaultSeqVarPredictor.predict_pp3bp4")
@patch("src.seqvar.default_predictor.DefaultSeqVarPredictor.predict_bp7")
def test_predict(
    mock_predict_bp7,
    mock_predict_pp3bp4,
    mock_predict_pp2bp1,
    mock_predict_pm4bp3,
    mock_predict_pm2ba1bs1bs2,
    mock_predict_pm1,
    mock_predict_ps1pm5,
    mock_predict_pvs1,
    default_predictor,
    auto_acmg_result,
    auto_acmg_data,
):
    """Test the full prediction flow of DefaultSeqVarPredictor."""

    # Setting return values for each criterion prediction
    mock_predict_pvs1.return_value = MagicMock()
    mock_predict_ps1pm5.return_value = (MagicMock(), MagicMock())
    mock_predict_pm1.return_value = MagicMock()
    mock_predict_pm2ba1bs1bs2.return_value = (
        MagicMock(),
        MagicMock(),
        MagicMock(),
        MagicMock(),
    )
    mock_predict_pm4bp3.return_value = (MagicMock(), MagicMock())
    mock_predict_pp2bp1.return_value = (MagicMock(), MagicMock())
    mock_predict_pp3bp4.return_value = (MagicMock(), MagicMock())
    mock_predict_bp7.return_value = MagicMock()

    result = default_predictor.predict()

    assert result == auto_acmg_result, "Predict should return the AutoACMGResult object."

    # Asserting that each prediction method was called
    mock_predict_pvs1.assert_called_once_with(
        default_predictor.seqvar, default_predictor.result.data
    )
    mock_predict_ps1pm5.assert_called_once_with(
        default_predictor.seqvar, default_predictor.result.data
    )
    mock_predict_pm1.assert_called_once_with(
        default_predictor.seqvar, default_predictor.result.data
    )
    mock_predict_pm2ba1bs1bs2.assert_called_once_with(
        default_predictor.seqvar, default_predictor.result.data
    )
    mock_predict_pm4bp3.assert_called_once_with(
        default_predictor.seqvar, default_predictor.result.data
    )
    mock_predict_pp2bp1.assert_called_once_with(
        default_predictor.seqvar, default_predictor.result.data
    )
    mock_predict_pp3bp4.assert_called_once_with(
        default_predictor.seqvar, default_predictor.result.data
    )
    mock_predict_bp7.assert_called_once_with(
        default_predictor.seqvar, default_predictor.result.data
    )

    # Asserting that the predictions were stored in the result
    assert result.criteria.pvs1 is not None, "PVS1 prediction should be stored in the result."
    assert result.criteria.ps1 is not None, "PS1 prediction should be stored in the result."
    assert result.criteria.pm1 is not None, "PM1 prediction should be stored in the result."
    assert result.criteria.pm2 is not None, "PM2 prediction should be stored in the result."
    assert result.criteria.bp3 is not None, "BP3 prediction should be stored in the result."
    assert result.criteria.bp7 is not None, "BP7 prediction should be stored in the result."
