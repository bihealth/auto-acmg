from unittest.mock import Mock

import pytest
from typer.testing import CliRunner

from src.defs.autopvs1 import PVS1PredictionStrucVarPath
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar, SeqVarResolver
from src.defs.strucvar import StrucVar, StrucVarResolver, StrucVarType
from src.pvs1.auto_pvs1 import AutoPVS1
from src.pvs1.seqvar_pvs1 import PVS1Prediction, PVS1PredictionSeqVarPath, SeqVarPVS1
from src.pvs1.strucvar_pvs1 import StrucVarPVS1

runner = CliRunner()


@pytest.fixture
def mock_seqvar():
    return SeqVar(
        genome_release=GenomeRelease.GRCh38, chrom="1", pos=100000, delete="A", insert="T"
    )


@pytest.fixture
def mock_strucvar():
    return StrucVar(
        genome_release=GenomeRelease.GRCh38,
        chrom="1",
        start=100000,
        stop=200000,
        sv_type=StrucVarType.DEL,
    )


@pytest.fixture
def mock_seqvar_pvs1(monkeypatch):
    mock_pvs1 = Mock(SeqVarPVS1)
    mock_pvs1.initialize.return_value = None
    mock_pvs1.verify_PVS1.return_value = None
    mock_pvs1.get_prediction.return_value = (PVS1Prediction.PVS1, PVS1PredictionSeqVarPath.NF1)
    monkeypatch.setattr("src.pvs1.auto_pvs1.SeqVarPVS1", lambda *args, **kwargs: mock_pvs1)
    return mock_pvs1


@pytest.fixture
def mock_seqvar_pvs1_failure(monkeypatch):
    mock_pvs1 = Mock(SeqVarPVS1)
    mock_pvs1.initialize.return_value = None
    mock_pvs1.verify_PVS1.return_value = None
    mock_pvs1.get_prediction.return_value = (None, None)
    monkeypatch.setattr("src.pvs1.auto_pvs1.SeqVarPVS1", lambda *args, **kwargs: mock_pvs1)
    return mock_pvs1


@pytest.fixture
def mock_strucvar_pvs1(monkeypatch):
    mock_pvs1 = Mock(StrucVarPVS1)
    mock_pvs1.initialize.return_value = None
    mock_pvs1.verify_PVS1.return_value = None
    mock_pvs1.get_prediction.return_value = (PVS1Prediction.PVS1, PVS1PredictionStrucVarPath.DEL1)
    monkeypatch.setattr("src.pvs1.auto_pvs1.StrucVarPVS1", lambda *args, **kwargs: mock_pvs1)
    return mock_pvs1


@pytest.fixture
def mock_strucvar_pvs1_failure(monkeypatch):
    mock_pvs1 = Mock(StrucVarPVS1)
    mock_pvs1.initialize.return_value = None
    mock_pvs1.verify_PVS1.return_value = None
    mock_pvs1.get_prediction.return_value = (None, None)
    monkeypatch.setattr("src.pvs1.auto_pvs1.StrucVarPVS1", lambda *args, **kwargs: mock_pvs1)
    return mock_pvs1


def test_autoPVS1_predict_seqvar_success(mock_seqvar_pvs1, mock_seqvar):
    """Test predict method with a successful response."""
    autoPVS1 = AutoPVS1(mock_seqvar)
    with runner.isolated_filesystem():
        autoPVS1.predict()
    assert mock_seqvar_pvs1.initialize.called
    assert mock_seqvar_pvs1.verify_PVS1.called
    assert mock_seqvar_pvs1.get_prediction.called
    assert autoPVS1.seqvar_prediction == PVS1Prediction.PVS1
    assert autoPVS1.seqvar_prediction_path == PVS1PredictionSeqVarPath.NF1


def test_autoPVS1_predict_seqvar_failure(mock_seqvar_pvs1_failure, mock_seqvar):
    """Test predict method with a failed response."""
    autoPVS1 = AutoPVS1(mock_seqvar)
    with runner.isolated_filesystem():
        autoPVS1.predict()
    assert mock_seqvar_pvs1_failure.initialize.called
    assert mock_seqvar_pvs1_failure.verify_PVS1.called
    assert mock_seqvar_pvs1_failure.get_prediction.called
    assert autoPVS1.seqvar_prediction == None
    assert autoPVS1.seqvar_prediction_path == None


def test_autoPVS1_predict_strucvar_success(mock_strucvar_pvs1, mock_strucvar):
    """Test predict method with a successful response."""
    autoPVS1 = AutoPVS1(mock_strucvar)
    with runner.isolated_filesystem():
        autoPVS1.predict()
    assert mock_strucvar_pvs1.initialize.called
    assert mock_strucvar_pvs1.verify_PVS1.called
    assert mock_strucvar_pvs1.get_prediction.called
    assert autoPVS1.strucvar_prediction == PVS1Prediction.PVS1
    assert autoPVS1.strucvar_prediction_path == PVS1PredictionStrucVarPath.DEL1


def test_autoPVS1_predict_strucvar_failure(mock_strucvar_pvs1_failure, mock_strucvar):
    """Test predict method with a failed response."""
    autoPVS1 = AutoPVS1(mock_strucvar)
    with runner.isolated_filesystem():
        autoPVS1.predict()
    assert mock_strucvar_pvs1_failure.initialize.called
    assert mock_strucvar_pvs1_failure.verify_PVS1.called
    assert mock_strucvar_pvs1_failure.get_prediction.called
    assert autoPVS1.strucvar_prediction == None
    assert autoPVS1.strucvar_prediction_path == None
