from unittest.mock import MagicMock, Mock

import pytest
from typer.testing import CliRunner

from src.auto_acmg import AutoACMG
from src.auto_ps1_pm5 import AutoPS1PM5
from src.defs.auto_acmg import PS1PM5
from src.defs.exceptions import AutoAcmgBaseException, ParseError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar, SeqVarResolver
from src.pvs1.auto_pvs1 import AutoPVS1
from src.pvs1.seqvar_pvs1 import PVS1Prediction, PVS1PredictionSeqVarPath

runner = CliRunner()


@pytest.fixture
def mock_seqvar():
    return SeqVar(
        genome_release=GenomeRelease.GRCh38, chrom="1", pos=100000, delete="A", insert="T"
    )


@pytest.fixture
def mock_seqvar_resolver(monkeypatch):
    mock_resolver = MagicMock(SeqVarResolver)
    mock_resolver.resolve_seqvar.return_value = SeqVar(
        genome_release=GenomeRelease.GRCh38, chrom="1", pos=100000, delete="A", insert="T"
    )
    monkeypatch.setattr("src.auto_acmg.SeqVarResolver", lambda *args, **kwargs: mock_resolver)
    return mock_resolver


@pytest.fixture
def mock_seqvar_resolver_parse_error(monkeypatch):
    mock_resolver = MagicMock(SeqVarResolver)
    mock_resolver.resolve_seqvar.side_effect = ParseError("Invalid variant format")
    monkeypatch.setattr("src.auto_acmg.SeqVarResolver", lambda *args, **kwargs: mock_resolver)
    return mock_resolver


@pytest.fixture
def mock_seqvar_resolver_failure(monkeypatch):
    mock_resolver = MagicMock(SeqVarResolver)
    mock_resolver.resolve_seqvar.return_value = None
    monkeypatch.setattr("src.auto_acmg.SeqVarResolver", lambda *args, **kwargs: mock_resolver)
    return mock_resolver


@pytest.fixture
def mock_auto_pvs1_success(monkeypatch):
    mock_pvs1 = MagicMock(AutoPVS1)
    mock_pvs1.predict.return_value = (PVS1Prediction.PVS1, PVS1PredictionSeqVarPath.NF1)
    monkeypatch.setattr("src.auto_acmg.AutoPVS1", lambda *args, **kwargs: mock_pvs1)
    return mock_pvs1


@pytest.fixture
def mock_auto_pvs1_failure(monkeypatch):
    mock_pvs1 = MagicMock(AutoPVS1)
    mock_pvs1.predict.side_effect = AutoAcmgBaseException("An error occurred")
    monkeypatch.setattr("src.auto_acmg.AutoPVS1", lambda *args, **kwargs: mock_pvs1)
    return mock_pvs1


@pytest.fixture
def mock_auto_ps1_pm5_success(monkeypatch):
    mock_ps1_pm5 = MagicMock(AutoPS1PM5)
    mock_ps1_pm5.predict.return_value = PS1PM5()
    monkeypatch.setattr("src.auto_acmg.AutoPS1PM5", lambda *args, **kwargs: mock_ps1_pm5)
    return mock_ps1_pm5


@pytest.fixture
def mock_auto_ps1_pm5_failure(monkeypatch):
    mock_ps1_pm5 = MagicMock(AutoPS1PM5)
    mock_ps1_pm5.predict.side_effect = AutoAcmgBaseException("An error occurred")
    monkeypatch.setattr("src.auto_acmg.AutoPS1PM5", lambda *args, **kwargs: mock_ps1_pm5)
    return mock_ps1_pm5


def test_auto_acmg_resolve_sequence_variant_success(mock_seqvar_resolver, mock_seqvar):
    """Test the resolve_variant method for a sequence variant."""
    auto_acmg = AutoACMG("NM_000038.3:c.797G>A", GenomeRelease.GRCh38)
    variant = auto_acmg.resolve_variant()
    assert variant == mock_seqvar


def test_auto_acmg_resolve_sequence_variant_parse_error(mock_seqvar_resolver_parse_error):
    """Test the resolve_variant method for a sequence variant with a parse error."""
    auto_acmg = AutoACMG("NM_000038.3:c.797G>A", GenomeRelease.GRCh38)
    variant = auto_acmg.resolve_variant()
    assert variant is None


def test_auto_acmg_resolve_sequence_variant_failure(mock_seqvar_resolver_failure):
    """Test the resolve_variant method for a sequence variant that fails to resolve."""
    auto_acmg = AutoACMG("NM_000038.3:c.797G>A", GenomeRelease.GRCh38)
    variant = auto_acmg.resolve_variant()
    assert variant is None


def test_auto_acmg_predict_seqvar_success(
    mock_seqvar_resolver, mock_auto_pvs1_success, mock_auto_ps1_pm5_success, mock_seqvar
):
    """Test the predict method for a sequence variant with a successful response."""
    auto_acmg = AutoACMG("NM_000038.3:c.797G>A", GenomeRelease.GRCh38)
    with runner.isolated_filesystem():
        auto_acmg.predict()
    assert mock_auto_pvs1_success.predict.called
    assert mock_auto_ps1_pm5_success.predict.called


def test_auto_acmg_predict_seqvar_resolve_failure(
    mock_seqvar_resolver_failure, mock_auto_pvs1_failure
):
    """Test the predict method for a sequence variant with a failure response due to resolve method."""
    auto_acmg = AutoACMG("NM_000038.3:c.797G>A", GenomeRelease.GRCh38)
    with runner.isolated_filesystem():
        auto_acmg.predict()
    assert mock_seqvar_resolver_failure.resolve_seqvar.called
    assert not mock_auto_pvs1_failure.predict.called


def test_auto_acmg_predict_seqvar_pvs1_failure(
    mock_seqvar_resolver, mock_auto_pvs1_failure, mock_auto_ps1_pm5_success, mock_seqvar
):
    """Test the predict method for a sequence variant with a failure response."""
    auto_acmg = AutoACMG("NM_000038.3:c.797G>A", GenomeRelease.GRCh38)
    with runner.isolated_filesystem():
        auto_acmg.predict()
    assert mock_seqvar_resolver.resolve_seqvar.called
    assert mock_auto_pvs1_failure.predict.called
    assert auto_acmg.seqvar == mock_seqvar
    assert auto_acmg.seqvar_pvs1_prediction_path == PVS1PredictionSeqVarPath.NotSet


def test_auto_acmg_predict_seqvar_ps1_pm5_failure(
    mock_seqvar_resolver, mock_auto_pvs1_success, mock_auto_ps1_pm5_failure, mock_seqvar
):
    """Test the predict method for a sequence variant with a failure response."""
    auto_acmg = AutoACMG("NM_000038.3:c.797G>A", GenomeRelease.GRCh38)
    with runner.isolated_filesystem():
        auto_acmg.predict()
    assert mock_seqvar_resolver.resolve_seqvar.called
    assert mock_auto_pvs1_success.predict.called
    assert mock_auto_ps1_pm5_failure.predict.called
    assert auto_acmg.seqvar == mock_seqvar
    assert auto_acmg.seqvar_ps1 is None
    assert auto_acmg.seqvar_pm5 is None
