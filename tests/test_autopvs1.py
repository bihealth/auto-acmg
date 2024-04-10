from unittest.mock import Mock

import pytest
from typer.testing import CliRunner

from src.autoPVS1 import AutoPVS1
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar, SeqVarResolver
from src.seqvar_pvs1 import PVS1Prediction, PVS1PredictionSeqVarPath, SeqVarPVS1

runner = CliRunner()


@pytest.fixture
def mock_seqvar_resolver(monkeypatch):
    mock_resolver = Mock(SeqVarResolver)
    mock_resolver.resolve_seqvar.return_value = SeqVar(
        genome_release=GenomeRelease.GRCh38, chrom="1", pos=100000, delete="A", insert="T"
    )
    monkeypatch.setattr("src.autoPVS1.SeqVarResolver", lambda: mock_resolver)
    return mock_resolver


@pytest.fixture
def mock_seqvar_pvs1(monkeypatch):
    mock_pvs1 = Mock(SeqVarPVS1)
    mock_pvs1.initialize.return_value = None
    mock_pvs1.verify_PVS1.return_value = None
    mock_pvs1.get_prediction.return_value = (PVS1Prediction.PVS1, PVS1PredictionSeqVarPath.NF1)
    monkeypatch.setattr("src.autoPVS1.SeqVarPVS1", lambda *args, **kwargs: mock_pvs1)
    return mock_pvs1


def test_autoPVS1_resolve_variant_success(mock_seqvar_resolver):
    """Test resolve_variant method with a successful response."""
    autoPVS1 = AutoPVS1("GRCh38-1-100000-A-T")
    with runner.isolated_filesystem():
        result = autoPVS1.resolve_variant()
    assert isinstance(result, SeqVar)
    # Assert mock_seqvar_resolver.resolve_seqvar was called
    assert mock_seqvar_resolver.resolve_seqvar.called


def test_autoPVS1_predict_success(mock_seqvar_resolver, mock_seqvar_pvs1):
    """Test predict method with a successful response."""
    autoPVS1 = AutoPVS1("GRCh38-1-100000-A-T")
    with runner.isolated_filesystem():
        autoPVS1.predict()
    # Assert mock_seqvar_pvs1.initialize was called
    assert mock_seqvar_pvs1.initialize.called
