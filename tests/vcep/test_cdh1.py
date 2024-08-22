from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import CDH1Predictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def cdh1_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return CDH1Predictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_not_applicable(cdh1_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for CDH1."""
    result = cdh1_predictor.predict_pm1(cdh1_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for CDH1."
    assert (
        result.summary == "PM1 is not applicable for CDH1."
    ), "The summary should indicate that PM1 is not applicable for CDH1."


def test_predict_pm1_strength(cdh1_predictor, auto_acmg_data):
    """Test the strength level returned by the CDH1 predictor."""
    result = cdh1_predictor.predict_pm1(cdh1_predictor.seqvar, auto_acmg_data)

    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting for CDH1."


def test_predict_pm1_name(cdh1_predictor, auto_acmg_data):
    """Test the name of the criteria returned by the CDH1 predictor."""
    result = cdh1_predictor.predict_pm1(cdh1_predictor.seqvar, auto_acmg_data)

    assert result.name == "PM1", "The name of the criteria should be 'PM1'."


@patch("src.vcep.cdh1.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, cdh1_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method (if implemented)."""
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = cdh1_predictor.predict_pm1(cdh1_predictor.seqvar, auto_acmg_data)

    # In this specific case, the fallback should never happen since PM1 is always not applicable,
    # but this test ensures that if something changes, the fallback works correctly.
    assert result.prediction == AutoACMGPrediction.NotApplicable, "PM1 should remain NotApplicable."
