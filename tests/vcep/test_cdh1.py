from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
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


@patch.object(DefaultPredictor, "verify_pm4bp3")
def test_predict_pm4bp3_cdh1(mock_verify_pm4bp3, cdh1_predictor, seqvar, auto_acmg_data):
    """Test the predict_pm4bp3 method for CDH1."""
    # Set the mock return value for the verify_pm4bp3 method
    mock_pred = MagicMock()
    mock_pred.PM4 = True
    mock_pred.BP3 = False
    mock_pred.PM4_strength = AutoACMGStrength.PathogenicSupporting
    mock_pred.BP3_strength = AutoACMGStrength.BenignSupporting
    mock_verify_pm4bp3.return_value = (mock_pred, "PM4 is met")

    # Call the method under test
    pm4_result, bp3_result = cdh1_predictor.predict_pm4bp3(seqvar, auto_acmg_data)

    # Check PM4 result
    assert isinstance(
        pm4_result, AutoACMGCriteria
    ), "The PM4 result should be of type AutoACMGCriteria."
    assert pm4_result.prediction == AutoACMGPrediction.Met, "PM4 should be Met."
    assert (
        pm4_result.strength == AutoACMGStrength.PathogenicSupporting
    ), "PM4 strength should be PathogenicSupporting."
    assert "PM4 is met" in pm4_result.summary, "The summary should indicate PM4 is met."

    # Check BP3 result
    assert isinstance(
        bp3_result, AutoACMGCriteria
    ), "The BP3 result should be of type AutoACMGCriteria."
    assert bp3_result.prediction == AutoACMGPrediction.NotApplicable, "BP3 should be NotApplicable."
    assert (
        bp3_result.strength == AutoACMGStrength.BenignSupporting
    ), "BP3 strength should be BenignSupporting."
    assert (
        "BP3 is not applicable for CDH1" in bp3_result.summary
    ), "The summary should indicate BP3 is not applicable."


@patch.object(DefaultPredictor, "verify_pm4bp3", return_value=(None, ""))
def test_predict_pm4bp3_fallback(mock_verify_pm4bp3, cdh1_predictor, seqvar, auto_acmg_data):
    """Test the fallback behavior for PM4 when verification fails."""
    # Call the method under test
    pm4_result, bp3_result = cdh1_predictor.predict_pm4bp3(seqvar, auto_acmg_data)

    # Check PM4 result
    assert (
        pm4_result.prediction == AutoACMGPrediction.Failed
    ), "PM4 should be Failed if verification fails."
    assert (
        "PM4 could not be verified." in pm4_result.summary
    ), "The summary should indicate PM4 could not be verified."

    # Check BP3 result
    assert bp3_result.prediction == AutoACMGPrediction.NotApplicable, "BP3 should be NotApplicable."
    assert (
        bp3_result.strength == AutoACMGStrength.BenignSupporting
    ), "BP3 strength should be BenignSupporting."
    assert (
        "BP3 is not applicable for CDH1" in bp3_result.summary
    ), "The summary should indicate BP3 is not applicable."


def test_predict_bp7_threshold_adjustment(cdh1_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted."""
    auto_acmg_data.thresholds.bp7_donor = 1  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = cdh1_predictor.predict_bp7(cdh1_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.bp7_donor == 7
    ), "The BP7 donor threshold should be adjusted to 7."
    assert (
        auto_acmg_data.thresholds.bp7_acceptor == 21
    ), "The BP7 acceptor threshold should be adjusted to 4."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(DefaultPredictor, "predict_bp7")
def test_predict_bp7_fallback_to_default(mock_super_predict_bp7, cdh1_predictor, auto_acmg_data):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = cdh1_predictor.predict_bp7(cdh1_predictor.seqvar, auto_acmg_data)

    # Verify the result and ensure the superclass method was called
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.NotMet, "BP7 should return NotMet as mocked."
    assert (
        result.strength == AutoACMGStrength.BenignSupporting
    ), "The strength should be BenignSupporting."
    assert (
        "Default BP7 prediction fallback." in result.summary
    ), "The summary should indicate the fallback."
    assert mock_super_predict_bp7.called, "super().predict_bp7 should have been called."
