from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import CerebralCreatineDeficiencySyndromesPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def cerebral_creatine_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return CerebralCreatineDeficiencySyndromesPredictor(
        seqvar=seqvar, result=result, config=MagicMock()
    )


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_not_applicable(cerebral_creatine_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for Cerebral Creatine Deficiency Syndromes."""
    result = cerebral_creatine_predictor.predict_pm1(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for Cerebral Creatine Deficiency Syndromes."
    assert (
        result.summary == "PM1 is not applicable for Cerebral Creatine Deficiency Syndromes."
    ), "The summary should indicate that PM1 is not applicable for Cerebral Creatine Deficiency Syndromes."


def test_predict_pm1_strength(cerebral_creatine_predictor, auto_acmg_data):
    """Test the strength level returned by the Cerebral Creatine Deficiency Syndromes predictor."""
    result = cerebral_creatine_predictor.predict_pm1(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate for Cerebral Creatine Deficiency Syndromes."


def test_predict_pm1_name(cerebral_creatine_predictor, auto_acmg_data):
    """Test the name of the criteria returned by the Cerebral Creatine Deficiency Syndromes predictor."""
    result = cerebral_creatine_predictor.predict_pm1(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    assert result.name == "PM1", "The name of the criteria should be 'PM1'."


@patch("src.vcep.cerebral_creatine_deficiency_syndromes.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, cerebral_creatine_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method (if implemented)."""
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = cerebral_creatine_predictor.predict_pm1(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    # In this specific case, the fallback should never happen since PM1 is always not applicable,
    # but this test ensures that if something changes, the fallback works correctly.
    assert result.prediction == AutoACMGPrediction.NotApplicable, "PM1 should remain NotApplicable."


def test_bp3_not_applicable(cerebral_creatine_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = cerebral_creatine_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


def test_predict_bp7_threshold_adjustment(cerebral_creatine_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted."""
    auto_acmg_data.thresholds.bp7_donor = 1  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value
    # Set the hgnc_id to GAMT (HGNC:4136) to trigger the threshold adjustment
    auto_acmg_data.hgnc_id = "HGNC:4136"

    # Call predict_bp7 method
    result = cerebral_creatine_predictor.predict_bp7(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

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
def test_predict_bp7_fallback_to_default(
    mock_super_predict_bp7, cerebral_creatine_predictor, auto_acmg_data
):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = cerebral_creatine_predictor.predict_bp7(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

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
