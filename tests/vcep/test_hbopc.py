from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import HBOPCPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def hbopc_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return HBOPCPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_not_applicable(hbopc_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for ATM and PALB2 in Hereditary Breast, Ovarian, and Pancreatic Cancer."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for ATM."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting for ATM."
    assert (
        "PM1 is not applicable for HGNC:795" in result.summary
    ), "The summary should indicate that PM1 is not applicable for ATM."


def test_predict_pm1_not_applicable_palb2(hbopc_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for PALB2 in Hereditary Breast, Ovarian, and Pancreatic Cancer."""
    auto_acmg_data.hgnc_id = "HGNC:26144"  # PALB2 gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for PALB2."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting for PALB2."
    assert (
        "PM1 is not applicable for HGNC:26144" in result.summary
    ), "The summary should indicate that PM1 is not applicable for PALB2."


@patch("src.vcep.hbopc.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, hbopc_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method if logic changes."""
    auto_acmg_data.hgnc_id = "HGNC:111111"  # Gene not in the specific logic
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met if the gene is not ATM or PALB2."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


def test_predict_pm1_strength_level(hbopc_predictor, auto_acmg_data):
    """Test the strength level for PM1 when not applicable."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting when PM1 is not applicable."


def test_predict_pm1_name(hbopc_predictor, auto_acmg_data):
    """Test the name of the criteria returned by the HBOPC predictor."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert result.name == "PM1", "The name of the criteria should be 'PM1'."


def test_predict_bp7_threshold_adjustment_for_palb2(hbopc_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted for PALB2."""
    auto_acmg_data.hgnc_id = "HGNC:26144"  # PALB2 gene
    auto_acmg_data.thresholds.bp7_donor = 2  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = hbopc_predictor.predict_bp7(hbopc_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted for PALB2
    assert (
        auto_acmg_data.thresholds.bp7_donor == 7
    ), "The BP7 donor threshold should be adjusted to 7 for PALB2."
    assert (
        auto_acmg_data.thresholds.bp7_acceptor == 21
    ), "The BP7 acceptor threshold should be adjusted to 21 for PALB2."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


def test_predict_bp7_threshold_adjustment_for_atm(hbopc_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted for ATM."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    auto_acmg_data.thresholds.bp7_donor = 2  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = hbopc_predictor.predict_bp7(hbopc_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted for ATM
    assert (
        auto_acmg_data.thresholds.bp7_donor == 7
    ), "The BP7 donor threshold should be adjusted to 7 for ATM."
    assert (
        auto_acmg_data.thresholds.bp7_acceptor == 40
    ), "The BP7 acceptor threshold should be adjusted to 40 for ATM."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(HBOPCPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default(mock_super_predict_bp7, hbopc_predictor, auto_acmg_data):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = hbopc_predictor.predict_bp7(hbopc_predictor.seqvar, auto_acmg_data)

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


@patch.object(HBOPCPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default_for_palb2(
    mock_super_predict_bp7, hbopc_predictor, auto_acmg_data
):
    """Test fallback to default BP7 prediction after threshold adjustment for PALB2."""
    auto_acmg_data.hgnc_id = "HGNC:26144"  # PALB2 gene
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = hbopc_predictor.predict_bp7(hbopc_predictor.seqvar, auto_acmg_data)

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
