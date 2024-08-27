from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.rett_angelman import RettAngelmanPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="18", pos=100, delete="A", insert="T")


@pytest.fixture
def rett_angelman_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return RettAngelmanPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_in_critical_region(rett_angelman_predictor, auto_acmg_data):
    """Test when the variant falls within a critical region for a Rett and Angelman-like Disorder gene."""
    auto_acmg_data.hgnc_id = "HGNC:11634"  # TCF4 gene
    auto_acmg_data.prot_pos = 580  # Within the critical bHLH domain (564-617)
    result = rett_angelman_predictor.predict_pm1(rett_angelman_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for critical region variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "critical residue in HGNC:11634" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_not_applicable(rett_angelman_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for SLC9A6."""
    auto_acmg_data.hgnc_id = "HGNC:11079"  # SLC9A6 gene
    result = rett_angelman_predictor.predict_pm1(rett_angelman_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be applicable for SLC9A6."
    assert "not applicable" in result.summary, "The summary should indicate non-applicability."


def test_predict_pm1_outside_critical_region(rett_angelman_predictor, auto_acmg_data):
    """Test when the variant does not fall within any critical region for Rett and Angelman-like Disorders genes."""
    auto_acmg_data.hgnc_id = "HGNC:11634"  # TCF4 gene
    auto_acmg_data.prot_pos = 700  # Position outside all critical regions
    result = rett_angelman_predictor.predict_pm1(rett_angelman_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for non-critical region variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not meet the PM1 criteria" in result.summary
    ), "The summary should indicate no critical region."


@patch("src.vcep.rett_angelman.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, rett_angelman_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the Rett and Angelman-like VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = rett_angelman_predictor.predict_pm1(rett_angelman_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


def test_predict_bp7_threshold_adjustment(rett_angelman_predictor, auto_acmg_data):
    """Test that the BP7 thresholds are correctly adjusted for Rett and Angelman-like Disorders."""
    auto_acmg_data.thresholds.phyloP100 = 1.0  # Initial phyloP100 threshold value

    # Call predict_bp7 method
    result = rett_angelman_predictor.predict_bp7(rett_angelman_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.phyloP100 == 0.1
    ), "The phyloP100 threshold should be adjusted to 0.1."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(RettAngelmanPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default(
    mock_super_predict_bp7, rett_angelman_predictor, auto_acmg_data
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
    result = rett_angelman_predictor.predict_bp7(rett_angelman_predictor.seqvar, auto_acmg_data)

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
