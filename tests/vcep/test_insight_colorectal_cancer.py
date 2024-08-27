from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.insight_colorectal_cancer import InsightColorectalCancerPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="5", pos=100, delete="A", insert="T")


@pytest.fixture
def insight_colorectal_cancer_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return InsightColorectalCancerPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_not_applicable_apc(insight_colorectal_cancer_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for APC."""
    auto_acmg_data.hgnc_id = "HGNC:583"  # APC gene
    result = insight_colorectal_cancer_predictor.predict_pm1(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for APC."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not applicable for HGNC:583" in result.summary
    ), "The summary should indicate PM1 is not applicable."


def test_predict_pm1_not_applicable_mlh1(insight_colorectal_cancer_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for MLH1."""
    auto_acmg_data.hgnc_id = "HGNC:7127"  # MLH1 gene
    result = insight_colorectal_cancer_predictor.predict_pm1(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for MLH1."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not applicable for HGNC:7127" in result.summary
    ), "The summary should indicate PM1 is not applicable."


def test_predict_pm1_not_applicable_msh2(insight_colorectal_cancer_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for MSH2."""
    auto_acmg_data.hgnc_id = "HGNC:7325"  # MSH2 gene
    result = insight_colorectal_cancer_predictor.predict_pm1(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for MSH2."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not applicable for HGNC:7325" in result.summary
    ), "The summary should indicate PM1 is not applicable."


def test_predict_pm1_not_applicable_msh6(insight_colorectal_cancer_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for MSH6."""
    auto_acmg_data.hgnc_id = "HGNC:7329"  # MSH6 gene
    result = insight_colorectal_cancer_predictor.predict_pm1(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for MSH6."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not applicable for HGNC:7329" in result.summary
    ), "The summary should indicate PM1 is not applicable."


def test_predict_pm1_not_applicable_pms2(insight_colorectal_cancer_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for PMS2."""
    auto_acmg_data.hgnc_id = "HGNC:9122"  # PMS2 gene
    result = insight_colorectal_cancer_predictor.predict_pm1(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for PMS2."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not applicable for HGNC:9122" in result.summary
    ), "The summary should indicate PM1 is not applicable."


@patch("src.vcep.insight_colorectal_cancer.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, insight_colorectal_cancer_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method if logic changes."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the specific logic
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = insight_colorectal_cancer_predictor.predict_pm1(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


def test_predict_bp7_threshold_adjustment(insight_colorectal_cancer_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted."""
    auto_acmg_data.thresholds.bp7_donor = 1  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = insight_colorectal_cancer_predictor.predict_bp7(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
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
    mock_super_predict_bp7, insight_colorectal_cancer_predictor, auto_acmg_data
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
    result = insight_colorectal_cancer_predictor.predict_bp7(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
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
