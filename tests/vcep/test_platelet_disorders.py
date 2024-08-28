from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.platelet_disorders import PlateletDisordersPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="17", pos=100, delete="A", insert="T")


@pytest.fixture
def platelet_disorders_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return PlateletDisordersPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_not_met_for_itga2b(platelet_disorders_predictor, auto_acmg_data):
    """Test when the variant is in the ITGA2B gene and PM1 is not met."""
    auto_acmg_data.hgnc_id = "HGNC:6138"  # ITGA2B gene
    result = platelet_disorders_predictor.predict_pm1(
        platelet_disorders_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.NotMet, "PM1 should not be met for ITGA2B."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not met for ITGA2B and ITGB3" in result.summary
    ), "The summary should explain the reason."


def test_predict_pm1_not_met_for_itgb3(platelet_disorders_predictor, auto_acmg_data):
    """Test when the variant is in the ITGB3 gene and PM1 is not met."""
    auto_acmg_data.hgnc_id = "HGNC:6156"  # ITGB3 gene
    result = platelet_disorders_predictor.predict_pm1(
        platelet_disorders_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.NotMet, "PM1 should not be met for ITGB3."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not met for ITGA2B and ITGB3" in result.summary
    ), "The summary should explain the reason."


@patch("src.vcep.platelet_disorders.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, platelet_disorders_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the Platelet Disorders VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = platelet_disorders_predictor.predict_pm1(
        platelet_disorders_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


def test_predict_pp2bp1(platelet_disorders_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Platelet Disorders."""

    # Call the method under test
    pp2_result, bp1_result = platelet_disorders_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    # Check PP2 result
    assert isinstance(
        pp2_result, AutoACMGCriteria
    ), "The PP2 result should be of type AutoACMGCriteria."
    assert (
        pp2_result.prediction == AutoACMGPrediction.NotApplicable
    ), "PP2 should be NotApplicable for ACADVL."
    assert (
        pp2_result.summary == "PP2 is not applicable for the gene."
    ), "The summary should indicate PP2 is not applicable."

    # Check BP1 result
    assert isinstance(
        bp1_result, AutoACMGCriteria
    ), "The BP1 result should be of type AutoACMGCriteria."
    assert (
        bp1_result.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should be NotApplicable for ACADVL."
    assert (
        bp1_result.summary == "BP1 is not applicable for the gene."
    ), "The summary should indicate BP1 is not applicable."
