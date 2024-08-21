from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import BrainMalformationsPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


@pytest.fixture
def brain_malformations_predictor(seqvar, auto_acmg_data):
    result = MagicMock()  # Mocking the AutoACMGResult object, which will be passed along
    return BrainMalformationsPredictor(seqvar=seqvar, result=result, config=MagicMock())


def test_predict_pm1_in_critical_domain(brain_malformations_predictor, auto_acmg_data):
    """Test when variant falls within a critical domain for the specified gene."""
    auto_acmg_data.prot_pos = 100  # Set protein position within a critical domain
    auto_acmg_data.transcript_id = "NM_005465.4"  # AKT3 transcript ID
    result = brain_malformations_predictor.predict_pm1(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for critical domain variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."


def test_predict_pm1_outside_critical_domain(brain_malformations_predictor, auto_acmg_data):
    """Test when variant does not fall within any critical domain for the specified gene."""
    auto_acmg_data.prot_pos = 500  # Set protein position outside all critical domains
    auto_acmg_data.transcript_id = "NM_005465.4"  # AKT3 transcript ID
    result = brain_malformations_predictor.predict_pm1(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for non-critical domain variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "does not fall within any critical domain" in result.summary
    ), "The summary should indicate no critical domain."


def test_predict_pm1_edge_case_start_boundary(brain_malformations_predictor, auto_acmg_data):
    """Test when variant falls exactly on the start boundary of a critical domain."""
    auto_acmg_data.prot_pos = 1382  # Start boundary of the MTOR kinase domain (1382-1982)
    auto_acmg_data.transcript_id = "NM_004958.3_MTOR"  # MTOR transcript ID
    result = brain_malformations_predictor.predict_pm1(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the start boundary of a critical domain."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."


def test_predict_pm1_edge_case_end_boundary(brain_malformations_predictor, auto_acmg_data):
    """Test when variant falls exactly on the end boundary of a critical domain."""
    auto_acmg_data.prot_pos = 1982  # End boundary of the MTOR kinase domain (1382-1982)
    auto_acmg_data.transcript_id = "NM_004958.3_MTOR"  # MTOR transcript ID
    result = brain_malformations_predictor.predict_pm1(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the end boundary of a critical domain."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."


@patch.object(DefaultPredictor, "predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_super_predict_pm1, brain_malformations_predictor, auto_acmg_data
):
    """Test when the transcript ID is not in the PM1_CLUSTER mapping, it should fallback to default PM1 prediction."""
    auto_acmg_data.transcript_id = "NM_XXXXXX.4_AAA"  # A non-existent transcript ID for fallback

    # Set the mock return value for the superclass's predict_pm1 method
    mock_super_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )

    result = brain_malformations_predictor.predict_pm1(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met if no specific cluster mapping is found."
    assert mock_super_predict_pm1.called, "super().predict_pm1 should have been called."
