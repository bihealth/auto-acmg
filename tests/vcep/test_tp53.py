from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.tp53 import TP53Predictor


@pytest.fixture
def seqvar():
    return SeqVar(
        genome_release=GenomeRelease.GRCh37, chrom="17", pos=7578406, delete="C", insert="T"
    )


@pytest.fixture
def tp53_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return TP53Predictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_in_moderate_residue(tp53_predictor, auto_acmg_data):
    """Test when the variant affects a critical residue in TP53."""
    auto_acmg_data.hgnc_id = "HGNC:11998"  # TP53 gene
    auto_acmg_data.prot_pos = 175  # A critical residue in TP53
    result = tp53_predictor.predict_pm1(tp53_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "affects a critical residue in TP53" in result.summary
    ), "The summary should indicate the PM1 criteria."


def test_predict_pm1_not_met(tp53_predictor, auto_acmg_data):
    """Test when the variant does not meet the PM1 criteria for TP53."""
    auto_acmg_data.hgnc_id = "HGNC:11998"  # TP53 gene
    auto_acmg_data.prot_pos = 300  # Position outside of critical residues
    result = tp53_predictor.predict_pm1(tp53_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for non-critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not meet the PM1 criteria for TP53" in result.summary
    ), "The summary should indicate the PM1 criteria."


@patch("src.vcep.tp53.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, tp53_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the TP53 VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default fallback for PM1.",
    )
    result = tp53_predictor.predict_pm1(tp53_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."


def test_predict_pm4bp3_tp53(tp53_predictor, seqvar, auto_acmg_data):
    """Test the predict_pm4bp3 method for TP53 VCEP."""
    # Call the method under test
    pm4_result, bp3_result = tp53_predictor.predict_pm4bp3(seqvar, auto_acmg_data)

    # Check PM4 result
    assert isinstance(
        pm4_result, AutoACMGCriteria
    ), "The PM4 result should be of type AutoACMGCriteria."
    assert (
        pm4_result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM4 should be NotApplicable for TP53."
    assert (
        "PM4 is not applicable for TP53 VCEP." in pm4_result.summary
    ), "The summary should indicate PM4 is not applicable for TP53."

    # Check BP3 result
    assert isinstance(
        bp3_result, AutoACMGCriteria
    ), "The BP3 result should be of type AutoACMGCriteria."
    assert (
        bp3_result.prediction == AutoACMGPrediction.NotApplicable
    ), "BP3 should be NotApplicable for TP53."
    assert (
        "BP3 is not applicable for TP53 VCEP." in bp3_result.summary
    ), "The summary should indicate BP3 is not applicable for TP53."


def test_predict_bp7_threshold_adjustment(tp53_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted."""
    auto_acmg_data.thresholds.bp7_donor = 1  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = tp53_predictor.predict_bp7(tp53_predictor.seqvar, auto_acmg_data)

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
def test_predict_bp7_fallback_to_default(mock_super_predict_bp7, tp53_predictor, auto_acmg_data):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = tp53_predictor.predict_bp7(tp53_predictor.seqvar, auto_acmg_data)

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
