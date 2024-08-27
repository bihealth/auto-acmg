from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.vhl import VHLPredictor


@pytest.fixture
def seqvar():
    return SeqVar(
        genome_release=GenomeRelease.GRCh37, chrom="3", pos=10191539, delete="C", insert="T"
    )


@pytest.fixture
def vhl_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return VHLPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_in_moderate_residue(vhl_predictor, auto_acmg_data):
    """Test when the variant affects a critical residue in VHL."""
    auto_acmg_data.hgnc_id = "HGNC:12687"  # VHL gene
    auto_acmg_data.prot_pos = 167  # A critical residue in VHL
    result = vhl_predictor.predict_pm1(vhl_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "affects a germline hotspot or key functional domain in VHL" in result.summary
    ), "The summary should indicate the PM1 criteria."


def test_predict_pm1_not_met(vhl_predictor, auto_acmg_data):
    """Test when the variant does not meet the PM1 criteria for VHL."""
    auto_acmg_data.hgnc_id = "HGNC:12687"  # VHL gene
    auto_acmg_data.prot_pos = 200  # Position outside of critical residues
    result = vhl_predictor.predict_pm1(vhl_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for non-critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not meet the PM1 criteria for VHL" in result.summary
    ), "The summary should indicate the PM1 criteria."


@patch("src.vcep.vhl.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, vhl_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the VHL VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default fallback for PM1.",
    )
    result = vhl_predictor.predict_pm1(vhl_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default fallback for PM1." in result.summary
    ), "The summary should indicate the default fallback."


def test_predict_bp7_threshold_adjustment(vhl_predictor, auto_acmg_data):
    """Test that the BP7 PhyloP threshold is correctly adjusted for VHL."""
    auto_acmg_data.thresholds.phyloP100 = 1.0  # Initial PhyloP threshold value

    # Call predict_bp7 method
    result = vhl_predictor.predict_bp7(vhl_predictor.seqvar, auto_acmg_data)

    # Check that the PhyloP threshold was adjusted
    assert (
        auto_acmg_data.thresholds.phyloP100 == 0.2
    ), "The BP7 PhyloP threshold should be adjusted to 0.2."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(VHLPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default(mock_super_predict_bp7, vhl_predictor, auto_acmg_data):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = vhl_predictor.predict_bp7(vhl_predictor.seqvar, auto_acmg_data)

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
