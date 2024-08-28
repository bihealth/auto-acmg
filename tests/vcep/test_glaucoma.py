from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.exceptions import MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import GlaucomaPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def glaucoma_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return GlaucomaPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_not_applicable(glaucoma_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for MYOC in Glaucoma."""
    auto_acmg_data.hgnc_id = "HGNC:7610"  # MYOC gene
    result = glaucoma_predictor.predict_pm1(glaucoma_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for MYOC."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting for MYOC."
    assert (
        "PM1 is not applicable for MYOC" in result.summary
    ), "The summary should indicate that PM1 is not applicable for MYOC."


def test_predict_pm1_strength_level(glaucoma_predictor, auto_acmg_data):
    """Test the strength level for PM1 when not applicable."""
    auto_acmg_data.hgnc_id = "HGNC:7610"  # MYOC gene
    result = glaucoma_predictor.predict_pm1(glaucoma_predictor.seqvar, auto_acmg_data)

    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting when PM1 is not applicable."


def test_predict_pm1_name(glaucoma_predictor, auto_acmg_data):
    """Test the name of the criteria returned by the Glaucoma predictor."""
    auto_acmg_data.hgnc_id = "HGNC:7610"  # MYOC gene
    result = glaucoma_predictor.predict_pm1(glaucoma_predictor.seqvar, auto_acmg_data)

    assert result.name == "PM1", "The name of the criteria should be 'PM1'."


def test_bp3_not_applicable(glaucoma_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = glaucoma_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


def test_is_conserved_with_valid_gerp_score(glaucoma_predictor, auto_acmg_data):
    """Test that the variant is correctly identified as conserved when the GERP score is valid."""
    auto_acmg_data.scores.cadd.gerp = 3.5  # Example GERP score above the threshold
    auto_acmg_data.thresholds.gerp = 2.0  # Set a lower threshold

    is_conserved = glaucoma_predictor._is_conserved(auto_acmg_data)

    assert is_conserved, "The variant should be identified as conserved with a high GERP score."


def test_is_conserved_with_low_gerp_score(glaucoma_predictor, auto_acmg_data):
    """Test that the variant is correctly identified as not conserved when the GERP score is low."""
    auto_acmg_data.scores.cadd.gerp = 1.5  # Example GERP score below the threshold
    auto_acmg_data.thresholds.gerp = 2.0  # Set a higher threshold

    is_conserved = glaucoma_predictor._is_conserved(auto_acmg_data)

    assert (
        not is_conserved
    ), "The variant should not be identified as conserved with a low GERP score."


def test_is_conserved_missing_gerp_score(glaucoma_predictor, auto_acmg_data):
    """Test that the method raises an error when the GERP score is missing."""
    auto_acmg_data.scores.cadd.gerp = None  # No GERP score

    with pytest.raises(MissingDataError, match="GERP score is missing."):
        glaucoma_predictor._is_conserved(auto_acmg_data)


def test_predict_bp7_threshold_adjustment(glaucoma_predictor, auto_acmg_data):
    """Test that the BP7 spliceAI thresholds are correctly adjusted for Glaucoma VCEP."""
    auto_acmg_data.thresholds.spliceAI_acceptor_gain = 0.05  # Initial threshold value
    auto_acmg_data.thresholds.spliceAI_acceptor_loss = 0.05  # Initial threshold value
    auto_acmg_data.thresholds.spliceAI_donor_gain = 0.05  # Initial threshold value
    auto_acmg_data.thresholds.spliceAI_donor_loss = 0.05  # Initial threshold value

    # Call predict_bp7 method
    result = glaucoma_predictor.predict_bp7(glaucoma_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
    ), "The BP7 acceptor gain threshold should be adjusted to 0.2."
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
    ), "The BP7 acceptor loss threshold should be adjusted to 0.2."
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
    ), "The BP7 donor gain threshold should be adjusted to 0.2."
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2
    ), "The BP7 donor loss threshold should be adjusted to 0.2."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(GlaucomaPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default(
    mock_super_predict_bp7, glaucoma_predictor, auto_acmg_data
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
    result = glaucoma_predictor.predict_bp7(glaucoma_predictor.seqvar, auto_acmg_data)

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
