from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import (
    PS1PM5,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
)
from src.defs.exceptions import MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor
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
    return AutoACMGSeqVarData()


def test_predict_pvs1_not_applicable(glaucoma_predictor, seqvar, auto_acmg_data):
    result = glaucoma_predictor.predict_pvs1(seqvar, auto_acmg_data)

    # Verify the outcome is as expected, always not applicable
    assert result.name == "PVS1"
    assert result.prediction == AutoACMGPrediction.NotApplicable
    assert result.strength == AutoACMGStrength.PathogenicVeryStrong
    assert result.summary == "PVS1 is not applicable for the gene."


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_overrides(mock_super_verify, glaucoma_predictor, seqvar, auto_acmg_data):
    """Test that the overridden verify_ps1pm5 method in GlaucomaPredictor works correctly."""
    # Set up the mock to return PS1 and PM5 as applicable initially
    mock_super_verify.return_value = (PS1PM5(PS1=True, PM5=True), "Initial evaluation")

    # Setup the data
    auto_acmg_data.consequence = MagicMock(mehari=["missense_variant"])
    auto_acmg_data.thresholds = MagicMock(
        spliceAI_acceptor_gain=0.5,
        spliceAI_acceptor_loss=0.5,
        spliceAI_donor_gain=0.5,
        spliceAI_donor_loss=0.5,
    )
    auto_acmg_data.scores = MagicMock(
        cadd=MagicMock(spliceAI_acceptor_gain=0.3, spliceAI_donor_gain=0.3)
    )

    # Run the method under test
    prediction, comment = glaucoma_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that the splicing effect leads to overriding PS1 and PM5 as not applicable
    assert not prediction.PS1, "PS1 should be marked as not applicable due to splicing effect."
    assert not prediction.PM5, "PM5 should be marked as not applicable due to splicing effect."
    assert "Variant affects splicing" in comment, "Comment should note the splicing effect."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
    ), "SpliceAI acceptor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
    ), "SpliceAI acceptor loss threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
    ), "SpliceAI donor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2
    ), "SpliceAI donor loss threshold should be adjusted to 0.2"


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_no_splicing_effect(
    mock_super_verify, glaucoma_predictor, seqvar, auto_acmg_data
):
    """Test that PS1 and PM5 remain applicable when there's no splicing effect."""
    # Set up the mock to return PS1 and PM5 as applicable initially
    mock_super_verify.return_value = (PS1PM5(PS1=True, PM5=True), "Initial evaluation")

    # Setup the data with no splicing effect
    auto_acmg_data.consequence = MagicMock(mehari=["missense_variant"])
    auto_acmg_data.thresholds = MagicMock(
        spliceAI_acceptor_gain=0.5,
        spliceAI_acceptor_loss=0.5,
        spliceAI_donor_gain=0.5,
        spliceAI_donor_loss=0.5,
    )
    auto_acmg_data.scores = MagicMock(
        cadd=MagicMock(
            spliceAI_acceptor_gain=0.1,
            spliceAI_acceptor_loss=0.1,
            spliceAI_donor_gain=0.1,
            spliceAI_donor_loss=0.1,
        )
    )

    # Run the method under test
    prediction, comment = glaucoma_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that PS1 and PM5 remain applicable
    assert prediction.PS1, "PS1 should remain applicable when there's no splicing effect."
    assert prediction.PM5, "PM5 should remain applicable when there's no splicing effect."
    assert "Initial evaluation" in comment, "Comment should reflect the initial evaluation."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
    ), "SpliceAI acceptor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
    ), "SpliceAI acceptor loss threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
    ), "SpliceAI donor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2
    ), "SpliceAI donor loss threshold should be adjusted to 0.2"


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_non_missense(mock_super_verify, glaucoma_predictor, seqvar, auto_acmg_data):
    """Test that PS1 and PM5 remain as per superclass for non-missense variants."""
    # Set up the mock to return PS1 and PM5 as not applicable initially
    mock_super_verify.return_value = (
        PS1PM5(PS1=False, PM5=False),
        "Not applicable for non-missense",
    )

    # Setup the data with a non-missense variant
    auto_acmg_data.consequence = MagicMock(mehari=["synonymous_variant"])

    # Run the method under test
    prediction, comment = glaucoma_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that PS1 and PM5 remain as per superclass evaluation
    assert not prediction.PS1, "PS1 should remain not applicable for non-missense variants."
    assert not prediction.PM5, "PM5 should remain not applicable for non-missense variants."
    assert (
        "Not applicable for non-missense" in comment
    ), "Comment should reflect superclass evaluation."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
    ), "SpliceAI acceptor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
    ), "SpliceAI acceptor loss threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
    ), "SpliceAI donor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2
    ), "SpliceAI donor loss threshold should be adjusted to 0.2"


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_exception_handling(
    mock_super_verify, glaucoma_predictor, seqvar, auto_acmg_data
):
    """Test verify_ps1pm5 method exception handling."""
    mock_super_verify.side_effect = Exception("Test exception")

    with pytest.raises(Exception) as exc_info:
        glaucoma_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    assert "Test exception" in str(exc_info.value), "Should raise the original exception"


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


def test_bs2_not_applicable(glaucoma_predictor, auto_acmg_data):
    """Test when BS2 is not applicable for Glaucoma VCEP."""
    result = glaucoma_predictor._bs2_not_applicable(auto_acmg_data)
    assert result is True, "BS2 should always be not applicable for Glaucoma VCEP."


@patch.object(
    DefaultSeqVarPredictor,
    "predict_pm2ba1bs1bs2",
    return_value=(
        AutoACMGCriteria(name="PM2"),
        AutoACMGCriteria(name="BA1"),
        AutoACMGCriteria(name="BS1"),
        AutoACMGCriteria(name="BS2"),
    ),
)
def test_predict_pm2ba1bs1bs2(mock_super_method, glaucoma_predictor, auto_acmg_data, seqvar):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = glaucoma_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.0001
    assert auto_acmg_data.thresholds.ba1_benign == 0.01
    assert auto_acmg_data.thresholds.bs1_benign == 0.001

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


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


def test_predict_pp2bp1(glaucoma_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Glaucoma predictor."""

    # Call the method under test
    pp2_result, bp1_result = glaucoma_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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
