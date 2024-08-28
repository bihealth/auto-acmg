from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.pten import PTENPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="10", pos=100, delete="A", insert="T")


@pytest.fixture
def pten_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return PTENPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_in_catalytic_motifs(pten_predictor, auto_acmg_data):
    """Test when the variant falls within the catalytic motifs of PTEN."""
    auto_acmg_data.hgnc_id = "HGNC:9588"  # PTEN gene
    auto_acmg_data.prot_pos = 92  # Within the catalytic motif (90-94)
    result = pten_predictor.predict_pm1(pten_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for a variant in the catalytic motif."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "catalytic motifs of PTEN" in result.summary
    ), "The summary should indicate the catalytic motif."


def test_predict_pm1_outside_catalytic_motifs(pten_predictor, auto_acmg_data):
    """Test when the variant does not fall within the catalytic motifs of PTEN."""
    auto_acmg_data.hgnc_id = "HGNC:9588"  # PTEN gene
    auto_acmg_data.prot_pos = 150  # Outside the catalytic motifs
    result = pten_predictor.predict_pm1(pten_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for a variant outside the catalytic motifs."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant does not meet the PM1 criteria for PTEN" in result.summary
    ), "The summary should indicate no criteria were met."


@patch("src.vcep.pten.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, pten_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the PTEN VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = pten_predictor.predict_pm1(pten_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


def test_bp3_not_applicable(pten_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = pten_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


def test_verify_pm4bp3_stop_loss(pten_predictor, seqvar, auto_acmg_data):
    """Test verify_pm4bp3 when the variant is a stop-loss mutation."""
    auto_acmg_data.consequence.cadd = "stop_loss"

    # Call the method under test
    prediction, comment = pten_predictor.verify_pm4bp3(seqvar, auto_acmg_data)

    assert prediction.PM4 is True, "PM4 should be met for stop-loss variants."
    assert prediction.BP3 is False, "BP3 should not be met for stop-loss variants."
    assert "Variant consequence is stop-loss. PM4 is met." in comment


@pytest.mark.skip("THis test should work")
def test_verify_pm4bp3_inframe_delins_in_catalytic_motif(pten_predictor, seqvar, auto_acmg_data):
    """Test verify_pm4bp3 when the variant is an in-frame indel in the catalytic motif."""
    auto_acmg_data.consequence.cadd = "inframe_deletion"
    auto_acmg_data.prot_pos = 167  # Assume this is within the catalytic motif

    # Call the method under test
    prediction, comment = pten_predictor.verify_pm4bp3(seqvar, auto_acmg_data)

    assert prediction.PM4 is True, "PM4 should be met for in-frame indels in the catalytic motif."
    assert (
        prediction.BP3 is False
    ), "BP3 should not be met for in-frame indels in the catalytic motif."
    assert "Impacting catalytic motif. PM4 is met." in comment


def test_verify_pm4bp3_inframe_delins_outside_catalytic_motif(
    pten_predictor, seqvar, auto_acmg_data
):
    """Test verify_pm4bp3 when the variant is an in-frame indel outside the catalytic motif."""
    auto_acmg_data.consequence.cadd = "inframe_deletion"
    auto_acmg_data.prot_pos = 200  # Assume this is outside the catalytic motif

    # Call the method under test
    prediction, comment = pten_predictor.verify_pm4bp3(seqvar, auto_acmg_data)

    assert (
        prediction.PM4 is False
    ), "PM4 should not be met for in-frame indels outside the catalytic motif."
    assert (
        prediction.BP3 is False
    ), "BP3 should not be met for in-frame indels outside the catalytic motif."
    assert "No impact on catalytic motif. PM4 is not met." in comment


def test_verify_pm4bp3_neither_indel_nor_stop_loss(pten_predictor, seqvar, auto_acmg_data):
    """Test verify_pm4bp3 when the variant is neither an in-frame indel nor a stop-loss."""
    auto_acmg_data.consequence.cadd = "missense_variant"

    # Call the method under test
    prediction, comment = pten_predictor.verify_pm4bp3(seqvar, auto_acmg_data)

    assert prediction.PM4 is False, "PM4 should not be met for non-indel, non-stop-loss variants."
    assert prediction.BP3 is False, "BP3 should not be met for non-indel, non-stop-loss variants."
    assert "consequence is not stop" in comment


@patch.object(PTENPredictor, "verify_pp2bp1")
def test_predict_pp2bp1_pten(mock_verify, pten_predictor, seqvar, auto_acmg_data):
    """Test PP2 and BP1 prediction for PTEN."""
    # Setting up the mock to return a specific result for the PP2 prediction
    mock_verify.return_value = (
        MagicMock(PP2=True, PP2_strength=AutoACMGStrength.PathogenicSupporting),
        "PP2 met for PTEN variant.",
    )

    auto_acmg_data.hgnc_id = "HGNC:9588"  # PTEN gene

    pp2, bp1 = pten_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert pp2.prediction == AutoACMGPrediction.Met, "PP2 should be Met for PTEN."
    assert (
        pp2.strength == AutoACMGStrength.PathogenicSupporting
    ), "PP2 strength should be PathogenicSupporting."
    assert "PP2 met for PTEN variant" in pp2.summary, "PP2 summary should confirm the met criteria."
    assert (
        bp1.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should be NotApplicable for PTEN."


@patch.object(PTENPredictor, "verify_pp2bp1")
def test_predict_pp2bp1_pten_failure(mock_verify, pten_predictor, seqvar, auto_acmg_data):
    """Test when PP2 prediction fails for PTEN."""
    # Simulating a failure scenario
    mock_verify.return_value = (None, "PP2 prediction failed due to an error.")

    auto_acmg_data.hgnc_id = "HGNC:9588"  # PTEN gene

    pp2, bp1 = pten_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert pp2.prediction == AutoACMGPrediction.Failed, "PP2 prediction should fail."
    assert (
        "PP2 prediction failed due to an error" in pp2.summary
    ), "PP2 summary should report the failure."
    assert (
        bp1.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should still be NotApplicable for PTEN."


def test_predict_bp7_threshold_adjustment(pten_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted."""
    auto_acmg_data.thresholds.bp7_donor = 1  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = pten_predictor.predict_bp7(pten_predictor.seqvar, auto_acmg_data)

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
def test_predict_bp7_fallback_to_default(mock_super_predict_bp7, pten_predictor, auto_acmg_data):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = pten_predictor.predict_bp7(pten_predictor.seqvar, auto_acmg_data)

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
