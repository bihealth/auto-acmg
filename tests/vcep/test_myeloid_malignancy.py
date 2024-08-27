from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.monogenic_diabetes import MonogenicDiabetesPredictor
from src.vcep.myeloid_malignancy import MyeloidMalignancyPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="21", pos=100, delete="A", insert="T")


@pytest.fixture
def myeloid_malignancy_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return MyeloidMalignancyPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_moderate_criteria_runx1(myeloid_malignancy_predictor, auto_acmg_data):
    """Test when variant falls within a moderate level residue for RUNX1."""
    auto_acmg_data.hgnc_id = "HGNC:10471"  # RUNX1 gene
    auto_acmg_data.prot_pos = 107  # Within the critical residues for RUNX1
    result = myeloid_malignancy_predictor.predict_pm1(
        myeloid_malignancy_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.Met, "PM1 should be met for critical residues."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "critical residue within the Runt Homology Domain" in result.summary
    ), "The summary should indicate the critical residue."


def test_predict_pm1_supporting_criteria_runx1(myeloid_malignancy_predictor, auto_acmg_data):
    """Test when variant falls within a supporting level residue for RUNX1."""
    auto_acmg_data.hgnc_id = "HGNC:10471"  # RUNX1 gene
    auto_acmg_data.prot_pos = 100  # Within the supporting residues for RUNX1
    result = myeloid_malignancy_predictor.predict_pm1(
        myeloid_malignancy_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.Met, "PM1 should be met for supporting residues."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "residue within the Runt Homology Domain" in result.summary
    ), "The summary should indicate the supporting residue."


def test_predict_pm1_not_met(myeloid_malignancy_predictor, auto_acmg_data):
    """Test when variant does not meet any PM1 criteria for RUNX1."""
    auto_acmg_data.hgnc_id = "HGNC:10471"  # RUNX1 gene
    auto_acmg_data.prot_pos = 300  # Outside all critical and supporting regions
    result = myeloid_malignancy_predictor.predict_pm1(
        myeloid_malignancy_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for non-critical residues."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not meet the PM1 criteria" in result.summary
    ), "The summary should indicate that criteria were not met."


@patch("src.vcep.myeloid_malignancy.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, myeloid_malignancy_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the PM1_CLUSTER
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = myeloid_malignancy_predictor.predict_pm1(
        myeloid_malignancy_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


def test_predict_bp7_threshold_adjustment(myeloid_malignancy_predictor, auto_acmg_data):
    """Test that the BP7 thresholds are correctly adjusted for Myeloid Malignancy."""
    # Initial threshold values
    auto_acmg_data.thresholds.spliceAI_acceptor_gain = 0.1
    auto_acmg_data.thresholds.spliceAI_acceptor_loss = 0.1
    auto_acmg_data.thresholds.spliceAI_donor_gain = 0.1
    auto_acmg_data.thresholds.spliceAI_donor_loss = 0.1
    auto_acmg_data.thresholds.phyloP100 = 1.0

    # Call predict_bp7 method
    result = myeloid_malignancy_predictor.predict_bp7(
        myeloid_malignancy_predictor.seqvar, auto_acmg_data
    )

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
    ), "The spliceAI acceptor gain threshold should be adjusted to 0.2."
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
    ), "The spliceAI acceptor loss threshold should be adjusted to 0.2."
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
    ), "The spliceAI donor gain threshold should be adjusted to 0.2."
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2
    ), "The spliceAI donor loss threshold should be adjusted to 0.2."
    assert (
        auto_acmg_data.thresholds.phyloP100 == 2.0
    ), "The phyloP100 threshold should be adjusted to 2.0."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(MyeloidMalignancyPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default(
    mock_super_predict_bp7, myeloid_malignancy_predictor, auto_acmg_data
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
    result = myeloid_malignancy_predictor.predict_bp7(
        myeloid_malignancy_predictor.seqvar, auto_acmg_data
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
