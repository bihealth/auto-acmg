from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
    GenomicStrand,
)
from src.defs.exceptions import AlgorithmError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import CoagulationFactorDeficiencyPredictor


@pytest.fixture
def seqvar():
    return SeqVar(
        genome_release=GenomeRelease.GRCh37, chrom="X", pos=100000, delete="A", insert="T"
    )


@pytest.fixture
def coagulation_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return CoagulationFactorDeficiencyPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    data = AutoACMGData()
    data.hgnc_id = "HGNC:3546"  # F8
    data.prot_pos = 391  # Test position
    data.exons = [MagicMock(altStartI=1, altEndI=1000000)]
    data.strand = GenomicStrand.Plus
    return data


def test_predict_pm1_strong_criteria(coagulation_predictor, auto_acmg_data):
    """Test when the variant falls within strong criteria for F8."""
    auto_acmg_data.prot_pos = 391  # Within the strong criteria for F8
    result = coagulation_predictor.predict_pm1(coagulation_predictor.seqvar, auto_acmg_data)

    assert result.prediction == AutoACMGPrediction.Met, "PM1 should be met at the Strong level."
    assert (
        result.strength == AutoACMGStrength.PathogenicStrong
    ), "The strength should be PathogenicStrong."
    assert (
        "PM1 is met at the Strong level" in result.summary
    ), "The summary should indicate the strong criteria."


def test_predict_pm1_moderate_criteria_residues(coagulation_predictor, auto_acmg_data):
    """Test when the variant falls within moderate criteria for F8 based on residues."""
    auto_acmg_data.prot_pos = 1667  # Within the moderate criteria for F8
    result = coagulation_predictor.predict_pm1(coagulation_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met at the Moderate level for residues."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "PM1 is met at the Moderate level" in result.summary
    ), "The summary should indicate the moderate criteria for residues."


def test_predict_pm1_moderate_criteria_regions(coagulation_predictor, auto_acmg_data):
    """Test when the variant falls within moderate criteria for F8 based on regions."""
    auto_acmg_data.prot_pos = 2270  # Within the FXa-binding region for F8
    result = coagulation_predictor.predict_pm1(coagulation_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met at the Moderate level for regions."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "PM1 is met at the Moderate level" in result.summary
    ), "The summary should indicate the moderate criteria for regions."


def test_predict_pm1_moderate_criteria_exons(coagulation_predictor, auto_acmg_data):
    """Test when the variant falls within moderate criteria for F9 based on exons."""
    auto_acmg_data.hgnc_id = "HGNC:3551"  # F9
    auto_acmg_data.strand = GenomicStrand.Plus
    auto_acmg_data.prot_pos = 1
    auto_acmg_data.exons = [
        MagicMock(altStartI=1, altEndI=200),
        MagicMock(altStartI=201, altEndI=400),
        MagicMock(altStartI=401, altEndI=600),
        MagicMock(altStartI=601, altEndI=1000000),
    ]

    result = coagulation_predictor.predict_pm1(coagulation_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met at the Moderate level for exons."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "PM1 is met at the Moderate level" in result.summary
    ), "The summary should indicate the moderate criteria for exons."


def test_predict_pm1_no_criteria_met(coagulation_predictor, auto_acmg_data):
    """Test when the variant does not meet any PM1 criteria."""
    auto_acmg_data.prot_pos = 5000  # Position outside any defined criteria
    result = coagulation_predictor.predict_pm1(coagulation_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met when no criteria are matched."
    assert (
        "Variant does not meet the PM1 criteria" in result.summary
    ), "The summary should indicate no criteria were met."


def test_predict_pm1_invalid_strand(coagulation_predictor, auto_acmg_data):
    """Test when an invalid strand is provided."""
    auto_acmg_data.strand = "invalid_strand"

    with pytest.raises(AlgorithmError):
        coagulation_predictor._get_affected_exon(auto_acmg_data, coagulation_predictor.seqvar)


def test_predict_bp7_threshold_adjustment_for_hgnc_3546(coagulation_predictor, auto_acmg_data):
    """Test that the BP7 thresholds are correctly adjusted for HGNC:3546 (F5)."""
    auto_acmg_data.hgnc_id = "HGNC:3546"  # F5 gene

    # Call predict_bp7 method
    result = coagulation_predictor.predict_bp7(coagulation_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.05
    ), "The spliceAI acceptor gain threshold should be adjusted to 0.05."
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.05
    ), "The spliceAI acceptor loss threshold should be adjusted to 0.05."
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.05
    ), "The spliceAI donor gain threshold should be adjusted to 0.05."
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.05
    ), "The spliceAI donor loss threshold should be adjusted to 0.05."
    assert (
        auto_acmg_data.thresholds.phyloP100 == 0.1
    ), "The phyloP100 threshold should be adjusted to 0.1."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


def test_predict_bp7_threshold_adjustment_for_hgnc_3551(coagulation_predictor, auto_acmg_data):
    """Test that the BP7 thresholds are correctly adjusted for HGNC:3551 (F7)."""
    auto_acmg_data.hgnc_id = "HGNC:3551"  # F7 gene

    # Call predict_bp7 method
    result = coagulation_predictor.predict_bp7(coagulation_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.01
    ), "The spliceAI acceptor gain threshold should be adjusted to 0.01."
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.01
    ), "The spliceAI acceptor loss threshold should be adjusted to 0.01."
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.01
    ), "The spliceAI donor gain threshold should be adjusted to 0.01."
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.01
    ), "The spliceAI donor loss threshold should be adjusted to 0.01."
    assert (
        auto_acmg_data.thresholds.phyloP100 == 0.1
    ), "The phyloP100 threshold should be adjusted to 0.1."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(DefaultPredictor, "predict_bp7")
def test_predict_bp7_fallback_to_default(
    mock_super_predict_bp7, coagulation_predictor, auto_acmg_data
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
    result = coagulation_predictor.predict_bp7(coagulation_predictor.seqvar, auto_acmg_data)

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
