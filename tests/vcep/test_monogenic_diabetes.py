from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.monogenic_diabetes import MonogenicDiabetesPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="12", pos=100, delete="A", insert="T")


@pytest.fixture
def monogenic_diabetes_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return MonogenicDiabetesPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


@pytest.fixture
def auto_acmg_data_gck():
    data = AutoACMGData(hgnc_id="HGNC:4195")
    data.consequence = MagicMock(mehari=[], cadd=None)
    return data


def test_predict_pm1_moderate_criteria_hnf1a(monogenic_diabetes_predictor, auto_acmg_data):
    """Test when variant falls within a moderate level residue for HNF1A."""
    auto_acmg_data.hgnc_id = "HGNC:11621"  # HNF1A gene
    auto_acmg_data.prot_pos = 130  # Within the critical residues for HNF1A
    result = monogenic_diabetes_predictor.predict_pm1(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.Met, "PM1 should be met for critical residues."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "critical for Monogenic Diabetes" in result.summary
    ), "The summary should indicate the critical residue."


def test_predict_pm1_supporting_criteria_hnf4a_domain(monogenic_diabetes_predictor, auto_acmg_data):
    """Test when variant falls within a supporting level domain for HNF4A."""
    auto_acmg_data.hgnc_id = "HGNC:5024"  # HNF4A gene
    auto_acmg_data.prot_pos = 100  # Within the critical domain for HNF4A
    result = monogenic_diabetes_predictor.predict_pm1(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.Met, "PM1 should be met for critical domains."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "falls within a critical region" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_supporting_criteria_hnf4a_promoter(
    monogenic_diabetes_predictor, auto_acmg_data
):
    """Test when variant falls within a supporting level promoter region for HNF4A."""
    auto_acmg_data.hgnc_id = "HGNC:5024"  # HNF4A gene
    auto_acmg_data.cds_pos = -140  # Within the critical promoter region for HNF4A
    result = monogenic_diabetes_predictor.predict_pm1(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for critical promoter regions."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "falls within a critical promoter region" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_moderate_criteria_gck(monogenic_diabetes_predictor, auto_acmg_data):
    """Test when variant falls within a moderate level residue for GCK."""
    auto_acmg_data.hgnc_id = "HGNC:4195"  # GCK gene
    auto_acmg_data.prot_pos = 151  # Within the critical residues for GCK
    result = monogenic_diabetes_predictor.predict_pm1(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.Met, "PM1 should be met for critical residues."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "critical for Monogenic Diabetes" in result.summary
    ), "The summary should indicate the critical residue."


def test_predict_pm1_not_met(monogenic_diabetes_predictor, auto_acmg_data):
    """Test when variant does not meet any PM1 criteria for Monogenic Diabetes genes."""
    auto_acmg_data.hgnc_id = "HGNC:4195"  # GCK gene
    auto_acmg_data.prot_pos = 500  # Outside all critical regions
    result = monogenic_diabetes_predictor.predict_pm1(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
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


@patch("src.vcep.monogenic_diabetes.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, monogenic_diabetes_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the PM1_CLUSTER
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = monogenic_diabetes_predictor.predict_pm1(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


def test_bp3_not_applicable(monogenic_diabetes_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = monogenic_diabetes_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


@patch.object(MonogenicDiabetesPredictor, "_is_missense")
def test_predict_pp2bp1_gck_missense(
    mock_is_missense, monogenic_diabetes_predictor, seqvar, auto_acmg_data_gck
):
    """Test predict_pp2bp1 for GCK where the variant is missense."""
    mock_is_missense.return_value = True
    pp2, bp1 = monogenic_diabetes_predictor.predict_pp2bp1(seqvar, auto_acmg_data_gck)

    assert (
        pp2.prediction == AutoACMGPrediction.Met
    ), "PP2 should be Met for a missense variant in GCK."
    assert (
        pp2.strength == AutoACMGStrength.PathogenicSupporting
    ), "PP2 strength should be PathogenicSupporting."
    assert "missense variant" in pp2.summary, "PP2 summary should confirm the missense nature."
    assert (
        bp1.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should be NotApplicable for GCK."


@patch.object(MonogenicDiabetesPredictor, "_is_missense")
def test_predict_pp2bp1_gck_non_missense(
    mock_is_missense, monogenic_diabetes_predictor, seqvar, auto_acmg_data_gck
):
    """Test predict_pp2bp1 for GCK where the variant is not missense."""
    mock_is_missense.return_value = False
    pp2, bp1 = monogenic_diabetes_predictor.predict_pp2bp1(seqvar, auto_acmg_data_gck)

    assert (
        pp2.prediction == AutoACMGPrediction.NotMet
    ), "PP2 should be NotMet for non-missense variants in GCK."
    assert (
        "not a missense variant" in pp2.summary
    ), "PP2 summary should confirm the non-missense nature."
    assert (
        bp1.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should still be NotApplicable for GCK."


def test_predict_bp7_threshold_adjustment(monogenic_diabetes_predictor, auto_acmg_data):
    """Test that the BP7 thresholds are correctly adjusted for Monogenic Diabetes."""
    # Initial threshold values
    auto_acmg_data.thresholds.spliceAI_acceptor_gain = 0.1
    auto_acmg_data.thresholds.spliceAI_acceptor_loss = 0.1
    auto_acmg_data.thresholds.spliceAI_donor_gain = 0.1
    auto_acmg_data.thresholds.spliceAI_donor_loss = 0.1
    auto_acmg_data.thresholds.phyloP100 = 1.0

    # Call predict_bp7 method
    result = monogenic_diabetes_predictor.predict_bp7(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
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


@patch.object(MonogenicDiabetesPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default(
    mock_super_predict_bp7, monogenic_diabetes_predictor, auto_acmg_data
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
    result = monogenic_diabetes_predictor.predict_bp7(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
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
