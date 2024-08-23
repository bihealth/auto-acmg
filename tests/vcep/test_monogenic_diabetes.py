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
