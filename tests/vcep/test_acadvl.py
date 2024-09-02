from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import ACADVLPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="17", pos=100, delete="A", insert="T")


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


@pytest.fixture
def acadvl_predictor(seqvar, auto_acmg_data):
    result = MagicMock()  # Mocking the AutoACMGResult object, which will be passed along
    return ACADVLPredictor(seqvar=seqvar, result=result, config=MagicMock())


def test_predict_pm1_in_critical_region(acadvl_predictor, auto_acmg_data):
    """Test when variant falls within a critical region for ACADVL."""
    auto_acmg_data.prot_pos = 220  # Set protein position within a critical region (214-223)
    result = acadvl_predictor.predict_pm1(acadvl_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for critical region variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "falls within a critical region" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_outside_critical_region(acadvl_predictor, auto_acmg_data):
    """Test when variant does not fall within any critical region for ACADVL."""
    auto_acmg_data.prot_pos = 300  # Set protein position outside all critical regions
    result = acadvl_predictor.predict_pm1(acadvl_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for non-critical region variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not fall within any critical region" in result.summary
    ), "The summary should indicate no critical region."


def test_predict_pm1_edge_case_start_boundary(acadvl_predictor, auto_acmg_data):
    """Test when variant falls exactly on the start boundary of a critical region."""
    auto_acmg_data.prot_pos = 214  # Start boundary of the first critical region (214-223)
    result = acadvl_predictor.predict_pm1(acadvl_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the start boundary of a critical region."
    assert (
        "falls within a critical region" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_edge_case_end_boundary(acadvl_predictor, auto_acmg_data):
    """Test when variant falls exactly on the end boundary of a critical region."""
    auto_acmg_data.prot_pos = 223  # End boundary of the first critical region (214-223)
    result = acadvl_predictor.predict_pm1(acadvl_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the end boundary of a critical region."
    assert (
        "falls within a critical region" in result.summary
    ), "The summary should indicate the critical region."


@patch.object(
    DefaultPredictor,
    "predict_pm2ba1bs1bs2",
    return_value=(
        AutoACMGCriteria(name="PM2"),
        AutoACMGCriteria(name="BA1"),
        AutoACMGCriteria(name="BS1"),
        AutoACMGCriteria(name="BS2"),
    ),
)
def test_predict_pm2ba1bs1bs2(mock_super_method, acadvl_predictor, auto_acmg_data, seqvar):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = acadvl_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.001
    assert auto_acmg_data.thresholds.ba1_benign == 0.007
    assert auto_acmg_data.thresholds.bs1_benign == 0.0035

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_bp3_not_applicable_acadvl(acadvl_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = acadvl_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable for ACADVL."


def test_predict_pp2bp1(acadvl_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for ACADVL predictor."""

    # Call the method under test
    pp2_result, bp1_result = acadvl_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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
