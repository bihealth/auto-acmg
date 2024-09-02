from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import HearingLossPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def hearing_loss_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return HearingLossPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_kcnq4_strong(hearing_loss_predictor, auto_acmg_data):
    """Test when PM1 is met at the Strong level for a variant in KCNQ4."""
    auto_acmg_data.hgnc_id = "HGNC:6298"  # KCNQ4 gene
    auto_acmg_data.prot_pos = 280  # Within the critical pore-forming intramembrane region
    result = hearing_loss_predictor.predict_pm1(hearing_loss_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.Met, "PM1 should be met at the Strong level."
    assert (
        result.strength == AutoACMGStrength.PathogenicStrong
    ), "The strength should be PathogenicStrong."
    assert (
        "critical pore-forming intramembrane region" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_not_applicable(hearing_loss_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for other hearing loss genes."""
    auto_acmg_data.hgnc_id = "HGNC:4284"  # GJB2 gene
    result = hearing_loss_predictor.predict_pm1(hearing_loss_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for GJB2."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not applicable" in result.summary
    ), "The summary should indicate non-applicability."


def test_predict_pm1_not_met(hearing_loss_predictor, auto_acmg_data):
    """Test when PM1 is not met for KCNQ4 but outside the critical region."""
    auto_acmg_data.hgnc_id = "HGNC:6298"  # KCNQ4 gene
    auto_acmg_data.prot_pos = 300  # Outside the critical pore-forming intramembrane region
    result = hearing_loss_predictor.predict_pm1(hearing_loss_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.NotMet, "PM1 should not be met for KCNQ4."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant does not meet the PM1 criteria" in result.summary
    ), "The summary should indicate no criteria were met."


@patch("src.vcep.hearing_loss.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, hearing_loss_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method if logic changes."""
    auto_acmg_data.hgnc_id = "HGNC:111111"  # Gene not in the specific logic
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = hearing_loss_predictor.predict_pm1(hearing_loss_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


def test_predict_pm1_edge_case_start_boundary(hearing_loss_predictor, auto_acmg_data):
    """Test when variant falls exactly on the start boundary of the critical region."""
    auto_acmg_data.hgnc_id = "HGNC:6298"  # KCNQ4 gene
    auto_acmg_data.prot_pos = (
        271  # Start boundary of the critical pore-forming intramembrane region
    )
    result = hearing_loss_predictor.predict_pm1(hearing_loss_predictor.seqvar, auto_acmg_data)

    assert result.prediction == AutoACMGPrediction.Met, "PM1 should be met on the start boundary."
    assert (
        "critical pore-forming intramembrane region" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_edge_case_end_boundary(hearing_loss_predictor, auto_acmg_data):
    """Test when variant falls exactly on the end boundary of the critical region."""
    auto_acmg_data.hgnc_id = "HGNC:6298"  # KCNQ4 gene
    auto_acmg_data.prot_pos = 292  # End boundary of the critical pore-forming intramembrane region
    result = hearing_loss_predictor.predict_pm1(hearing_loss_predictor.seqvar, auto_acmg_data)

    assert result.prediction == AutoACMGPrediction.Met, "PM1 should be met on the end boundary."
    assert (
        "critical pore-forming intramembrane region" in result.summary
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
def test_predict_pm2ba1bs1bs2(mock_super_method, hearing_loss_predictor, auto_acmg_data, seqvar):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = hearing_loss_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.00007
    assert auto_acmg_data.thresholds.ba1_benign == 0.001
    assert auto_acmg_data.thresholds.bs1_benign == 0.0007

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_predict_pp2bp1(hearing_loss_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Hearing Loss predictor."""

    # Call the method under test
    pp2_result, bp1_result = hearing_loss_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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
