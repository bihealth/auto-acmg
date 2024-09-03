from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import CardiomyopathyPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def cardiomyopathy_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object, which will be passed along
    return CardiomyopathyPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


@patch.object(DefaultPredictor, "predict_pvs1")
def test_predict_pvs1_not_applicable(
    mock_super_predict_pvs1, cardiomyopathy_predictor, seqvar, auto_acmg_data
):
    # Set the HGNC ID to one that makes PVS1 not applicable
    auto_acmg_data.hgnc_id = "HGNC:7577"
    result = cardiomyopathy_predictor.predict_pvs1(seqvar, auto_acmg_data)

    # Verify the outcome is as expected for the specific HGNC ID
    assert result.name == "PVS1", "The criterion name should be 'PVS1'"
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "The prediction should be NotApplicable"
    assert (
        result.strength == AutoACMGStrength.PathogenicVeryStrong
    ), "The strength should be PathogenicVeryStrong"
    assert (
        result.summary == "PVS1 is not applicable for the gene."
    ), "The summary should indicate that PVS1 is not applicable"
    mock_super_predict_pvs1.assert_not_called()  # Ensure superclass method is not called


@patch.object(
    DefaultPredictor,
    "predict_pvs1",
    return_value=AutoACMGCriteria(
        name="PVS1",
        prediction=AutoACMGPrediction.Met,
        strength=AutoACMGStrength.PathogenicVeryStrong,
        summary="Default behavior.",
    ),
)
def test_predict_pvs1_calls_superclass(
    mock_super_predict_pvs1, cardiomyopathy_predictor, seqvar, auto_acmg_data
):
    # Set the HGNC ID to one not affecting the PVS1 applicability
    auto_acmg_data.hgnc_id = "HGNC:Random"
    result = cardiomyopathy_predictor.predict_pvs1(seqvar, auto_acmg_data)

    # Verify that the superclass method is called
    mock_super_predict_pvs1.assert_called_once_with(seqvar, auto_acmg_data)
    # Check the response from the superclass method
    assert result.name == "PVS1"
    assert result.prediction == AutoACMGPrediction.Met
    assert result.strength == AutoACMGStrength.PathogenicVeryStrong
    assert result.summary == "Default behavior."


def test_predict_pm1_in_critical_domain(cardiomyopathy_predictor, auto_acmg_data):
    """Test when the variant falls within a critical domain for MYH7."""
    auto_acmg_data.prot_pos = 500  # Within the critical domain (167-931) for MYH7
    auto_acmg_data.hgnc_id = "HGNC:7577"  # MYH7 gene
    result = cardiomyopathy_predictor.predict_pm1(cardiomyopathy_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for a variant in a critical domain."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."


def test_predict_pm1_outside_critical_domain(cardiomyopathy_predictor, auto_acmg_data):
    """Test when the variant does not fall within any critical domain."""
    auto_acmg_data.prot_pos = 1000  # Outside the critical domain (167-931) for MYH7
    auto_acmg_data.hgnc_id = "HGNC:7577"  # MYH7 gene
    result = cardiomyopathy_predictor.predict_pm1(cardiomyopathy_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for a variant outside any critical domain."
    assert (
        "does not fall within any critical domain" in result.summary
    ), "The summary should indicate the lack of critical domain."


def test_predict_pm1_not_applicable(cardiomyopathy_predictor, auto_acmg_data):
    """Test when the PM1 criterion is not applicable."""
    auto_acmg_data.hgnc_id = "HGNC:12010"  # Example for TPM1, where PM1 is not applicable
    result = cardiomyopathy_predictor.predict_pm1(cardiomyopathy_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for TPM1."
    assert "Not applicable" in result.summary, "The summary should indicate non-applicability."


@patch("src.vcep.cardiomyopathy.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, cardiomyopathy_predictor, auto_acmg_data
):
    """Test when the transcript ID is not in the PM1_CLUSTER mapping, it should fallback to default PM1 prediction."""
    auto_acmg_data.hgnc_id = "HGNC:111111111111111"  # Not in the PM1_CLUSTER mapping
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = cardiomyopathy_predictor.predict_pm1(cardiomyopathy_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        result.summary == "Default PM1 prediction fallback."
    ), "The summary should be from the default fallback."


def test_predict_pm1_edge_case_start_boundary(cardiomyopathy_predictor, auto_acmg_data):
    """Test when variant falls exactly on the start boundary of a critical domain."""
    auto_acmg_data.prot_pos = 167  # Start boundary of the MYH7 critical domain (167-931)
    auto_acmg_data.hgnc_id = "HGNC:7577"  # MYH7 gene
    result = cardiomyopathy_predictor.predict_pm1(cardiomyopathy_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the start boundary of a critical domain."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."


def test_predict_pm1_edge_case_end_boundary(cardiomyopathy_predictor, auto_acmg_data):
    """Test when variant falls exactly on the end boundary of a critical domain."""
    auto_acmg_data.prot_pos = 931  # End boundary of the MYH7 critical domain (167-931)
    auto_acmg_data.hgnc_id = "HGNC:7577"  # MYH7 gene
    result = cardiomyopathy_predictor.predict_pm1(cardiomyopathy_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the end boundary of a critical domain."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."


def test_bs2_not_applicable(cardiomyopathy_predictor, auto_acmg_data):
    """Test BS2 is not applicable for Cardiomyopathy."""
    result = cardiomyopathy_predictor._bs2_not_applicable(auto_acmg_data)
    assert result is True, "BS2 should always be not applicable for Cardiomyopathy."


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
def test_predict_pm2ba1bs1bs2_myh7(
    mock_super_method, cardiomyopathy_predictor, auto_acmg_data, seqvar
):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01
    auto_acmg_data.hgnc_id = "HGNC:7551"  # MYH7 gene

    result = cardiomyopathy_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.00004
    assert auto_acmg_data.thresholds.ba1_benign == 0.001
    assert auto_acmg_data.thresholds.bs1_benign == 0.0002

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


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
def test_predict_pm2ba1bs1bs2_tpm1(
    mock_super_method, cardiomyopathy_predictor, auto_acmg_data, seqvar
):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01
    auto_acmg_data.hgnc_id = "HGNC:12010"

    result = cardiomyopathy_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.00004
    assert auto_acmg_data.thresholds.ba1_benign == 0.001
    assert auto_acmg_data.thresholds.bs1_benign == 0.0001

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_bp3_not_applicable(cardiomyopathy_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = cardiomyopathy_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


def test_predict_pp2bp1(cardiomyopathy_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for CardiomyopathyPredictor."""

    # Call the method under test
    pp2_result, bp1_result = cardiomyopathy_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


def test_predict_bp7_threshold_adjustment(cardiomyopathy_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted for Cardiomyopathy."""
    auto_acmg_data.thresholds.bp7_donor = 2  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = cardiomyopathy_predictor.predict_bp7(cardiomyopathy_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.bp7_donor == 7
    ), "The BP7 donor threshold should be adjusted to 7."
    assert (
        auto_acmg_data.thresholds.bp7_acceptor == 4
    ), "The BP7 acceptor threshold should be adjusted to 4."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(DefaultPredictor, "predict_bp7")
def test_predict_bp7_fallback_to_default(
    mock_super_predict_bp7, cardiomyopathy_predictor, auto_acmg_data
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
    result = cardiomyopathy_predictor.predict_bp7(cardiomyopathy_predictor.seqvar, auto_acmg_data)

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
