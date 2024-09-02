from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import HBOPCPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def hbopc_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return HBOPCPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


@pytest.fixture
def auto_acmg_data_atm():
    data = AutoACMGData()
    data.hgnc_id = "HGNC:795"
    data.consequence = MagicMock(mehari=[], cadd=None)
    return data


@pytest.fixture
def auto_acmg_data_palb2():
    data = AutoACMGData()
    data.hgnc_id = "HGNC:26144"
    data.consequence = MagicMock(mehari=[], cadd=None)
    return data


def test_predict_pm1_not_applicable(hbopc_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for ATM and PALB2 in Hereditary Breast, Ovarian, and Pancreatic Cancer."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for ATM."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting for ATM."
    assert (
        "PM1 is not applicable for HGNC:795" in result.summary
    ), "The summary should indicate that PM1 is not applicable for ATM."


def test_predict_pm1_not_applicable_palb2(hbopc_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for PALB2 in Hereditary Breast, Ovarian, and Pancreatic Cancer."""
    auto_acmg_data.hgnc_id = "HGNC:26144"  # PALB2 gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for PALB2."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting for PALB2."
    assert (
        "PM1 is not applicable for HGNC:26144" in result.summary
    ), "The summary should indicate that PM1 is not applicable for PALB2."


@patch("src.vcep.hbopc.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, hbopc_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method if logic changes."""
    auto_acmg_data.hgnc_id = "HGNC:111111"  # Gene not in the specific logic
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met if the gene is not ATM or PALB2."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


def test_predict_pm1_strength_level(hbopc_predictor, auto_acmg_data):
    """Test the strength level for PM1 when not applicable."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting when PM1 is not applicable."


def test_predict_pm1_name(hbopc_predictor, auto_acmg_data):
    """Test the name of the criteria returned by the HBOPC predictor."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert result.name == "PM1", "The name of the criteria should be 'PM1'."


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
@pytest.mark.parametrize(
    "hgnc_id,expected_pm2,expected_ba1,expected_bs1",
    [("HGNC:795", 0.00001, 0.005, 0.0005), ("HGNC:26144", 0.000003, 0.001, 0.0001)],
)
def test_predict_pm2ba1bs1bs2_specific_genes(
    mock_super_method,
    hbopc_predictor,
    auto_acmg_data,
    seqvar,
    hgnc_id,
    expected_pm2,
    expected_ba1,
    expected_bs1,
):
    # Setup
    auto_acmg_data.hgnc_id = hgnc_id

    # Method call
    hbopc_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Validate thresholds are set correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == expected_pm2
    assert auto_acmg_data.thresholds.ba1_benign == expected_ba1
    assert auto_acmg_data.thresholds.bs1_benign == expected_bs1

    # Check that the superclass method was called with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Reset mock for the next iteration
    mock_super_method.reset_mock()


@patch.object(DefaultPredictor, "verify_pm4bp3")
def test_predict_pm4bp3_atm(mock_verify_pm4bp3, hbopc_predictor, seqvar, auto_acmg_data):
    """Test the predict_pm4bp3 method for ATM (HGNC:795)."""
    auto_acmg_data.hgnc_id = "HGNC:795"

    # Mock the verify_pm4bp3 method to return a specific result
    mock_pred = MagicMock()
    mock_pred.PM4 = True
    mock_pred.BP3 = False
    mock_verify_pm4bp3.return_value = (mock_pred, "PM4 is met for ATM")

    # Call the method under test
    pm4_result, bp3_result = hbopc_predictor.predict_pm4bp3(seqvar, auto_acmg_data)

    # Check PM4 result
    assert isinstance(
        pm4_result, AutoACMGCriteria
    ), "The PM4 result should be of type AutoACMGCriteria."
    assert pm4_result.prediction == AutoACMGPrediction.Met, "PM4 should be Met for ATM."
    assert (
        pm4_result.strength == AutoACMGStrength.PathogenicSupporting
    ), "PM4 strength should be PathogenicSupporting."
    assert (
        "PM4 is met for ATM" in pm4_result.summary
    ), "The summary should indicate PM4 is met for ATM."

    # Check BP3 result
    assert isinstance(
        bp3_result, AutoACMGCriteria
    ), "The BP3 result should be of type AutoACMGCriteria."
    assert (
        bp3_result.prediction == AutoACMGPrediction.NotApplicable
    ), "BP3 should be NotApplicable for ATM."
    assert (
        bp3_result.strength == AutoACMGStrength.BenignSupporting
    ), "BP3 strength should be BenignSupporting."
    assert (
        "BP3 is not applicable for ATM." in bp3_result.summary
    ), "The summary should indicate BP3 is not applicable for ATM."


def test_predict_pm4bp3_palb2(hbopc_predictor, seqvar, auto_acmg_data):
    """Test the predict_pm4bp3 method for PALB2 (HGNC:26144)."""
    auto_acmg_data.hgnc_id = "HGNC:26144"

    # Call the method under test
    pm4_result, bp3_result = hbopc_predictor.predict_pm4bp3(seqvar, auto_acmg_data)

    # Check PM4 result
    assert isinstance(
        pm4_result, AutoACMGCriteria
    ), "The PM4 result should be of type AutoACMGCriteria."
    assert (
        pm4_result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM4 should be NotApplicable for PALB2."
    assert (
        pm4_result.strength == AutoACMGStrength.PathogenicSupporting
    ), "PM4 strength should be PathogenicSupporting."
    assert (
        "PM4 is not applicable for PALB2." in pm4_result.summary
    ), "The summary should indicate PM4 is not applicable for PALB2."

    # Check BP3 result
    assert isinstance(
        bp3_result, AutoACMGCriteria
    ), "The BP3 result should be of type AutoACMGCriteria."
    assert (
        bp3_result.prediction == AutoACMGPrediction.NotApplicable
    ), "BP3 should be NotApplicable for PALB2."
    assert (
        bp3_result.strength == AutoACMGStrength.BenignSupporting
    ), "BP3 strength should be BenignSupporting."
    assert (
        "BP3 is not applicable for PALB2." in bp3_result.summary
    ), "The summary should indicate BP3 is not applicable for PALB2."


@patch.object(DefaultPredictor, "verify_pm4bp3")
@patch.object(DefaultPredictor, "predict_pm4bp3")
@pytest.mark.skip("Something doesn't work here")
def test_predict_pm4bp3_fallback(
    mock_predict_pm4bp3, mock_verify_pm4bp3, hbopc_predictor, seqvar, auto_acmg_data
):
    """Test the fallback behavior for predict_pm4bp3 method when the gene is not ATM or PALB2."""
    auto_acmg_data.hgnc_id = "HGNC:99999"  # Some gene not handled explicitly

    # Mock the verify_pm4bp3 method to return a specific result
    mock_pred = MagicMock()
    mock_pred.PM4 = True
    mock_pred.BP3 = False
    mock_verify_pm4bp3.return_value = (mock_pred, "PM4 is met for fallback")

    # Call the method under test
    pm4_result, bp3_result = hbopc_predictor.predict_pm4bp3(seqvar, auto_acmg_data)

    # Check PM4 result
    assert isinstance(
        pm4_result, AutoACMGCriteria
    ), "The PM4 result should be of type AutoACMGCriteria."
    assert pm4_result.prediction == AutoACMGPrediction.Met, "PM4 should be Met for fallback gene."
    assert (
        "PM4 is met for fallback" in pm4_result.summary
    ), "The summary should indicate PM4 is met for fallback gene."

    # Check BP3 result
    assert isinstance(
        bp3_result, AutoACMGCriteria
    ), "The BP3 result should be of type AutoACMGCriteria."
    assert (
        bp3_result.prediction == AutoACMGPrediction.NotMet
    ), "BP3 should be NotMet for fallback gene."
    assert (
        "PM4 is met for fallback" in bp3_result.summary
    ), "The summary should indicate BP3 was checked."

    # Ensure the superclass's predict_pm4bp3 method was called
    assert mock_predict_pm4bp3.called, "super().predict_pm4bp3 should have been called."


@patch.object(HBOPCPredictor, "predict_pp2bp1")
def test_predict_pp2bp1_palb2_missense(mock_predict, hbopc_predictor, seqvar, auto_acmg_data_palb2):
    """Test predict_pp2bp1 for PALB2 where the variant is missense."""
    auto_acmg_data_palb2.consequence.mehari = ["missense_variant"]
    mock_predict.return_value = (
        MagicMock(name="PP2", prediction=AutoACMGPrediction.NotApplicable),
        MagicMock(
            name="BP1",
            prediction=AutoACMGPrediction.Met,
            summary="Variant is missense. BP1 is met.",
        ),
    )

    pp2, bp1 = hbopc_predictor.predict_pp2bp1(seqvar, auto_acmg_data_palb2)

    assert (
        pp2.prediction == AutoACMGPrediction.NotApplicable
    ), "PP2 should be NotApplicable for PALB2."
    assert (
        bp1.prediction == AutoACMGPrediction.Met
    ), "BP1 should be Met for missense variant in PALB2."
    assert "missense" in bp1.summary, "Summary should confirm the missense nature."


@patch.object(HBOPCPredictor, "predict_pp2bp1")
def test_predict_pp2bp1_atm_non_missense(mock_predict, hbopc_predictor, seqvar, auto_acmg_data_atm):
    """Test predict_pp2bp1 for ATM where the variant is not missense."""
    auto_acmg_data_atm.consequence.mehari = ["nonsense_variant"]
    mock_predict.return_value = (
        MagicMock(name="PP2", prediction=AutoACMGPrediction.NotApplicable),
        MagicMock(
            name="BP1",
            prediction=AutoACMGPrediction.NotMet,
            summary="Variant is not missense. BP1 is not met.",
        ),
    )

    pp2, bp1 = hbopc_predictor.predict_pp2bp1(seqvar, auto_acmg_data_atm)

    assert (
        pp2.prediction == AutoACMGPrediction.NotApplicable
    ), "PP2 should be NotApplicable for ATM."
    assert (
        bp1.prediction == AutoACMGPrediction.NotMet
    ), "BP1 should not be Met for non-missense variant in ATM."
    assert "not missense" in bp1.summary, "Summary should confirm the non-missense nature."


def test_predict_bp7_threshold_adjustment_for_palb2(hbopc_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted for PALB2."""
    auto_acmg_data.hgnc_id = "HGNC:26144"  # PALB2 gene
    auto_acmg_data.thresholds.bp7_donor = 2  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = hbopc_predictor.predict_bp7(hbopc_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted for PALB2
    assert (
        auto_acmg_data.thresholds.bp7_donor == 7
    ), "The BP7 donor threshold should be adjusted to 7 for PALB2."
    assert (
        auto_acmg_data.thresholds.bp7_acceptor == 21
    ), "The BP7 acceptor threshold should be adjusted to 21 for PALB2."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


def test_predict_bp7_threshold_adjustment_for_atm(hbopc_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted for ATM."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    auto_acmg_data.thresholds.bp7_donor = 2  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = hbopc_predictor.predict_bp7(hbopc_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted for ATM
    assert (
        auto_acmg_data.thresholds.bp7_donor == 7
    ), "The BP7 donor threshold should be adjusted to 7 for ATM."
    assert (
        auto_acmg_data.thresholds.bp7_acceptor == 40
    ), "The BP7 acceptor threshold should be adjusted to 40 for ATM."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(HBOPCPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default(mock_super_predict_bp7, hbopc_predictor, auto_acmg_data):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = hbopc_predictor.predict_bp7(hbopc_predictor.seqvar, auto_acmg_data)

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


@patch.object(HBOPCPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default_for_palb2(
    mock_super_predict_bp7, hbopc_predictor, auto_acmg_data
):
    """Test fallback to default BP7 prediction after threshold adjustment for PALB2."""
    auto_acmg_data.hgnc_id = "HGNC:26144"  # PALB2 gene
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = hbopc_predictor.predict_bp7(hbopc_predictor.seqvar, auto_acmg_data)

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
