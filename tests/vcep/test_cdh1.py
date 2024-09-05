from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    GenomicStrand,
)
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor
from src.vcep import CDH1Predictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def cdh1_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return CDH1Predictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


def test_predict_ps1pm5_not_applicable(cdh1_predictor, seqvar, auto_acmg_data):
    """Test PS1 is always not applicable and PM5 is evaluated based on specific conditions."""
    auto_acmg_data.consequence = MagicMock(mehari=["frameshift_variant"])
    auto_acmg_data.tx_pos_utr = 10
    auto_acmg_data.hgnc_id = "HGNC:1748"
    auto_acmg_data.strand = GenomicStrand.Plus
    auto_acmg_data.exons = [
        MagicMock(altStartI=1, altEndI=100, altCdsStartI=1, altCdsEndI=100),
        MagicMock(altStartI=200, altEndI=300, altCdsStartI=200, altCdsEndI=300),
    ]
    # Need for undergo_nmd
    cdh1_predictor.comment_pvs1 = ""
    ps1, pm5 = cdh1_predictor.predict_ps1pm5(seqvar, auto_acmg_data)

    # Check PS1
    assert ps1.name == "PS1"
    assert ps1.prediction == AutoACMGPrediction.NotApplicable
    assert ps1.strength == AutoACMGStrength.PathogenicSupporting
    assert "PS1 is not applicable for CDH1." in ps1.summary

    # Check PM5
    assert pm5.name == "PM5"
    assert pm5.prediction == AutoACMGPrediction.Met
    assert pm5.strength == AutoACMGStrength.PathogenicSupporting
    assert "Nonsense or frameshift variant predicted to undergo NMD. PM5 is met." in pm5.summary


def test_predict_ps1pm5_not_met(cdh1_predictor, seqvar, auto_acmg_data):
    """Test PM5 is not met if the variant does not lead to NMD."""
    auto_acmg_data.consequence = MagicMock(mehari=["stop_gained"])
    auto_acmg_data.tx_pos_utr = 5000
    auto_acmg_data.hgnc_id = "HGNC:1748"
    auto_acmg_data.strand = GenomicStrand.Plus
    auto_acmg_data.exons = [
        MagicMock(altStartI=1, altEndI=100, altCdsStartI=1, altCdsEndI=100),
        MagicMock(altStartI=200, altEndI=300, altCdsStartI=200, altCdsEndI=300),
    ]
    # Need for undergo_nmd
    cdh1_predictor.comment_pvs1 = ""
    ps1, pm5 = cdh1_predictor.predict_ps1pm5(seqvar, auto_acmg_data)

    # Check PM5 not met
    assert pm5.prediction == AutoACMGPrediction.NotMet
    assert (
        "Consequence is not frameshift or nonsense or variant is not predicted to undergo NMD. PM5 is not met."
        in pm5.summary
    )


def test_predict_pm1_not_applicable(cdh1_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for CDH1."""
    result = cdh1_predictor.predict_pm1(cdh1_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for CDH1."
    assert (
        result.summary == "PM1 is not applicable for CDH1."
    ), "The summary should indicate that PM1 is not applicable for CDH1."


def test_predict_pm1_strength(cdh1_predictor, auto_acmg_data):
    """Test the strength level returned by the CDH1 predictor."""
    result = cdh1_predictor.predict_pm1(cdh1_predictor.seqvar, auto_acmg_data)

    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting for CDH1."


def test_predict_pm1_name(cdh1_predictor, auto_acmg_data):
    """Test the name of the criteria returned by the CDH1 predictor."""
    result = cdh1_predictor.predict_pm1(cdh1_predictor.seqvar, auto_acmg_data)

    assert result.name == "PM1", "The name of the criteria should be 'PM1'."


@patch("src.vcep.cdh1.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, cdh1_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method (if implemented)."""
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = cdh1_predictor.predict_pm1(cdh1_predictor.seqvar, auto_acmg_data)

    # In this specific case, the fallback should never happen since PM1 is always not applicable,
    # but this test ensures that if something changes, the fallback works correctly.
    assert result.prediction == AutoACMGPrediction.NotApplicable, "PM1 should remain NotApplicable."


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
def test_predict_pm2ba1bs1bs2(mock_super_method, cdh1_predictor, auto_acmg_data, seqvar):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = cdh1_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.00001
    assert auto_acmg_data.thresholds.ba1_benign == 0.002
    assert auto_acmg_data.thresholds.bs1_benign == 0.001

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


@patch.object(DefaultSeqVarPredictor, "verify_pm4bp3")
def test_predict_pm4bp3_cdh1(mock_verify_pm4bp3, cdh1_predictor, seqvar, auto_acmg_data):
    """Test the predict_pm4bp3 method for CDH1."""
    # Set the mock return value for the verify_pm4bp3 method
    mock_pred = MagicMock()
    mock_pred.PM4 = True
    mock_pred.BP3 = False
    mock_pred.PM4_strength = AutoACMGStrength.PathogenicSupporting
    mock_pred.BP3_strength = AutoACMGStrength.BenignSupporting
    mock_verify_pm4bp3.return_value = (mock_pred, "PM4 is met")

    # Call the method under test
    pm4_result, bp3_result = cdh1_predictor.predict_pm4bp3(seqvar, auto_acmg_data)

    # Check PM4 result
    assert isinstance(
        pm4_result, AutoACMGCriteria
    ), "The PM4 result should be of type AutoACMGCriteria."
    assert pm4_result.prediction == AutoACMGPrediction.Met, "PM4 should be Met."
    assert (
        pm4_result.strength == AutoACMGStrength.PathogenicSupporting
    ), "PM4 strength should be PathogenicSupporting."
    assert "PM4 is met" in pm4_result.summary, "The summary should indicate PM4 is met."

    # Check BP3 result
    assert isinstance(
        bp3_result, AutoACMGCriteria
    ), "The BP3 result should be of type AutoACMGCriteria."
    assert bp3_result.prediction == AutoACMGPrediction.NotApplicable, "BP3 should be NotApplicable."
    assert (
        bp3_result.strength == AutoACMGStrength.BenignSupporting
    ), "BP3 strength should be BenignSupporting."
    assert (
        "BP3 is not applicable for CDH1" in bp3_result.summary
    ), "The summary should indicate BP3 is not applicable."


@patch.object(DefaultSeqVarPredictor, "verify_pm4bp3", return_value=(None, ""))
def test_predict_pm4bp3_fallback(mock_verify_pm4bp3, cdh1_predictor, seqvar, auto_acmg_data):
    """Test the fallback behavior for PM4 when verification fails."""
    # Call the method under test
    pm4_result, bp3_result = cdh1_predictor.predict_pm4bp3(seqvar, auto_acmg_data)

    # Check PM4 result
    assert (
        pm4_result.prediction == AutoACMGPrediction.Failed
    ), "PM4 should be Failed if verification fails."
    assert (
        "PM4 could not be verified." in pm4_result.summary
    ), "The summary should indicate PM4 could not be verified."

    # Check BP3 result
    assert bp3_result.prediction == AutoACMGPrediction.NotApplicable, "BP3 should be NotApplicable."
    assert (
        bp3_result.strength == AutoACMGStrength.BenignSupporting
    ), "BP3 strength should be BenignSupporting."
    assert (
        "BP3 is not applicable for CDH1" in bp3_result.summary
    ), "The summary should indicate BP3 is not applicable."


def test_predict_pp2bp1(cdh1_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for CDH1 predictor."""

    # Call the method under test
    pp2_result, bp1_result = cdh1_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


def test_predict_bp7_threshold_adjustment(cdh1_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted."""
    auto_acmg_data.thresholds.bp7_donor = 1  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = cdh1_predictor.predict_bp7(cdh1_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.bp7_donor == 7
    ), "The BP7 donor threshold should be adjusted to 7."
    assert (
        auto_acmg_data.thresholds.bp7_acceptor == 21
    ), "The BP7 acceptor threshold should be adjusted to 4."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(DefaultSeqVarPredictor, "predict_bp7")
def test_predict_bp7_fallback_to_default(mock_super_predict_bp7, cdh1_predictor, auto_acmg_data):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = cdh1_predictor.predict_bp7(cdh1_predictor.seqvar, auto_acmg_data)

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
