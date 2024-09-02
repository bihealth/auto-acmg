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


def test_bs2_not_applicable(myeloid_malignancy_predictor, auto_acmg_data):
    """Test BS2 is not applicable for Myeloid Malignancy."""
    result = myeloid_malignancy_predictor._bs2_not_applicable(auto_acmg_data)
    assert result is True, "BS2 should always be not applicable"


@patch.object(MyeloidMalignancyPredictor, "_get_af", return_value=0.1)
@patch.object(MyeloidMalignancyPredictor, "_ba1_exception", return_value=False)
def test_verify_pm2ba1bs1bs2(
    mock_get_af, mock_ba1_exception, myeloid_malignancy_predictor, auto_acmg_data, seqvar
):
    # Setup: Adjusting the thresholds to test under different conditions
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.02
    auto_acmg_data.thresholds.pm2_pathogenic = 0.001

    # Call the method under test
    result, comment = myeloid_malignancy_predictor.verify_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assertions to validate the expected behavior
    assert result.BA1 is True, "Expected PM2 to be True based on the mocked allele frequency"

    # Assert changed thresholds
    assert (
        auto_acmg_data.thresholds.ba1_benign == 0.0015
    ), "BA1 threshold should be adjusted to 0.0015"
    assert (
        auto_acmg_data.thresholds.bs1_benign == 0.00015
    ), "BS1 threshold should be adjusted to 0.00015"


def test_bp3_not_applicable(myeloid_malignancy_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = myeloid_malignancy_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


def test_predict_pp2bp1(myeloid_malignancy_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Myeloid Malignancy."""

    # Call the method under test
    pp2_result, bp1_result = myeloid_malignancy_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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
