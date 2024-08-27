from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import DICER1Predictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="14", pos=100, delete="A", insert="T")


@pytest.fixture
def dicer1_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return DICER1Predictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_moderate_criteria_residue(dicer1_predictor, auto_acmg_data):
    """Test when the variant affects a metal ion-binding residue in DICER1."""
    auto_acmg_data.prot_pos = 1705  # Metal ion-binding residue
    auto_acmg_data.hgnc_id = "HGNC:17098"  # DICER1 gene
    result = dicer1_predictor.predict_pm1(dicer1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for metal ion-binding residues."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "affects a metal ion-binding residue" in result.summary
    ), "The summary should indicate the metal ion-binding residue."


def test_predict_pm1_supporting_criteria_domain(dicer1_predictor, auto_acmg_data):
    """Test when the variant affects a residue within the RNase IIIb domain but outside metal ion-binding residues."""
    auto_acmg_data.prot_pos = (
        1720  # Within the RNase IIIb domain but not a metal ion-binding residue
    )
    auto_acmg_data.hgnc_id = "HGNC:17098"  # DICER1 gene
    result = dicer1_predictor.predict_pm1(dicer1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for residues within the RNase IIIb domain."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "affects a residue in the RNase IIIb domain" in result.summary
    ), "The summary should indicate the RNase IIIb domain."


def test_predict_pm1_outside_critical_region(dicer1_predictor, auto_acmg_data):
    """Test when the variant does not affect any critical residues or domains in DICER1."""
    auto_acmg_data.prot_pos = 1600  # Outside any critical region or domain
    auto_acmg_data.hgnc_id = "HGNC:17098"  # DICER1 gene
    result = dicer1_predictor.predict_pm1(dicer1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for non-critical regions."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not affect a critical domain" in result.summary
    ), "The summary should indicate no critical domain was affected."


def test_predict_pm1_not_applicable_gene(dicer1_predictor, auto_acmg_data):
    """Test when PM1 is not applicable because the gene is not DICER1."""
    auto_acmg_data.hgnc_id = "HGNC:99999"  # Gene not in PM1_CLUSTER mapping
    result = dicer1_predictor.predict_pm1(dicer1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met if gene is not in PM1_CLUSTER."
    assert (
        "does not affect a critical domain for DICER1" in result.summary
    ), "The summary should indicate gene is not applicable."


def test_predict_pm1_edge_case_start_boundary_moderate(dicer1_predictor, auto_acmg_data):
    """Test when the variant falls exactly on the start boundary of a moderate criteria residue."""
    auto_acmg_data.prot_pos = 1344  # Start boundary of a metal ion-binding residue
    auto_acmg_data.hgnc_id = "HGNC:17098"  # DICER1 gene
    result = dicer1_predictor.predict_pm1(dicer1_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for the start boundary of a metal ion-binding residue."
    assert (
        "affects a metal ion-binding residue" in result.summary
    ), "The summary should indicate the metal ion-binding residue."


def test_predict_pm1_edge_case_end_boundary_supporting(dicer1_predictor, auto_acmg_data):
    """Test when the variant falls exactly on the end boundary of a supporting criteria domain."""
    auto_acmg_data.prot_pos = 1846  # End boundary of the RNase IIIb domain
    auto_acmg_data.hgnc_id = "HGNC:17098"  # DICER1 gene
    result = dicer1_predictor.predict_pm1(dicer1_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for the end boundary of the RNase IIIb domain."
    assert (
        "affects a residue in the RNase IIIb domain" in result.summary
    ), "The summary should indicate the RNase IIIb domain."


def test_predict_bp7_threshold_adjustment(dicer1_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted."""
    auto_acmg_data.thresholds.bp7_donor = 1  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = dicer1_predictor.predict_bp7(dicer1_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.bp7_donor == 7
    ), "The BP7 donor threshold should be adjusted to 7."
    assert (
        auto_acmg_data.thresholds.bp7_acceptor == 21
    ), "The BP7 acceptor threshold should be adjusted to 4."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(DefaultPredictor, "predict_bp7")
def test_predict_bp7_fallback_to_default(mock_super_predict_bp7, dicer1_predictor, auto_acmg_data):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = dicer1_predictor.predict_bp7(dicer1_predictor.seqvar, auto_acmg_data)

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
