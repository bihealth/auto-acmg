from unittest.mock import MagicMock, patch

import pytest

from src.defs.annonars_variant import VariantResult
from src.defs.auto_acmg import (
    PS1PM5,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
)
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor
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
    return AutoACMGSeqVarData()


@pytest.fixture
def variant_result():
    return VariantResult()


def test_is_pathogenic_with_pathogenic_variant(dicer1_predictor, variant_result):
    """Test _is_pathogenic method with a pathogenic variant."""
    # Create a mock ClinvarRecord with a pathogenic classification
    clinvar_record = MagicMock()
    clinvar_record.classifications.germlineClassification.description = "Pathogenic"
    variant_result.clinvar = MagicMock(records=[clinvar_record])

    result = dicer1_predictor._is_pathogenic(variant_result)
    assert result is True, "Should return True for a pathogenic variant"


def test_is_pathogenic_with_non_pathogenic_variant(dicer1_predictor, variant_result):
    """Test _is_pathogenic method with a non-pathogenic variant."""
    # Create a mock ClinvarRecord with a non-pathogenic classification
    clinvar_record = MagicMock()
    clinvar_record.classifications.germlineClassification.description = "Benign"
    variant_result.clinvar = MagicMock(records=[clinvar_record])

    result = dicer1_predictor._is_pathogenic(variant_result)
    assert result is False, "Should return False for a non-pathogenic variant"


def test_is_pathogenic_with_no_clinvar_data(dicer1_predictor, variant_result):
    """Test _is_pathogenic method with no ClinVar data."""
    variant_result.clinvar = None

    result = dicer1_predictor._is_pathogenic(variant_result)
    assert result is False, "Should return False when no ClinVar data is available"


def test_is_pathogenic_with_empty_clinvar_records(dicer1_predictor, variant_result):
    """Test _is_pathogenic method with empty ClinVar records."""
    variant_result.clinvar = MagicMock(records=[])

    result = dicer1_predictor._is_pathogenic(variant_result)
    assert result is False, "Should return False when ClinVar records are empty"


def test_is_pathogenic_with_multiple_classifications(dicer1_predictor, variant_result):
    """Test _is_pathogenic method with multiple ClinVar classifications."""
    # Create mock ClinvarRecords with different classifications
    clinvar_record1 = MagicMock()
    clinvar_record1.classifications.germlineClassification.description = "Benign"
    clinvar_record2 = MagicMock()
    clinvar_record2.classifications.germlineClassification.description = "Pathogenic"
    variant_result.clinvar = MagicMock(records=[clinvar_record1, clinvar_record2])

    result = dicer1_predictor._is_pathogenic(variant_result)
    assert result is True, "Should return True if any classification is Pathogenic"


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_overrides(mock_super_verify, dicer1_predictor, seqvar, auto_acmg_data):
    """Test that the overridden verify_ps1pm5 method in DICER1Predictor works correctly."""
    # Set up the mock to return PS1 and PM5 as applicable initially
    mock_super_verify.return_value = (PS1PM5(PS1=True, PM5=True), "Initial evaluation")

    # Setup the data
    auto_acmg_data.consequence = MagicMock(mehari=["missense_variant"])
    auto_acmg_data.thresholds = MagicMock(
        spliceAI_acceptor_gain=0.5,
        spliceAI_acceptor_loss=0.5,
        spliceAI_donor_gain=0.5,
        spliceAI_donor_loss=0.5,
    )
    auto_acmg_data.scores = MagicMock(
        cadd=MagicMock(spliceAI_acceptor_gain=0.6, spliceAI_donor_gain=0.6)
    )

    # Run the method under test
    prediction, comment = dicer1_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that the splicing effect leads to overriding PS1 and PM5 as not applicable
    assert not prediction.PS1, "PS1 should be marked as not applicable due to splicing effect."
    assert not prediction.PM5, "PM5 should be marked as not applicable due to splicing effect."
    assert "Variant affects splicing" in comment, "Comment should note the splicing effect."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)


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
def test_predict_pm2ba1bs1bs2(mock_super_method, dicer1_predictor, auto_acmg_data, seqvar):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = dicer1_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.000005
    assert auto_acmg_data.thresholds.ba1_benign == 0.003
    assert auto_acmg_data.thresholds.bs1_benign == 0.0003

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_predict_pp2bp1(dicer1_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for DICER1 predictor."""

    # Call the method under test
    pp2_result, bp1_result = dicer1_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


def test_bp3_not_applicable(dicer1_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = dicer1_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


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


@patch.object(DefaultSeqVarPredictor, "predict_bp7")
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
