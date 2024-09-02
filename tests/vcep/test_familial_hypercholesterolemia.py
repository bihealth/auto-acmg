from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import FamilialHypercholesterolemiaPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="19", pos=100, delete="A", insert="T")


@pytest.fixture
def familial_hypercholesterolemia_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return FamilialHypercholesterolemiaPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_in_critical_residue(familial_hypercholesterolemia_predictor, auto_acmg_data):
    """Test when variant affects a critical residue in LDLR."""
    auto_acmg_data.hgnc_id = "HGNC:6547"  # LDLR gene
    auto_acmg_data.prot_pos = 75  # Critical residue in LDLR
    result = familial_hypercholesterolemia_predictor.predict_pm1(
        familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for a variant affecting a critical residue."
    assert (
        "affects a critical cysteine residue" in result.summary
    ), "The summary should indicate the critical residue."


def test_predict_pm1_in_exon_4(familial_hypercholesterolemia_predictor, auto_acmg_data):
    """Test when variant affects the 4th exon in LDLR."""
    auto_acmg_data.hgnc_id = "HGNC:6547"  # LDLR gene
    with patch.object(FamilialHypercholesterolemiaPredictor, "_get_affected_exon", return_value=4):
        result = familial_hypercholesterolemia_predictor.predict_pm1(
            familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
        )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for a variant affecting the 4th exon."
    assert (
        "affects the 4th exon" in result.summary
    ), "The summary should indicate the affected exon."


def test_predict_pm1_outside_critical_region(
    familial_hypercholesterolemia_predictor, auto_acmg_data
):
    """Test when variant does not fall within any critical region or exon 4 for LDLR."""
    auto_acmg_data.hgnc_id = "HGNC:6547"  # LDLR gene
    auto_acmg_data.prot_pos = 500  # Position outside critical residues and not in exon 4
    with patch.object(FamilialHypercholesterolemiaPredictor, "_get_affected_exon", return_value=5):
        result = familial_hypercholesterolemia_predictor.predict_pm1(
            familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
        )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for a variant outside critical residues and exon 4."
    assert (
        "Variant does not meet the PM1 criteria for LDLR" in result.summary
    ), "The summary should indicate that no criteria were met."


@patch("src.vcep.familial_hypercholesterolemia.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, familial_hypercholesterolemia_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method for genes not in PM1_CLUSTER_LDLR."""
    auto_acmg_data.hgnc_id = "HGNC:99999"  # Gene not in PM1_CLUSTER_LDLR
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = familial_hypercholesterolemia_predictor.predict_pm1(
        familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met if the gene is not in the PM1_CLUSTER_LDLR."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should be from the default fallback."


def test_predict_pm1_edge_case_start_boundary(
    familial_hypercholesterolemia_predictor, auto_acmg_data
):
    """Test when variant falls exactly on the start boundary of a critical residue."""
    auto_acmg_data.hgnc_id = "HGNC:6547"  # LDLR gene
    auto_acmg_data.prot_pos = 27  # Start boundary of the critical residues (27)
    result = familial_hypercholesterolemia_predictor.predict_pm1(
        familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the start boundary of a critical residue."
    assert (
        "affects a critical cysteine residue" in result.summary
    ), "The summary should indicate the critical residue."


def test_predict_pm1_edge_case_end_boundary(
    familial_hypercholesterolemia_predictor, auto_acmg_data
):
    """Test when variant falls exactly on the end boundary of a critical residue."""
    auto_acmg_data.hgnc_id = "HGNC:6547"  # LDLR gene
    auto_acmg_data.prot_pos = 711  # End boundary of the critical residues (711)
    result = familial_hypercholesterolemia_predictor.predict_pm1(
        familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the end boundary of a critical residue."
    assert (
        "affects a critical cysteine residue" in result.summary
    ), "The summary should indicate the critical residue."


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
def test_predict_pm2ba1bs1bs2(
    mock_super_method, familial_hypercholesterolemia_predictor, auto_acmg_data, seqvar
):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = familial_hypercholesterolemia_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.0002
    assert auto_acmg_data.thresholds.ba1_benign == 0.005
    assert auto_acmg_data.thresholds.bs1_benign == 0.002

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_bp3_not_applicable(familial_hypercholesterolemia_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = familial_hypercholesterolemia_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


def test_predict_pp2bp1(familial_hypercholesterolemia_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Familial Hypercholesterolemia predictor."""

    # Call the method under test
    pp2_result, bp1_result = familial_hypercholesterolemia_predictor.predict_pp2bp1(
        seqvar, auto_acmg_data
    )

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
