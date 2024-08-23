from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import FBN1Predictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="15", pos=100, delete="A", insert="T")


@pytest.fixture
def fbn1_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return FBN1Predictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_strong_residue(fbn1_predictor, auto_acmg_data):
    """Test when variant affects a strong critical residue in FBN1."""
    auto_acmg_data.hgnc_id = "HGNC:3603"  # FBN1 gene
    auto_acmg_data.prot_pos = 250  # Strong residue in FBN1
    result = fbn1_predictor.predict_pm1(fbn1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for a variant affecting a strong critical residue."
    assert (
        result.strength == AutoACMGStrength.PathogenicStrong
    ), "The strength should be PathogenicStrong for strong residues."
    assert (
        "affects a critical cysteine residue" in result.summary
    ), "The summary should indicate the strong critical residue."


def test_predict_pm1_moderate_residue(fbn1_predictor, auto_acmg_data):
    """Test when variant affects a moderate critical residue in FBN1."""
    auto_acmg_data.hgnc_id = "HGNC:3603"  # FBN1 gene
    auto_acmg_data.prot_pos = 100  # Moderate residue in FBN1
    result = fbn1_predictor.predict_pm1(fbn1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for a variant affecting a moderate critical residue."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate for moderate residues."
    assert (
        "Variant affects a residue in FBN1" in result.summary
    ), "The summary should indicate the moderate critical residue."


def test_predict_pm1_outside_critical_residues(fbn1_predictor, auto_acmg_data):
    """Test when variant does not fall within any critical residues in FBN1."""
    auto_acmg_data.hgnc_id = "HGNC:3603"  # FBN1 gene
    auto_acmg_data.prot_pos = 5000  # Position outside critical residues
    result = fbn1_predictor.predict_pm1(fbn1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for a variant outside critical residues."
    assert (
        "Variant does not meet the PM1 criteria for FBN1" in result.summary
    ), "The summary should indicate no criteria were met."


@patch("src.vcep.fbn1.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, fbn1_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for genes not in PM1_CLUSTER."""
    auto_acmg_data.hgnc_id = "HGNC:99999"  # Gene not in PM1_CLUSTER
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = fbn1_predictor.predict_pm1(fbn1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met if the gene is not in the PM1_CLUSTER."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should be from the default fallback."


def test_predict_pm1_edge_case_strong_boundary(fbn1_predictor, auto_acmg_data):
    """Test when variant falls exactly on the boundary of a strong critical residue."""
    auto_acmg_data.hgnc_id = "HGNC:3603"  # FBN1 gene
    auto_acmg_data.prot_pos = 2686  # Boundary of a strong critical residue in FBN1
    result = fbn1_predictor.predict_pm1(fbn1_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the boundary of a strong critical residue."
    assert (
        result.strength == AutoACMGStrength.PathogenicStrong
    ), "The strength should be PathogenicStrong for strong residues."
    assert (
        "affects a critical cysteine residue" in result.summary
    ), "The summary should indicate the strong critical residue."


def test_predict_pm1_edge_case_moderate_boundary(fbn1_predictor, auto_acmg_data):
    """Test when variant falls exactly on the boundary of a moderate critical residue."""
    auto_acmg_data.hgnc_id = "HGNC:3603"  # FBN1 gene
    auto_acmg_data.prot_pos = 2390  # Boundary of a moderate critical residue in FBN1
    result = fbn1_predictor.predict_pm1(fbn1_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the boundary of a moderate critical residue."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate for moderate residues."
    assert (
        "Variant affects a residue in FBN1" in result.summary
    ), "The summary should indicate the moderate critical residue."
