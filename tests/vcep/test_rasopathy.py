from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.rasopathy import RASopathyPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="12", pos=100, delete="A", insert="T")


@pytest.fixture
def rasopathy_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return RASopathyPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_in_critical_region(rasopathy_predictor, auto_acmg_data):
    """Test when the variant falls within a critical region for a RASopathy gene."""
    auto_acmg_data.hgnc_id = "HGNC:7989"  # NRAS gene
    auto_acmg_data.prot_pos = 15  # Within the critical P-loop region (10-17)
    result = rasopathy_predictor.predict_pm1(rasopathy_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for critical region variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "critical residue in HGNC:7989" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_in_critical_exon(rasopathy_predictor, auto_acmg_data):
    """Test when the variant falls within a critical exon for a RASopathy gene."""
    auto_acmg_data.hgnc_id = "HGNC:1097"  # BRAF gene
    auto_acmg_data.prot_pos = 480  # Position in the P-loop region, but within exon 11
    rasopathy_predictor._get_affected_exon = MagicMock(return_value=11)  # Mock the affected exon
    result = rasopathy_predictor.predict_pm1(rasopathy_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for critical exon variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert "critical exon 11" in result.summary, "The summary should indicate the critical exon."


def test_predict_pm1_not_applicable(rasopathy_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for certain RASopathy genes."""
    auto_acmg_data.hgnc_id = "HGNC:15454"  # SHOC2 gene
    result = rasopathy_predictor.predict_pm1(rasopathy_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be applicable for SHOC2."
    assert "not applicable" in result.summary, "The summary should indicate non-applicability."


def test_predict_pm1_outside_critical_region(rasopathy_predictor, auto_acmg_data):
    """Test when the variant does not fall within any critical region for RASopathy genes."""
    auto_acmg_data.hgnc_id = "HGNC:7989"  # NRAS gene
    auto_acmg_data.prot_pos = 200  # Position outside all critical regions
    rasopathy_predictor._get_affected_exon = MagicMock(return_value=1)  # Mock the affected exon
    result = rasopathy_predictor.predict_pm1(rasopathy_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for non-critical region variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not meet the PM1 criteria" in result.summary
    ), "The summary should indicate no critical region."


@patch("src.vcep.rasopathy.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, rasopathy_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the RASopathy VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = rasopathy_predictor.predict_pm1(rasopathy_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."
