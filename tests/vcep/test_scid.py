from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.scid import SCIDPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="10", pos=100, delete="A", insert="T")


@pytest.fixture
def scid_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return SCIDPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_in_moderate_domain(scid_predictor, auto_acmg_data):
    """Test when the variant falls within a moderate domain for a SCID gene."""
    auto_acmg_data.hgnc_id = "HGNC:12765"  # FOXN1 gene
    auto_acmg_data.prot_pos = 300  # Within the DNA binding forkhead domain (270-367)
    result = scid_predictor.predict_pm1(scid_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for critical domain variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "critical residue in HGNC:12765" in result.summary
    ), "The summary should indicate the critical domain."


def test_predict_pm1_in_supporting_domain(scid_predictor, auto_acmg_data):
    """Test when the variant falls within a supporting domain for a SCID gene."""
    auto_acmg_data.hgnc_id = "HGNC:9831"  # RAG1 gene
    auto_acmg_data.prot_pos = 600  # Within the core domain (387-1011)
    result = scid_predictor.predict_pm1(scid_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for supporting domain variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "supporting region of HGNC:9831" in result.summary
    ), "The summary should indicate the supporting domain."


def test_predict_pm1_strong_residue(scid_predictor, auto_acmg_data):
    """Test when the variant affects a strong residue for a SCID gene."""
    auto_acmg_data.hgnc_id = "HGNC:6010"  # IL2RG gene
    auto_acmg_data.prot_pos = 62  # A conserved cysteine residue
    result = scid_predictor.predict_pm1(scid_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for strong residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicStrong
    ), "The strength should be PathogenicStrong."
    assert (
        "strong residue in HGNC:6010" in result.summary
    ), "The summary should indicate the strong residue."


def test_predict_pm1_not_applicable(scid_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for ADA, DCLRE1C, or IL7R."""
    auto_acmg_data.hgnc_id = "HGNC:186"  # ADA gene
    result = scid_predictor.predict_pm1(scid_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be applicable for ADA."
    assert "not applicable" in result.summary, "The summary should indicate non-applicability."


@patch("src.vcep.scid.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, scid_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the SCID VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default fallback for PM1.",
    )

    result = scid_predictor.predict_pm1(scid_predictor.seqvar, auto_acmg_data)
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default fallback for PM1." in result.summary
    ), "The summary should indicate the default fallback."
