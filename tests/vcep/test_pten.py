from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.pten import PTENPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="10", pos=100, delete="A", insert="T")


@pytest.fixture
def pten_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return PTENPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_in_catalytic_motifs(pten_predictor, auto_acmg_data):
    """Test when the variant falls within the catalytic motifs of PTEN."""
    auto_acmg_data.hgnc_id = "HGNC:9588"  # PTEN gene
    auto_acmg_data.prot_pos = 92  # Within the catalytic motif (90-94)
    result = pten_predictor.predict_pm1(pten_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for a variant in the catalytic motif."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "catalytic motifs of PTEN" in result.summary
    ), "The summary should indicate the catalytic motif."


def test_predict_pm1_outside_catalytic_motifs(pten_predictor, auto_acmg_data):
    """Test when the variant does not fall within the catalytic motifs of PTEN."""
    auto_acmg_data.hgnc_id = "HGNC:9588"  # PTEN gene
    auto_acmg_data.prot_pos = 150  # Outside the catalytic motifs
    result = pten_predictor.predict_pm1(pten_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for a variant outside the catalytic motifs."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant does not meet the PM1 criteria for PTEN" in result.summary
    ), "The summary should indicate no criteria were met."


@patch("src.vcep.pten.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, pten_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the PTEN VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = pten_predictor.predict_pm1(pten_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."
