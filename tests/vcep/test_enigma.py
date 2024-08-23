from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import ENIGMAPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="17", pos=100, delete="A", insert="T")


@pytest.fixture
def enigma_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return ENIGMAPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_not_applicable_brca1(enigma_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for BRCA1."""
    auto_acmg_data.hgnc_id = "HGNC:1100"  # BRCA1 gene
    result = enigma_predictor.predict_pm1(enigma_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for BRCA1."
    assert (
        result.summary == "PM1 is not applicable for BRCA1 and BRCA2."
    ), "The summary should indicate that PM1 is not applicable for BRCA1."


def test_predict_pm1_not_applicable_brca2(enigma_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for BRCA2."""
    auto_acmg_data.hgnc_id = "HGNC:1101"  # BRCA2 gene
    result = enigma_predictor.predict_pm1(enigma_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for BRCA2."
    assert (
        result.summary == "PM1 is not applicable for BRCA1 and BRCA2."
    ), "The summary should indicate that PM1 is not applicable for BRCA2."


@patch("src.vcep.enigma.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, enigma_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for genes other than BRCA1 and BRCA2."""
    auto_acmg_data.hgnc_id = "HGNC:99999"  # Gene not BRCA1 or BRCA2
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = enigma_predictor.predict_pm1(enigma_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for genes other than BRCA1 and BRCA2."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should be from the default fallback."


def test_predict_pm1_name(enigma_predictor, auto_acmg_data):
    """Test the name of the criteria returned by the ENIGMA predictor."""
    auto_acmg_data.hgnc_id = "HGNC:1100"  # BRCA1 gene
    result = enigma_predictor.predict_pm1(enigma_predictor.seqvar, auto_acmg_data)

    assert result.name == "PM1", "The name of the criteria should be 'PM1'."


def test_predict_pm1_strength(enigma_predictor, auto_acmg_data):
    """Test the strength level returned by the ENIGMA predictor."""
    auto_acmg_data.hgnc_id = "HGNC:1101"  # BRCA2 gene
    result = enigma_predictor.predict_pm1(enigma_predictor.seqvar, auto_acmg_data)

    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate for ENIGMA."
