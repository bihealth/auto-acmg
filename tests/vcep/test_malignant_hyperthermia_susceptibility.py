from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.malignant_hyperthermia_susceptibility import MalignantHyperthermiaPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="19", pos=100, delete="A", insert="T")


@pytest.fixture
def malignant_hyperthermia_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return MalignantHyperthermiaPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_met_for_moderate_region(malignant_hyperthermia_predictor, auto_acmg_data):
    """Test when the variant falls within a moderate region for RYR1."""
    auto_acmg_data.prot_pos = 300  # Within the moderate region (1-552) for RYR1
    auto_acmg_data.hgnc_id = "HGNC:10483"  # RYR1 gene
    result = malignant_hyperthermia_predictor.predict_pm1(
        malignant_hyperthermia_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for variants in moderate regions."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant falls within a critical region in RYR1 between positions 1-552" in result.summary
    ), "The summary should indicate the affected region."


def test_predict_pm1_met_for_supporting_region(malignant_hyperthermia_predictor, auto_acmg_data):
    """Test when the variant falls within a supporting region for RYR1."""
    auto_acmg_data.prot_pos = 4800  # Within the supporting region (4631-4991) for RYR1
    auto_acmg_data.hgnc_id = "HGNC:10483"  # RYR1 gene
    result = malignant_hyperthermia_predictor.predict_pm1(
        malignant_hyperthermia_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for variants in supporting regions."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "Variant falls within a critical region in RYR1 between positions 4631-4991"
        in result.summary
    ), "The summary should indicate the affected region."


def test_predict_pm1_not_met(malignant_hyperthermia_predictor, auto_acmg_data):
    """Test when the variant does not fall within any critical regions for RYR1."""
    auto_acmg_data.prot_pos = 5000  # Not within any critical regions
    auto_acmg_data.hgnc_id = "HGNC:10483"  # RYR1 gene
    result = malignant_hyperthermia_predictor.predict_pm1(
        malignant_hyperthermia_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for non-critical region variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant does not meet the PM1 criteria for Malignant Hyperthermia Susceptibility"
        in result.summary
    ), "The summary should indicate the lack of criteria met."


@patch("src.vcep.malignant_hyperthermia_susceptibility.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, malignant_hyperthermia_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method if logic changes."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the specific logic
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = malignant_hyperthermia_predictor.predict_pm1(
        malignant_hyperthermia_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."
