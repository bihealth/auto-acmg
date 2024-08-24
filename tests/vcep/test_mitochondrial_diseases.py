from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.mitochondrial_diseases import MitochondrialDiseasesPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="X", pos=100, delete="A", insert="T")


@pytest.fixture
def mitochondrial_diseases_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return MitochondrialDiseasesPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_not_applicable_ethe1(mitochondrial_diseases_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for ETHE1."""
    auto_acmg_data.hgnc_id = "HGNC:23287"  # ETHE1 gene
    result = mitochondrial_diseases_predictor.predict_pm1(
        mitochondrial_diseases_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be applicable for ETHE1."
    assert (
        "PM1 is not applicable for HGNC:23287" in result.summary
    ), "The summary should indicate non-applicability."


def test_predict_pm1_not_applicable_pdah1(mitochondrial_diseases_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for PDHA1."""
    auto_acmg_data.hgnc_id = "HGNC:8806"  # PDHA1 gene
    result = mitochondrial_diseases_predictor.predict_pm1(
        mitochondrial_diseases_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be applicable for PDHA1."
    assert (
        "PM1 is not applicable for HGNC:8806" in result.summary
    ), "The summary should indicate non-applicability."


def test_predict_pm1_not_applicable_polg(mitochondrial_diseases_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for POLG."""
    auto_acmg_data.hgnc_id = "HGNC:9179"  # POLG gene
    result = mitochondrial_diseases_predictor.predict_pm1(
        mitochondrial_diseases_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be applicable for POLG."
    assert (
        "PM1 is not applicable for HGNC:9179" in result.summary
    ), "The summary should indicate non-applicability."


def test_predict_pm1_not_applicable_slc19a3(mitochondrial_diseases_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for SLC19A3."""
    auto_acmg_data.hgnc_id = "HGNC:16266"  # SLC19A3 gene
    result = mitochondrial_diseases_predictor.predict_pm1(
        mitochondrial_diseases_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be applicable for SLC19A3."
    assert (
        "PM1 is not applicable for HGNC:16266" in result.summary
    ), "The summary should indicate non-applicability."


@patch("src.vcep.mitochondrial_diseases.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, mitochondrial_diseases_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method if logic changes."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the specific logic
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = mitochondrial_diseases_predictor.predict_pm1(
        mitochondrial_diseases_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."
