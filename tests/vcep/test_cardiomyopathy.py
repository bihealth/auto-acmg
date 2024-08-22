from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import CardiomyopathyPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def cardiomyopathy_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object, which will be passed along
    return CardiomyopathyPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_in_critical_domain(cardiomyopathy_predictor, auto_acmg_data):
    """Test when the variant falls within a critical domain for MYH7."""
    auto_acmg_data.prot_pos = 500  # Within the critical domain (167-931) for MYH7
    auto_acmg_data.hgnc_id = "HGNC:7577"  # MYH7 gene
    result = cardiomyopathy_predictor.predict_pm1(cardiomyopathy_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for a variant in a critical domain."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."


def test_predict_pm1_outside_critical_domain(cardiomyopathy_predictor, auto_acmg_data):
    """Test when the variant does not fall within any critical domain."""
    auto_acmg_data.prot_pos = 1000  # Outside the critical domain (167-931) for MYH7
    auto_acmg_data.hgnc_id = "HGNC:7577"  # MYH7 gene
    result = cardiomyopathy_predictor.predict_pm1(cardiomyopathy_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for a variant outside any critical domain."
    assert (
        "does not fall within any critical domain" in result.summary
    ), "The summary should indicate the lack of critical domain."


def test_predict_pm1_not_applicable(cardiomyopathy_predictor, auto_acmg_data):
    """Test when the PM1 criterion is not applicable."""
    auto_acmg_data.hgnc_id = "HGNC:12010"  # Example for TPM1, where PM1 is not applicable
    result = cardiomyopathy_predictor.predict_pm1(cardiomyopathy_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for TPM1."
    assert "Not applicable" in result.summary, "The summary should indicate non-applicability."


@patch("src.vcep.cardiomyopathy.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, cardiomyopathy_predictor, auto_acmg_data
):
    """Test when the transcript ID is not in the PM1_CLUSTER mapping, it should fallback to default PM1 prediction."""
    auto_acmg_data.hgnc_id = "HGNC:111111111111111"  # Not in the PM1_CLUSTER mapping
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = cardiomyopathy_predictor.predict_pm1(cardiomyopathy_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        result.summary == "Default PM1 prediction fallback."
    ), "The summary should be from the default fallback."


def test_predict_pm1_edge_case_start_boundary(cardiomyopathy_predictor, auto_acmg_data):
    """Test when variant falls exactly on the start boundary of a critical domain."""
    auto_acmg_data.prot_pos = 167  # Start boundary of the MYH7 critical domain (167-931)
    auto_acmg_data.hgnc_id = "HGNC:7577"  # MYH7 gene
    result = cardiomyopathy_predictor.predict_pm1(cardiomyopathy_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the start boundary of a critical domain."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."


def test_predict_pm1_edge_case_end_boundary(cardiomyopathy_predictor, auto_acmg_data):
    """Test when variant falls exactly on the end boundary of a critical domain."""
    auto_acmg_data.prot_pos = 931  # End boundary of the MYH7 critical domain (167-931)
    auto_acmg_data.hgnc_id = "HGNC:7577"  # MYH7 gene
    result = cardiomyopathy_predictor.predict_pm1(cardiomyopathy_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the end boundary of a critical domain."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."
