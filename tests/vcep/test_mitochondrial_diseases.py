from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
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
    mock_super_method, mitochondrial_diseases_predictor, auto_acmg_data, seqvar
):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = mitochondrial_diseases_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.00005
    assert auto_acmg_data.thresholds.ba1_benign == 0.001
    assert auto_acmg_data.thresholds.bs1_benign == 0.0005

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_predict_pp2bp1(mitochondrial_diseases_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Mitochondrial Diseases predictor."""

    # Call the method under test
    pp2_result, bp1_result = mitochondrial_diseases_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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
