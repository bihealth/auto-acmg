from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import HBOPCPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def hbopc_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return HBOPCPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_not_applicable(hbopc_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for ATM and PALB2 in Hereditary Breast, Ovarian, and Pancreatic Cancer."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for ATM."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting for ATM."
    assert (
        "PM1 is not applicable for HGNC:795" in result.summary
    ), "The summary should indicate that PM1 is not applicable for ATM."


def test_predict_pm1_not_applicable_palb2(hbopc_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for PALB2 in Hereditary Breast, Ovarian, and Pancreatic Cancer."""
    auto_acmg_data.hgnc_id = "HGNC:26144"  # PALB2 gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for PALB2."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting for PALB2."
    assert (
        "PM1 is not applicable for HGNC:26144" in result.summary
    ), "The summary should indicate that PM1 is not applicable for PALB2."


@patch("src.vcep.hbopc.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, hbopc_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method if logic changes."""
    auto_acmg_data.hgnc_id = "HGNC:111111"  # Gene not in the specific logic
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met if the gene is not ATM or PALB2."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


def test_predict_pm1_strength_level(hbopc_predictor, auto_acmg_data):
    """Test the strength level for PM1 when not applicable."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting when PM1 is not applicable."


def test_predict_pm1_name(hbopc_predictor, auto_acmg_data):
    """Test the name of the criteria returned by the HBOPC predictor."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert result.name == "PM1", "The name of the criteria should be 'PM1'."
