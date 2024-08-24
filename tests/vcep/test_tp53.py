from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.tp53 import TP53Predictor


@pytest.fixture
def seqvar():
    return SeqVar(
        genome_release=GenomeRelease.GRCh37, chrom="17", pos=7578406, delete="C", insert="T"
    )


@pytest.fixture
def tp53_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return TP53Predictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_in_moderate_residue(tp53_predictor, auto_acmg_data):
    """Test when the variant affects a critical residue in TP53."""
    auto_acmg_data.hgnc_id = "HGNC:11998"  # TP53 gene
    auto_acmg_data.prot_pos = 175  # A critical residue in TP53
    result = tp53_predictor.predict_pm1(tp53_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "affects a critical residue in TP53" in result.summary
    ), "The summary should indicate the PM1 criteria."


def test_predict_pm1_not_met(tp53_predictor, auto_acmg_data):
    """Test when the variant does not meet the PM1 criteria for TP53."""
    auto_acmg_data.hgnc_id = "HGNC:11998"  # TP53 gene
    auto_acmg_data.prot_pos = 300  # Position outside of critical residues
    result = tp53_predictor.predict_pm1(tp53_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for non-critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not meet the PM1 criteria for TP53" in result.summary
    ), "The summary should indicate the PM1 criteria."


@patch("src.vcep.tp53.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, tp53_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the TP53 VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default fallback for PM1.",
    )
    result = tp53_predictor.predict_pm1(tp53_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
