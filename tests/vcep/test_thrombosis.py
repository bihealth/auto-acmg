from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.thrombosis import ThrombosisPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def thrombosis_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return ThrombosisPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_in_cysteine_residue(thrombosis_predictor, auto_acmg_data):
    """Test when the variant affects a critical cysteine residue in SERPINC1."""
    auto_acmg_data.hgnc_id = "HGNC:775"  # SERPINC1 gene
    auto_acmg_data.prot_pos = 40  # A cysteine residue involved in disulfide bridges
    result = thrombosis_predictor.predict_pm1(thrombosis_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for cysteine residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "meets the PM1 criteria for SERPINC1" in result.summary
    ), "The summary should indicate the PM1 criteria."


def test_predict_pm1_in_heparin_binding_residue(thrombosis_predictor, auto_acmg_data):
    """Test when the variant affects a residue in the heparin binding site of SERPINC1."""
    auto_acmg_data.hgnc_id = "HGNC:775"  # SERPINC1 gene
    auto_acmg_data.prot_pos = 39  # A heparin binding site residue
    result = thrombosis_predictor.predict_pm1(thrombosis_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for heparin binding site variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "meets the PM1 criteria for SERPINC1" in result.summary
    ), "The summary should indicate the PM1 criteria."


def test_predict_pm1_in_reactive_site_residue(thrombosis_predictor, auto_acmg_data):
    """Test when the variant affects a residue in the reactive site of SERPINC1."""
    auto_acmg_data.hgnc_id = "HGNC:775"  # SERPINC1 gene
    auto_acmg_data.prot_pos = 414  # A reactive site residue
    result = thrombosis_predictor.predict_pm1(thrombosis_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for reactive site variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "meets the PM1 criteria for SERPINC1" in result.summary
    ), "The summary should indicate the PM1 criteria."


def test_predict_pm1_not_met(thrombosis_predictor, auto_acmg_data):
    """Test when the variant does not meet the PM1 criteria for SERPINC1."""
    auto_acmg_data.hgnc_id = "HGNC:775"  # SERPINC1 gene
    auto_acmg_data.prot_pos = 500  # Position outside of critical regions
    result = thrombosis_predictor.predict_pm1(thrombosis_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for non-critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not meet the PM1 criteria for SERPINC1" in result.summary
    ), "The summary should indicate the PM1 criteria."


@patch("src.vcep.thrombosis.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, thrombosis_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the Thrombosis VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default fallback for PM1.",
    )

    result = thrombosis_predictor.predict_pm1(thrombosis_predictor.seqvar, auto_acmg_data)
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default fallback for PM1." in result.summary
    ), "The summary should indicate the default fallback."