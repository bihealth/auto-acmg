from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.leber_congenital_amaurosis import LeberCongenitalAmaurosisPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def leber_congenital_amaurosis_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return LeberCongenitalAmaurosisPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_met_for_critical_residue(leber_congenital_amaurosis_predictor, auto_acmg_data):
    """Test when the variant falls within the critical residues for RPE65."""
    auto_acmg_data.prot_pos = 180  # Critical residue for RPE65
    auto_acmg_data.hgnc_id = "HGNC:10294"  # RPE65 gene
    result = leber_congenital_amaurosis_predictor.predict_pm1(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant affects a residue in RPE65 at position 180" in result.summary
    ), "The summary should indicate the affected residue."


def test_predict_pm1_met_for_residues_range(leber_congenital_amaurosis_predictor, auto_acmg_data):
    """Test when the variant falls within the residue range for RPE65."""
    auto_acmg_data.prot_pos = 110  # Within the residue range 107-125 for RPE65
    auto_acmg_data.hgnc_id = "HGNC:10294"  # RPE65 gene
    result = leber_congenital_amaurosis_predictor.predict_pm1(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for residues within the specified range."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant affects a residue in RPE65 at position 110" in result.summary
    ), "The summary should indicate the affected residue."


def test_predict_pm1_not_met(leber_congenital_amaurosis_predictor, auto_acmg_data):
    """Test when the variant does not fall within any critical residues or ranges."""
    auto_acmg_data.prot_pos = 300  # Not within any critical residues or ranges
    auto_acmg_data.hgnc_id = "HGNC:10294"  # RPE65 gene
    result = leber_congenital_amaurosis_predictor.predict_pm1(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for non-critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant does not meet the PM1 criteria" in result.summary
    ), "The summary should indicate the lack of criteria met."


@patch("src.vcep.leber_congenital_amaurosis.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, leber_congenital_amaurosis_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method if logic changes."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the specific logic
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = leber_congenital_amaurosis_predictor.predict_pm1(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."
