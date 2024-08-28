from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.lysosomal_diseases import LysosomalDiseasesPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="17", pos=100, delete="A", insert="T")


@pytest.fixture
def lysosomal_diseases_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return LysosomalDiseasesPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_met_for_critical_residue(lysosomal_diseases_predictor, auto_acmg_data):
    """Test when the variant falls within the critical residues for GAA."""
    auto_acmg_data.prot_pos = 282  # Critical residue for GAA
    auto_acmg_data.hgnc_id = "HGNC:4065"  # GAA gene
    result = lysosomal_diseases_predictor.predict_pm1(
        lysosomal_diseases_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant affects a residue in GAA at position 282" in result.summary
    ), "The summary should indicate the affected residue."


def test_predict_pm1_not_met(lysosomal_diseases_predictor, auto_acmg_data):
    """Test when the variant does not fall within any critical residues for GAA."""
    auto_acmg_data.prot_pos = 700  # Not within any critical residues
    auto_acmg_data.hgnc_id = "HGNC:4065"  # GAA gene
    result = lysosomal_diseases_predictor.predict_pm1(
        lysosomal_diseases_predictor.seqvar, auto_acmg_data
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


@patch("src.vcep.lysosomal_diseases.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, lysosomal_diseases_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method if logic changes."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the specific logic
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = lysosomal_diseases_predictor.predict_pm1(
        lysosomal_diseases_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


def test_bp3_not_applicable(lysosomal_diseases_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = lysosomal_diseases_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"
