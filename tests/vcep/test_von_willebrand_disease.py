from unittest.mock import MagicMock

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.von_willebrand_disease import VonWillebrandDiseasePredictor


@pytest.fixture
def seqvar():
    return SeqVar(
        genome_release=GenomeRelease.GRCh37, chrom="12", pos=6135437, delete="A", insert="T"
    )


@pytest.fixture
def vwf_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return VonWillebrandDiseasePredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_not_applicable(vwf_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for von Willebrand Disease."""
    auto_acmg_data.hgnc_id = "HGNC:12726"  # VWF gene
    result = vwf_predictor.predict_pm1(vwf_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be NotApplicable for VWF."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not applicable for" in result.summary
    ), "The summary should indicate PM1 is not applicable."


def test_bp3_not_applicable(vwf_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = vwf_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"
