from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import GlaucomaPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def glaucoma_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return GlaucomaPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_not_applicable(glaucoma_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for MYOC in Glaucoma."""
    auto_acmg_data.hgnc_id = "HGNC:7610"  # MYOC gene
    result = glaucoma_predictor.predict_pm1(glaucoma_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for MYOC."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting for MYOC."
    assert (
        "PM1 is not applicable for MYOC" in result.summary
    ), "The summary should indicate that PM1 is not applicable for MYOC."


def test_predict_pm1_strength_level(glaucoma_predictor, auto_acmg_data):
    """Test the strength level for PM1 when not applicable."""
    auto_acmg_data.hgnc_id = "HGNC:7610"  # MYOC gene
    result = glaucoma_predictor.predict_pm1(glaucoma_predictor.seqvar, auto_acmg_data)

    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting when PM1 is not applicable."


def test_predict_pm1_name(glaucoma_predictor, auto_acmg_data):
    """Test the name of the criteria returned by the Glaucoma predictor."""
    auto_acmg_data.hgnc_id = "HGNC:7610"  # MYOC gene
    result = glaucoma_predictor.predict_pm1(glaucoma_predictor.seqvar, auto_acmg_data)

    assert result.name == "PM1", "The name of the criteria should be 'PM1'."
