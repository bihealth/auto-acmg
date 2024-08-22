from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import CerebralCreatineDeficiencySyndromesPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def cerebral_creatine_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return CerebralCreatineDeficiencySyndromesPredictor(
        seqvar=seqvar, result=result, config=MagicMock()
    )


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_not_applicable(cerebral_creatine_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for Cerebral Creatine Deficiency Syndromes."""
    result = cerebral_creatine_predictor.predict_pm1(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for Cerebral Creatine Deficiency Syndromes."
    assert (
        result.summary == "PM1 is not applicable for Cerebral Creatine Deficiency Syndromes."
    ), "The summary should indicate that PM1 is not applicable for Cerebral Creatine Deficiency Syndromes."


def test_predict_pm1_strength(cerebral_creatine_predictor, auto_acmg_data):
    """Test the strength level returned by the Cerebral Creatine Deficiency Syndromes predictor."""
    result = cerebral_creatine_predictor.predict_pm1(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate for Cerebral Creatine Deficiency Syndromes."


def test_predict_pm1_name(cerebral_creatine_predictor, auto_acmg_data):
    """Test the name of the criteria returned by the Cerebral Creatine Deficiency Syndromes predictor."""
    result = cerebral_creatine_predictor.predict_pm1(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    assert result.name == "PM1", "The name of the criteria should be 'PM1'."


@patch("src.vcep.cerebral_creatine_deficiency_syndromes.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, cerebral_creatine_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method (if implemented)."""
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = cerebral_creatine_predictor.predict_pm1(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    # In this specific case, the fallback should never happen since PM1 is always not applicable,
    # but this test ensures that if something changes, the fallback works correctly.
    assert result.prediction == AutoACMGPrediction.NotApplicable, "PM1 should remain NotApplicable."
