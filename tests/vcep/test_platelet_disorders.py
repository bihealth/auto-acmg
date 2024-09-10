from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
)
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor
from src.vcep import PlateletDisordersPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="17", pos=100, delete="A", insert="T")


@pytest.fixture
def platelet_disorders_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return PlateletDisordersPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


def test_predict_pm1_not_met_for_itga2b(platelet_disorders_predictor, auto_acmg_data):
    """Test when the variant is in the ITGA2B gene and PM1 is not met."""
    auto_acmg_data.hgnc_id = "HGNC:6138"  # ITGA2B gene
    result = platelet_disorders_predictor.predict_pm1(
        platelet_disorders_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met for ITGA2B."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not met for ITGA2B and ITGB3" in result.summary
    ), "The summary should explain the reason."


def test_predict_pm1_not_met_for_itgb3(platelet_disorders_predictor, auto_acmg_data):
    """Test when the variant is in the ITGB3 gene and PM1 is not met."""
    auto_acmg_data.hgnc_id = "HGNC:6156"  # ITGB3 gene
    result = platelet_disorders_predictor.predict_pm1(
        platelet_disorders_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.NotApplicable, "PM1 should not be met for ITGB3."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not met for ITGA2B and ITGB3" in result.summary
    ), "The summary should explain the reason."


@patch("src.vcep.platelet_disorders.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, platelet_disorders_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the Platelet Disorders VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = platelet_disorders_predictor.predict_pm1(
        platelet_disorders_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


@patch.object(
    DefaultSeqVarPredictor,
    "predict_pm2ba1bs1bs2",
    return_value=(
        AutoACMGCriteria(name="PM2"),
        AutoACMGCriteria(name="BA1"),
        AutoACMGCriteria(name="BS1"),
        AutoACMGCriteria(name="BS2"),
    ),
)
def test_predict_pm2ba1bs1bs2(
    mock_super_method, platelet_disorders_predictor, auto_acmg_data, seqvar
):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = platelet_disorders_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.0001
    assert auto_acmg_data.thresholds.ba1_benign == 0.0024
    assert auto_acmg_data.thresholds.bs1_benign == 0.00158

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_predict_pp2bp1(platelet_disorders_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Platelet Disorders."""

    # Call the method under test
    pp2_result, bp1_result = platelet_disorders_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


def test_predict_pp3bp4_revel_strategy(platelet_disorders_predictor, auto_acmg_data):
    """Test that REVEL is set as the strategy for PP3/BP4 prediction."""
    pp3_result, bp4_result = platelet_disorders_predictor.predict_pp3bp4(
        platelet_disorders_predictor.seqvar, auto_acmg_data
    )

    assert auto_acmg_data.thresholds.pp3bp4_strategy == "revel"
    assert auto_acmg_data.thresholds.revel_pathogenic == 0.7
    assert auto_acmg_data.thresholds.revel_benign == 0.25


@patch("src.vcep.platelet_disorders.DefaultSeqVarPredictor.predict_pp3bp4")
def test_predict_pp3bp4_calls_superclass(
    mock_super_predict_pp3bp4, platelet_disorders_predictor, auto_acmg_data
):
    """Test that the superclass method is called with the correct parameters."""
    mock_super_predict_pp3bp4.return_value = (
        AutoACMGCriteria(
            name="PP3",
            prediction=AutoACMGPrediction.Applicable,
            strength=AutoACMGStrength.PathogenicSupporting,
        ),
        AutoACMGCriteria(
            name="BP4",
            prediction=AutoACMGPrediction.NotApplicable,
            strength=AutoACMGStrength.BenignSupporting,
        ),
    )

    pp3_result, bp4_result = platelet_disorders_predictor.predict_pp3bp4(
        platelet_disorders_predictor.seqvar, auto_acmg_data
    )

    mock_super_predict_pp3bp4.assert_called_once_with(
        platelet_disorders_predictor.seqvar, auto_acmg_data
    )
    assert pp3_result.prediction == AutoACMGPrediction.Applicable
    assert bp4_result.prediction == AutoACMGPrediction.NotApplicable


@pytest.mark.parametrize(
    "revel_score, expected_pp3, expected_bp4",
    [
        (
            0.8,
            AutoACMGPrediction.Applicable,
            AutoACMGPrediction.NotApplicable,
        ),  # High REVEL score
        (
            0.5,
            AutoACMGPrediction.NotApplicable,
            AutoACMGPrediction.NotApplicable,
        ),  # Intermediate REVEL score
        (
            0.2,
            AutoACMGPrediction.NotApplicable,
            AutoACMGPrediction.Applicable,
        ),  # Low REVEL score
    ],
)
def test_predict_pp3bp4_revel_scenarios(
    platelet_disorders_predictor,
    auto_acmg_data,
    revel_score,
    expected_pp3,
    expected_bp4,
):
    """Test different REVEL score scenarios."""
    auto_acmg_data.scores.dbnsfp.revel = revel_score

    with patch(
        "src.vcep.platelet_disorders.DefaultSeqVarPredictor.predict_pp3bp4"
    ) as mock_super_predict_pp3bp4:
        mock_super_predict_pp3bp4.return_value = (
            AutoACMGCriteria(
                name="PP3",
                prediction=expected_pp3,
                strength=AutoACMGStrength.PathogenicSupporting,
            ),
            AutoACMGCriteria(
                name="BP4",
                prediction=expected_bp4,
                strength=AutoACMGStrength.BenignSupporting,
            ),
        )

        pp3_result, bp4_result = platelet_disorders_predictor.predict_pp3bp4(
            platelet_disorders_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == expected_pp3
        assert bp4_result.prediction == expected_bp4


def test_predict_pp3bp4_strength(platelet_disorders_predictor, auto_acmg_data):
    """Test that the strength of PP3 and BP4 is correctly set."""
    with patch(
        "src.vcep.platelet_disorders.DefaultSeqVarPredictor.predict_pp3bp4"
    ) as mock_super_predict_pp3bp4:
        mock_super_predict_pp3bp4.return_value = (
            AutoACMGCriteria(
                name="PP3",
                prediction=AutoACMGPrediction.Applicable,
                strength=AutoACMGStrength.PathogenicSupporting,
            ),
            AutoACMGCriteria(
                name="BP4",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.BenignSupporting,
            ),
        )

        pp3_result, bp4_result = platelet_disorders_predictor.predict_pp3bp4(
            platelet_disorders_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.strength == AutoACMGStrength.PathogenicSupporting
        assert bp4_result.strength == AutoACMGStrength.BenignSupporting


def test_predict_pp3bp4_no_revel_score(platelet_disorders_predictor, auto_acmg_data):
    """Test behavior when no REVEL score is available."""
    auto_acmg_data.scores.dbnsfp.revel = None

    with patch(
        "src.vcep.platelet_disorders.DefaultSeqVarPredictor.predict_pp3bp4"
    ) as mock_super_predict_pp3bp4:
        mock_super_predict_pp3bp4.return_value = (
            AutoACMGCriteria(
                name="PP3",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
            ),
            AutoACMGCriteria(
                name="BP4",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.BenignSupporting,
            ),
        )

        pp3_result, bp4_result = platelet_disorders_predictor.predict_pp3bp4(
            platelet_disorders_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == AutoACMGPrediction.NotApplicable
        assert bp4_result.prediction == AutoACMGPrediction.NotApplicable


def test_predict_pp3bp4_error_handling(platelet_disorders_predictor, auto_acmg_data):
    """Test error handling in predict_pp3bp4 method."""
    with patch(
        "src.vcep.platelet_disorders.DefaultSeqVarPredictor.predict_pp3bp4",
        side_effect=Exception("Test error"),
    ):
        with pytest.raises(Exception) as exc_info:
            platelet_disorders_predictor.predict_pp3bp4(
                platelet_disorders_predictor.seqvar, auto_acmg_data
            )

        assert str(exc_info.value) == "Test error"
