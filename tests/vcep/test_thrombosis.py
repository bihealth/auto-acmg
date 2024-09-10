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
from src.vcep import ThrombosisPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def thrombosis_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return ThrombosisPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


# ---------------- PM1 ----------------


def test_predict_pm1_in_cysteine_residue(thrombosis_predictor, auto_acmg_data):
    """Test when the variant affects a critical cysteine residue in SERPINC1."""
    auto_acmg_data.hgnc_id = "HGNC:775"  # SERPINC1 gene
    auto_acmg_data.prot_pos = 40  # A cysteine residue involved in disulfide bridges
    result = thrombosis_predictor.predict_pm1(thrombosis_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
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
        result.prediction == AutoACMGPrediction.Applicable
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
        result.prediction == AutoACMGPrediction.Applicable
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
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met for non-critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not meet the PM1 criteria for SERPINC1" in result.summary
    ), "The summary should indicate the PM1 criteria."


@patch("src.vcep.thrombosis.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, thrombosis_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the Thrombosis VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default fallback for PM1.",
    )

    result = thrombosis_predictor.predict_pm1(thrombosis_predictor.seqvar, auto_acmg_data)
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default fallback for PM1." in result.summary
    ), "The summary should indicate the default fallback."


# ---------------- PM2, BA1, BS1, BS2 ----------------


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
def test_predict_pm2ba1bs1bs2(mock_super_method, thrombosis_predictor, auto_acmg_data, seqvar):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = thrombosis_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.00002
    assert auto_acmg_data.thresholds.ba1_benign == 0.002
    assert auto_acmg_data.thresholds.bs1_benign == 0.0002

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


# ---------------- PM4 & BP3 ----------------


def test_bp3_not_applicable(thrombosis_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = thrombosis_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


# ---------------- PP2 & BP1 ----------------


def test_predict_pp2bp1(thrombosis_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Thrombosis."""

    # Call the method under test
    pp2_result, bp1_result = thrombosis_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


# ---------------- PP3 & BP4 ----------------


def test_verify_pp3bp4_thresholds(thrombosis_predictor, auto_acmg_data):
    """Test that the thresholds for PP3/BP4 prediction are correctly set."""
    thrombosis_predictor.verify_pp3bp4(thrombosis_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.revel_pathogenic == 0.6
    assert auto_acmg_data.thresholds.revel_benign == 0.3


@patch.object(ThrombosisPredictor, "_is_pathogenic_score")
@patch.object(ThrombosisPredictor, "_is_benign_score")
@patch.object(ThrombosisPredictor, "_affect_spliceAI")
def test_verify_pp3bp4_prediction_logic(
    mock_affect_spliceAI,
    mock_is_benign_score,
    mock_is_pathogenic_score,
    thrombosis_predictor,
    auto_acmg_data,
):
    """Test the prediction logic for PP3 and BP4."""
    mock_is_pathogenic_score.return_value = True
    mock_is_benign_score.return_value = False
    mock_affect_spliceAI.side_effect = [
        True,
        False,
    ]  # First call True, second call False

    prediction, comment = thrombosis_predictor.verify_pp3bp4(
        thrombosis_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 is True
    assert prediction.BP4 is False


@pytest.mark.parametrize(
    "revel_score, spliceAI_scores, expected_pp3, expected_bp4",
    [
        (0.7, [0.6, 0.6, 0.6, 0.6], True, False),  # High REVEL score, high SpliceAI
        (0.2, [0.1, 0.1, 0.1, 0.1], False, True),  # Low REVEL score, low SpliceAI
        # (0.5, [0.3, 0.3, 0.3, 0.3], False, False),  # Intermediate scores
        (0.7, [0.1, 0.1, 0.1, 0.1], True, False),  # High REVEL score, low SpliceAI
        (0.2, [0.6, 0.6, 0.6, 0.6], True, False),  # Low REVEL score, high SpliceAI
    ],
)
def test_verify_pp3bp4_various_scenarios(
    thrombosis_predictor,
    auto_acmg_data,
    revel_score,
    spliceAI_scores,
    expected_pp3,
    expected_bp4,
):
    """Test different scenarios for PP3 and BP4 prediction."""
    auto_acmg_data.scores.dbnsfp.revel = revel_score
    auto_acmg_data.scores.cadd.spliceAI_acceptor_gain = spliceAI_scores[0]
    auto_acmg_data.scores.cadd.spliceAI_acceptor_loss = spliceAI_scores[1]
    auto_acmg_data.scores.cadd.spliceAI_donor_gain = spliceAI_scores[2]
    auto_acmg_data.scores.cadd.spliceAI_donor_loss = spliceAI_scores[3]

    prediction, _ = thrombosis_predictor.verify_pp3bp4(thrombosis_predictor.seqvar, auto_acmg_data)

    assert prediction.PP3 == expected_pp3
    assert prediction.BP4 == expected_bp4


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_missing_scores(thrombosis_predictor, auto_acmg_data):
    """Test behavior when scores are missing."""
    auto_acmg_data.scores.dbnsfp.revel = None

    prediction, comment = thrombosis_predictor.verify_pp3bp4(
        thrombosis_predictor.seqvar, auto_acmg_data
    )

    assert prediction is None
    assert "An error occurred during prediction" in comment


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_error_handling(thrombosis_predictor, auto_acmg_data):
    """Test error handling in verify_pp3bp4 method."""
    with patch.object(
        ThrombosisPredictor,
        "_is_pathogenic_score",
        side_effect=Exception("Test error"),
    ):
        prediction, comment = thrombosis_predictor.verify_pp3bp4(
            thrombosis_predictor.seqvar, auto_acmg_data
        )

        assert prediction is None
        assert "An error occurred during prediction" in comment
        assert "Test error" in comment


def test_verify_pp3bp4_spliceai_thresholds(thrombosis_predictor, auto_acmg_data):
    """Test that SpliceAI thresholds are correctly used during PP3/BP4 prediction."""
    with (
        patch.object(ThrombosisPredictor, "_is_pathogenic_score", return_value=False),
        patch.object(ThrombosisPredictor, "_is_benign_score", return_value=False),
        patch.object(ThrombosisPredictor, "_affect_spliceAI", return_value=False),
    ):
        thrombosis_predictor.verify_pp3bp4(thrombosis_predictor.seqvar, auto_acmg_data)

        # Check that default SpliceAI thresholds are used
        assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.1
        assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.1
        assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.1
        assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.1
