from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import (
    PS1PM5,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
)
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor
from src.vcep import TP53Predictor


@pytest.fixture
def seqvar():
    return SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="17",
        pos=7578406,
        delete="C",
        insert="T",
    )


@pytest.fixture
def tp53_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return TP53Predictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_overrides(mock_super_verify, tp53_predictor, seqvar, auto_acmg_data):
    """Test that the overridden verify_ps1pm5 method in TP53Predictor works correctly."""
    # Set up the mock to return PS1 and PM5 as applicable initially
    mock_super_verify.return_value = (PS1PM5(PS1=True, PM5=True), "Initial evaluation")

    # Setup the data
    auto_acmg_data.consequence = MagicMock(mehari=["missense_variant"])
    auto_acmg_data.thresholds = MagicMock(
        spliceAI_acceptor_gain=0.5,
        spliceAI_acceptor_loss=0.5,
        spliceAI_donor_gain=0.5,
        spliceAI_donor_loss=0.5,
    )
    auto_acmg_data.scores = MagicMock(
        cadd=MagicMock(spliceAI_acceptor_gain=0.3, spliceAI_donor_gain=0.3)
    )

    # Run the method under test
    prediction, comment = tp53_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that the splicing effect leads to overriding PS1 and PM5 as not applicable
    assert not prediction.PS1, "PS1 should be marked as not applicable due to splicing effect."
    assert not prediction.PM5, "PM5 should be marked as not applicable due to splicing effect."
    assert "Variant affects splicing" in comment, "Comment should note the splicing effect."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
    ), "SpliceAI acceptor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
    ), "SpliceAI acceptor loss threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
    ), "SpliceAI donor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2
    ), "SpliceAI donor loss threshold should be adjusted to 0.2"


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_no_splicing_effect(
    mock_super_verify, tp53_predictor, seqvar, auto_acmg_data
):
    """Test that PS1 and PM5 remain applicable when there's no splicing effect."""
    # Set up the mock to return PS1 and PM5 as applicable initially
    mock_super_verify.return_value = (PS1PM5(PS1=True, PM5=True), "Initial evaluation")

    # Setup the data with no splicing effect
    auto_acmg_data.consequence = MagicMock(mehari=["missense_variant"])
    auto_acmg_data.thresholds = MagicMock(
        spliceAI_acceptor_gain=0.5,
        spliceAI_acceptor_loss=0.5,
        spliceAI_donor_gain=0.5,
        spliceAI_donor_loss=0.5,
    )
    auto_acmg_data.scores = MagicMock(
        cadd=MagicMock(
            spliceAI_acceptor_gain=0.1,
            spliceAI_acceptor_loss=0.1,
            spliceAI_donor_gain=0.1,
            spliceAI_donor_loss=0.1,
        )
    )

    # Run the method under test
    prediction, comment = tp53_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that PS1 and PM5 remain applicable
    assert prediction.PS1, "PS1 should remain applicable when there's no splicing effect."
    assert prediction.PM5, "PM5 should remain applicable when there's no splicing effect."
    assert "Initial evaluation" in comment, "Comment should reflect the initial evaluation."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
    ), "SpliceAI acceptor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
    ), "SpliceAI acceptor loss threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
    ), "SpliceAI donor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2
    ), "SpliceAI donor loss threshold should be adjusted to 0.2"


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_non_missense(mock_super_verify, tp53_predictor, seqvar, auto_acmg_data):
    """Test that PS1 and PM5 remain as per superclass for non-missense variants."""
    # Set up the mock to return PS1 and PM5 as not applicable initially
    mock_super_verify.return_value = (
        PS1PM5(PS1=False, PM5=False),
        "Not applicable for non-missense",
    )

    # Setup the data with a non-missense variant
    auto_acmg_data.consequence = MagicMock(mehari=["synonymous_variant"])

    # Run the method under test
    prediction, comment = tp53_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that PS1 and PM5 remain as per superclass evaluation
    assert not prediction.PS1, "PS1 should remain not applicable for non-missense variants."
    assert not prediction.PM5, "PM5 should remain not applicable for non-missense variants."
    assert (
        "Not applicable for non-missense" in comment
    ), "Comment should reflect superclass evaluation."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
    ), "SpliceAI acceptor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
    ), "SpliceAI acceptor loss threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
    ), "SpliceAI donor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2
    ), "SpliceAI donor loss threshold should be adjusted to 0.2"


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_exception_handling(
    mock_super_verify, tp53_predictor, seqvar, auto_acmg_data
):
    """Test verify_ps1pm5 method exception handling."""
    mock_super_verify.side_effect = Exception("Test exception")

    with pytest.raises(Exception) as exc_info:
        tp53_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    assert "Test exception" in str(exc_info.value), "Should raise the original exception"


def test_predict_pm1_in_moderate_residue(tp53_predictor, auto_acmg_data):
    """Test when the variant affects a critical residue in TP53."""
    auto_acmg_data.hgnc_id = "HGNC:11998"  # TP53 gene
    auto_acmg_data.prot_pos = 175  # A critical residue in TP53
    result = tp53_predictor.predict_pm1(tp53_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
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
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met for non-critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not meet the PM1 criteria for TP53" in result.summary
    ), "The summary should indicate the PM1 criteria."


@patch("src.vcep.tp53.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, tp53_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the TP53 VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default fallback for PM1.",
    )
    result = tp53_predictor.predict_pm1(tp53_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met in the default fallback."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Default fallback for PM1." in result.summary
    ), "The summary should indicate the fallback."


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
def test_predict_pm2ba1bs1bs2(mock_super_method, tp53_predictor, auto_acmg_data, seqvar):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = tp53_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.00003
    assert auto_acmg_data.thresholds.ba1_benign == 0.001
    assert auto_acmg_data.thresholds.bs1_benign == 0.0003

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_predict_pm4bp3_tp53(tp53_predictor, seqvar, auto_acmg_data):
    """Test the predict_pm4bp3 method for TP53 VCEP."""
    # Call the method under test
    pm4_result, bp3_result = tp53_predictor.predict_pm4bp3(seqvar, auto_acmg_data)

    # Check PM4 result
    assert isinstance(
        pm4_result, AutoACMGCriteria
    ), "The PM4 result should be of type AutoACMGCriteria."
    assert (
        pm4_result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM4 should be NotApplicable for TP53."
    assert (
        "PM4 is not applicable for TP53 VCEP." in pm4_result.summary
    ), "The summary should indicate PM4 is not applicable for TP53."

    # Check BP3 result
    assert isinstance(
        bp3_result, AutoACMGCriteria
    ), "The BP3 result should be of type AutoACMGCriteria."
    assert (
        bp3_result.prediction == AutoACMGPrediction.NotApplicable
    ), "BP3 should be NotApplicable for TP53."
    assert (
        "BP3 is not applicable for TP53 VCEP." in bp3_result.summary
    ), "The summary should indicate BP3 is not applicable for TP53."


def test_predict_pp2bp1(tp53_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for TP53."""
    # Call the method under test
    pp2_result, bp1_result = tp53_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


def test_predict_bp7_threshold_adjustment(tp53_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted."""
    auto_acmg_data.thresholds.bp7_donor = 1  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = tp53_predictor.predict_bp7(tp53_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.bp7_donor == 7
    ), "The BP7 donor threshold should be adjusted to 7."
    assert (
        auto_acmg_data.thresholds.bp7_acceptor == 21
    ), "The BP7 acceptor threshold should be adjusted to 4."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(DefaultSeqVarPredictor, "predict_bp7")
def test_predict_bp7_fallback_to_default(mock_super_predict_bp7, tp53_predictor, auto_acmg_data):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = tp53_predictor.predict_bp7(tp53_predictor.seqvar, auto_acmg_data)

    # Verify the result and ensure the superclass method was called
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "BP7 should return NotMet as mocked."
    assert (
        result.strength == AutoACMGStrength.BenignSupporting
    ), "The strength should be BenignSupporting."
    assert (
        "Default BP7 prediction fallback." in result.summary
    ), "The summary should indicate the fallback."
    assert mock_super_predict_bp7.called, "super().predict_bp7 should have been called."


def test_verify_pp3bp4_thresholds(tp53_predictor, auto_acmg_data):
    """Test that the thresholds for PP3/BP4 prediction are correctly set."""
    tp53_predictor.verify_pp3bp4(tp53_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.bayesDel_noAF_pathogenic == 0.16
    assert auto_acmg_data.thresholds.bayesDel_noAF_benign == 0.16
    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.1
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.1
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.1
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.1


@patch.object(TP53Predictor, "_is_pathogenic_score")
@patch.object(TP53Predictor, "_is_benign_score")
@patch.object(TP53Predictor, "_affect_spliceAI")
def test_verify_pp3bp4_prediction_logic(
    mock_affect_spliceAI,
    mock_is_benign_score,
    mock_is_pathogenic_score,
    tp53_predictor,
    auto_acmg_data,
):
    """Test the prediction logic for PP3 and BP4."""
    mock_is_pathogenic_score.return_value = True
    mock_is_benign_score.return_value = False
    mock_affect_spliceAI.side_effect = [
        True,
        False,
    ]  # First call True, second call False

    prediction, comment = tp53_predictor.verify_pp3bp4(tp53_predictor.seqvar, auto_acmg_data)

    assert prediction.PP3 is True
    assert prediction.BP4 is False


@pytest.mark.parametrize(
    "bayesDel_score, spliceAI_scores, expected_pp3, expected_bp4",
    [
        (0.2, [0.3, 0.3, 0.3, 0.3], True, False),  # High BayesDel score, high SpliceAI
        (
            0.1,
            [0.05, 0.05, 0.05, 0.05],
            False,
            True,
        ),  # Low BayesDel score, low SpliceAI
        # (0.16, [0.15, 0.15, 0.15, 0.15], False, False),  # Intermediate scores
        (
            0.2,
            [0.05, 0.05, 0.05, 0.05],
            True,
            False,
        ),  # High BayesDel score, low SpliceAI
        (0.1, [0.3, 0.3, 0.3, 0.3], True, False),  # Low BayesDel score, high SpliceAI
    ],
)
def test_verify_pp3bp4_various_scenarios(
    tp53_predictor,
    auto_acmg_data,
    bayesDel_score,
    spliceAI_scores,
    expected_pp3,
    expected_bp4,
):
    """Test different scenarios for PP3 and BP4 prediction."""
    auto_acmg_data.scores.dbnsfp.bayesDel_noAF = bayesDel_score
    auto_acmg_data.scores.cadd.spliceAI_acceptor_gain = spliceAI_scores[0]
    auto_acmg_data.scores.cadd.spliceAI_acceptor_loss = spliceAI_scores[1]
    auto_acmg_data.scores.cadd.spliceAI_donor_gain = spliceAI_scores[2]
    auto_acmg_data.scores.cadd.spliceAI_donor_loss = spliceAI_scores[3]

    prediction, _ = tp53_predictor.verify_pp3bp4(tp53_predictor.seqvar, auto_acmg_data)

    assert prediction.PP3 == expected_pp3
    assert prediction.BP4 == expected_bp4


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_missing_scores(tp53_predictor, auto_acmg_data):
    """Test behavior when scores are missing."""
    auto_acmg_data.scores.dbnsfp.bayesDel_noAF = None

    prediction, comment = tp53_predictor.verify_pp3bp4(tp53_predictor.seqvar, auto_acmg_data)

    assert prediction is None
    assert "An error occurred during prediction" in comment


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_error_handling(tp53_predictor, auto_acmg_data):
    """Test error handling in verify_pp3bp4 method."""
    with patch.object(
        TP53Predictor,
        "_is_pathogenic_score",
        side_effect=Exception("Test error"),
    ):
        prediction, comment = tp53_predictor.verify_pp3bp4(tp53_predictor.seqvar, auto_acmg_data)

        assert prediction is None
        assert "An error occurred during prediction" in comment
        assert "Test error" in comment


def test_verify_pp3bp4_spliceai_thresholds(tp53_predictor, auto_acmg_data):
    """Test that SpliceAI thresholds are correctly adjusted during PP3/BP4 prediction."""
    with (
        patch.object(TP53Predictor, "_is_pathogenic_score", return_value=False),
        patch.object(TP53Predictor, "_is_benign_score", return_value=False),
        patch.object(TP53Predictor, "_affect_spliceAI", return_value=False),
    ):
        tp53_predictor.verify_pp3bp4(tp53_predictor.seqvar, auto_acmg_data)

        # Check that thresholds were adjusted for BP4
        assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.1
        assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.1
        assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.1
        assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.1
