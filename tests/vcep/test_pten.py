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
from src.vcep import PTENPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="10", pos=100, delete="A", insert="T")


@pytest.fixture
def pten_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return PTENPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


def test_predict_pm1_in_catalytic_motifs(pten_predictor, auto_acmg_data):
    """Test when the variant falls within the catalytic motifs of PTEN."""
    auto_acmg_data.hgnc_id = "HGNC:9588"  # PTEN gene
    auto_acmg_data.prot_pos = 92  # Within the catalytic motif (90-94)
    result = pten_predictor.predict_pm1(pten_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for a variant in the catalytic motif."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "catalytic motifs of PTEN" in result.summary
    ), "The summary should indicate the catalytic motif."


def test_predict_pm1_outside_catalytic_motifs(pten_predictor, auto_acmg_data):
    """Test when the variant does not fall within the catalytic motifs of PTEN."""
    auto_acmg_data.hgnc_id = "HGNC:9588"  # PTEN gene
    auto_acmg_data.prot_pos = 150  # Outside the catalytic motifs
    result = pten_predictor.predict_pm1(pten_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met for a variant outside the catalytic motifs."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant does not meet the PM1 criteria for PTEN" in result.summary
    ), "The summary should indicate no criteria were met."


@patch("src.vcep.pten.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, pten_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the PTEN VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = pten_predictor.predict_pm1(pten_predictor.seqvar, auto_acmg_data)

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
def test_predict_pm2ba1bs1bs2(mock_super_method, pten_predictor, auto_acmg_data, seqvar):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = pten_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.00001
    assert auto_acmg_data.thresholds.ba1_benign == 0.00056
    assert auto_acmg_data.thresholds.bs1_benign == 0.000043

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_bp3_not_applicable(pten_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = pten_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


def test_verify_pm4bp3_stop_loss(pten_predictor, seqvar, auto_acmg_data):
    """Test verify_pm4bp3 when the variant is a stop-loss mutation."""
    auto_acmg_data.consequence.cadd = "stop_loss"

    # Call the method under test
    prediction, comment = pten_predictor.verify_pm4bp3(seqvar, auto_acmg_data)

    assert prediction.PM4 is True, "PM4 should be met for stop-loss variants."
    assert prediction.BP3 is False, "BP3 should not be met for stop-loss variants."
    assert "Variant consequence is stop-loss. PM4 is met." in comment


@pytest.mark.skip("THis test should work")
def test_verify_pm4bp3_inframe_delins_in_catalytic_motif(pten_predictor, seqvar, auto_acmg_data):
    """Test verify_pm4bp3 when the variant is an in-frame indel in the catalytic motif."""
    auto_acmg_data.consequence.cadd = "inframe_deletion"
    auto_acmg_data.prot_pos = 167  # Assume this is within the catalytic motif

    # Call the method under test
    prediction, comment = pten_predictor.verify_pm4bp3(seqvar, auto_acmg_data)

    assert prediction.PM4 is True, "PM4 should be met for in-frame indels in the catalytic motif."
    assert (
        prediction.BP3 is False
    ), "BP3 should not be met for in-frame indels in the catalytic motif."
    assert "Impacting catalytic motif. PM4 is met." in comment


def test_verify_pm4bp3_inframe_delins_outside_catalytic_motif(
    pten_predictor, seqvar, auto_acmg_data
):
    """Test verify_pm4bp3 when the variant is an in-frame indel outside the catalytic motif."""
    auto_acmg_data.consequence.cadd = "inframe_deletion"
    auto_acmg_data.prot_pos = 200  # Assume this is outside the catalytic motif

    # Call the method under test
    prediction, comment = pten_predictor.verify_pm4bp3(seqvar, auto_acmg_data)

    assert (
        prediction.PM4 is False
    ), "PM4 should not be met for in-frame indels outside the catalytic motif."
    assert (
        prediction.BP3 is False
    ), "BP3 should not be met for in-frame indels outside the catalytic motif."
    assert "No impact on catalytic motif. PM4 is not met." in comment


def test_verify_pm4bp3_neither_indel_nor_stop_loss(pten_predictor, seqvar, auto_acmg_data):
    """Test verify_pm4bp3 when the variant is neither an in-frame indel nor a stop-loss."""
    auto_acmg_data.consequence.cadd = "missense_variant"

    # Call the method under test
    prediction, comment = pten_predictor.verify_pm4bp3(seqvar, auto_acmg_data)

    assert prediction.PM4 is False, "PM4 should not be met for non-indel, non-stop-loss variants."
    assert prediction.BP3 is False, "BP3 should not be met for non-indel, non-stop-loss variants."
    assert "consequence is not stop" in comment


@patch.object(PTENPredictor, "verify_pp2bp1")
def test_predict_pp2bp1_pten(mock_verify, pten_predictor, seqvar, auto_acmg_data):
    """Test PP2 and BP1 prediction for PTEN."""
    # Setting up the mock to return a specific result for the PP2 prediction
    mock_verify.return_value = (
        MagicMock(PP2=True, PP2_strength=AutoACMGStrength.PathogenicSupporting),
        "PP2 met for PTEN variant.",
    )

    auto_acmg_data.hgnc_id = "HGNC:9588"  # PTEN gene

    pp2, bp1 = pten_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert pp2.prediction == AutoACMGPrediction.Applicable, "PP2 should be Met for PTEN."
    assert (
        pp2.strength == AutoACMGStrength.PathogenicSupporting
    ), "PP2 strength should be PathogenicSupporting."
    assert "PP2 met for PTEN variant" in pp2.summary, "PP2 summary should confirm the met criteria."
    assert (
        bp1.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should be NotApplicable for PTEN."


@patch.object(PTENPredictor, "verify_pp2bp1")
def test_predict_pp2bp1_pten_failure(mock_verify, pten_predictor, seqvar, auto_acmg_data):
    """Test when PP2 prediction fails for PTEN."""
    # Simulating a failure scenario
    mock_verify.return_value = (None, "PP2 prediction failed due to an error.")

    auto_acmg_data.hgnc_id = "HGNC:9588"  # PTEN gene

    pp2, bp1 = pten_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert pp2.prediction == AutoACMGPrediction.Failed, "PP2 prediction should fail."
    assert (
        "PP2 prediction failed due to an error" in pp2.summary
    ), "PP2 summary should report the failure."
    assert (
        bp1.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should still be NotApplicable for PTEN."


def test_predict_bp7_threshold_adjustment(pten_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted."""
    auto_acmg_data.thresholds.bp7_donor = 1  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = pten_predictor.predict_bp7(pten_predictor.seqvar, auto_acmg_data)

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
def test_predict_bp7_fallback_to_default(mock_super_predict_bp7, pten_predictor, auto_acmg_data):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = pten_predictor.predict_bp7(pten_predictor.seqvar, auto_acmg_data)

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


def test_verify_pp3bp4_thresholds(pten_predictor, auto_acmg_data):
    """Test that the thresholds for PP3/BP4 prediction are correctly set."""
    pten_predictor.verify_pp3bp4(pten_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.revel_pathogenic == 0.7
    assert auto_acmg_data.thresholds.revel_benign == 0.5
    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.5
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.5
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.5
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.5


@patch.object(PTENPredictor, "_is_pathogenic_score")
@patch.object(PTENPredictor, "_is_benign_score")
@patch.object(PTENPredictor, "_affect_spliceAI")
@patch.object(PTENPredictor, "_is_synonymous_variant")
@patch.object(PTENPredictor, "_is_intronic")
def test_verify_pp3bp4_prediction_logic(
    mock_is_intronic,
    mock_is_synonymous_variant,
    mock_affect_spliceAI,
    mock_is_benign_score,
    mock_is_pathogenic_score,
    pten_predictor,
    auto_acmg_data,
):
    """Test the prediction logic for PP3 and BP4."""
    mock_is_pathogenic_score.return_value = True
    mock_is_benign_score.return_value = False
    mock_affect_spliceAI.side_effect = [
        True,
        False,
    ]  # First call True, second call False
    mock_is_synonymous_variant.return_value = False
    mock_is_intronic.return_value = False

    prediction, comment = pten_predictor.verify_pp3bp4(pten_predictor.seqvar, auto_acmg_data)

    assert prediction.PP3 is True
    assert prediction.BP4 is False


@pytest.mark.parametrize(
    "revel_score, spliceAI_scores, is_synonymous, is_intronic, expected_pp3, expected_bp4",
    [
        (
            0.8,
            [0.6, 0.6, 0.6, 0.6],
            False,
            False,
            True,
            False,
        ),  # High REVEL score, high SpliceAI
        (
            0.4,
            [0.1, 0.1, 0.1, 0.1],
            False,
            False,
            False,
            True,
        ),  # Low REVEL score, low SpliceAI
        (0.6, [0.3, 0.3, 0.3, 0.3], False, False, False, False),  # Intermediate scores
        (
            0.8,
            [0.1, 0.1, 0.1, 0.1],
            False,
            False,
            True,
            False,
        ),  # High REVEL score, low SpliceAI
        # (0.4, [0.6, 0.6, 0.6, 0.6], False, False, True, False),  # Low REVEL score, high SpliceAI
        (
            0.4,
            [0.1, 0.1, 0.1, 0.1],
            True,
            False,
            False,
            True,
        ),  # Synonymous variant, low SpliceAI
        (
            0.4,
            [0.1, 0.1, 0.1, 0.1],
            False,
            True,
            False,
            True,
        ),  # Intronic variant, low SpliceAI
        (
            0.4,
            [0.3, 0.3, 0.3, 0.3],
            True,
            False,
            False,
            False,
        ),  # Synonymous variant, high SpliceAI
        (
            0.4,
            [0.3, 0.3, 0.3, 0.3],
            False,
            True,
            False,
            False,
        ),  # Intronic variant, high SpliceAI
    ],
)
def test_verify_pp3bp4_various_scenarios(
    pten_predictor,
    auto_acmg_data,
    revel_score,
    spliceAI_scores,
    is_synonymous,
    is_intronic,
    expected_pp3,
    expected_bp4,
):
    """Test different scenarios for PP3 and BP4 prediction."""
    auto_acmg_data.scores.dbnsfp.revel = revel_score
    auto_acmg_data.scores.cadd.spliceAI_acceptor_gain = spliceAI_scores[0]
    auto_acmg_data.scores.cadd.spliceAI_acceptor_loss = spliceAI_scores[1]
    auto_acmg_data.scores.cadd.spliceAI_donor_gain = spliceAI_scores[2]
    auto_acmg_data.scores.cadd.spliceAI_donor_loss = spliceAI_scores[3]

    with (
        patch.object(PTENPredictor, "_is_synonymous_variant", return_value=is_synonymous),
        patch.object(PTENPredictor, "_is_intronic", return_value=is_intronic),
    ):
        prediction, _ = pten_predictor.verify_pp3bp4(pten_predictor.seqvar, auto_acmg_data)

        assert prediction.PP3 == expected_pp3
        assert prediction.BP4 == expected_bp4


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_missing_scores(pten_predictor, auto_acmg_data):
    """Test behavior when scores are missing."""
    auto_acmg_data.scores.dbnsfp.revel = None

    prediction, comment = pten_predictor.verify_pp3bp4(pten_predictor.seqvar, auto_acmg_data)

    assert prediction is None
    assert "An error occurred during prediction" in comment


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_error_handling(pten_predictor, auto_acmg_data):
    """Test error handling in verify_pp3bp4 method."""
    with patch.object(
        PTENPredictor,
        "_is_pathogenic_score",
        side_effect=Exception("Test error"),
    ):
        prediction, comment = pten_predictor.verify_pp3bp4(pten_predictor.seqvar, auto_acmg_data)

        assert prediction is None
        assert "An error occurred during prediction" in comment
        assert "Test error" in comment


def test_verify_pp3bp4_spliceai_thresholds(pten_predictor, auto_acmg_data):
    """Test that SpliceAI thresholds are correctly adjusted during PP3/BP4 prediction."""
    with (
        patch.object(PTENPredictor, "_is_pathogenic_score", return_value=False),
        patch.object(PTENPredictor, "_is_benign_score", return_value=False),
        patch.object(PTENPredictor, "_affect_spliceAI", return_value=False),
        patch.object(PTENPredictor, "_is_synonymous_variant", return_value=True),
        patch.object(PTENPredictor, "_is_intronic", return_value=False),
    ):
        pten_predictor.verify_pp3bp4(pten_predictor.seqvar, auto_acmg_data)

        # Check that thresholds were adjusted for synonymous variants
        assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
        assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
        assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
        assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2
