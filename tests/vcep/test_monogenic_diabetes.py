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
from src.vcep import MonogenicDiabetesPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="12", pos=100, delete="A", insert="T")


@pytest.fixture
def monogenic_diabetes_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return MonogenicDiabetesPredictor(seqvar=seqvar, result=result)


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


@pytest.fixture
def auto_acmg_data_gck():
    data = AutoACMGSeqVarData(hgnc_id="HGNC:4195")
    data.consequence = MagicMock(mehari=[], cadd=None)
    return data


# -------------- PM1 --------------


def test_predict_pm1_moderate_criteria_hnf1a(monogenic_diabetes_predictor, auto_acmg_data):
    """Test when variant falls within a moderate level residue for HNF1A."""
    auto_acmg_data.hgnc_id = "HGNC:11621"  # HNF1A gene
    auto_acmg_data.prot_pos = 130  # Within the critical residues for HNF1A
    result = monogenic_diabetes_predictor.predict_pm1(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for critical residues."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "critical for Monogenic Diabetes" in result.summary
    ), "The summary should indicate the critical residue."


def test_predict_pm1_supporting_criteria_hnf4a_domain(monogenic_diabetes_predictor, auto_acmg_data):
    """Test when variant falls within a supporting level domain for HNF4A."""
    auto_acmg_data.hgnc_id = "HGNC:5024"  # HNF4A gene
    auto_acmg_data.prot_pos = 100  # Within the critical domain for HNF4A
    result = monogenic_diabetes_predictor.predict_pm1(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for critical domains."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "falls within a critical region" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_supporting_criteria_hnf4a_promoter(
    monogenic_diabetes_predictor, auto_acmg_data
):
    """Test when variant falls within a supporting level promoter region for HNF4A."""
    auto_acmg_data.hgnc_id = "HGNC:5024"  # HNF4A gene
    auto_acmg_data.cds_pos = -140  # Within the critical promoter region for HNF4A
    result = monogenic_diabetes_predictor.predict_pm1(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for critical promoter regions."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "falls within a critical promoter region" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_moderate_criteria_gck(monogenic_diabetes_predictor, auto_acmg_data):
    """Test when variant falls within a moderate level residue for GCK."""
    auto_acmg_data.hgnc_id = "HGNC:4195"  # GCK gene
    auto_acmg_data.prot_pos = 151  # Within the critical residues for GCK
    result = monogenic_diabetes_predictor.predict_pm1(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for critical residues."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "critical for Monogenic Diabetes" in result.summary
    ), "The summary should indicate the critical residue."


def test_predict_pm1_not_met(monogenic_diabetes_predictor, auto_acmg_data):
    """Test when variant does not meet any PM1 criteria for Monogenic Diabetes genes."""
    auto_acmg_data.hgnc_id = "HGNC:4195"  # GCK gene
    auto_acmg_data.prot_pos = 500  # Outside all critical regions
    result = monogenic_diabetes_predictor.predict_pm1(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met for non-critical residues."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not meet the PM1 criteria" in result.summary
    ), "The summary should indicate that criteria were not met."


@patch("src.vcep.monogenic_diabetes.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, monogenic_diabetes_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the PM1_CLUSTER
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = monogenic_diabetes_predictor.predict_pm1(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


# -------------- PM2, BA1, BS1, BS2 --------------


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
@pytest.mark.parametrize(
    "hgnc_id,expected_bs1",
    [(None, 0.000033), ("HGNC:4195", 0.00004)],  # Default case  # GCK specific case
)
def test_predict_pm2ba1bs1bs2_adjustments(
    mock_super_method,
    monogenic_diabetes_predictor,
    auto_acmg_data,
    seqvar,
    hgnc_id,
    expected_bs1,
):
    # Setup
    auto_acmg_data.hgnc_id = hgnc_id

    # Call the method under test
    monogenic_diabetes_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds are set correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.000003
    assert auto_acmg_data.thresholds.ba1_benign == 0.0001
    assert auto_acmg_data.thresholds.bs1_benign == expected_bs1

    # Validate that the superclass method was called correctly with modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Reset mock for the next iteration
    mock_super_method.reset_mock()


# -------------- PM4 & BP3 --------------


def test_bp3_not_applicable(monogenic_diabetes_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = monogenic_diabetes_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


# -------------- PP2 & BP1 --------------


@patch.object(MonogenicDiabetesPredictor, "_is_missense")
def test_predict_pp2bp1_gck_missense(
    mock_is_missense, monogenic_diabetes_predictor, seqvar, auto_acmg_data_gck
):
    """Test predict_pp2bp1 for GCK where the variant is missense."""
    mock_is_missense.return_value = True
    pp2, bp1 = monogenic_diabetes_predictor.predict_pp2bp1(seqvar, auto_acmg_data_gck)

    assert (
        pp2.prediction == AutoACMGPrediction.Applicable
    ), "PP2 should be Met for a missense variant in GCK."
    assert (
        pp2.strength == AutoACMGStrength.PathogenicSupporting
    ), "PP2 strength should be PathogenicSupporting."
    assert "missense variant" in pp2.summary, "PP2 summary should confirm the missense nature."
    assert (
        bp1.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should be NotApplicable for GCK."


@patch.object(MonogenicDiabetesPredictor, "_is_missense")
def test_predict_pp2bp1_gck_non_missense(
    mock_is_missense, monogenic_diabetes_predictor, seqvar, auto_acmg_data_gck
):
    """Test predict_pp2bp1 for GCK where the variant is not missense."""
    mock_is_missense.return_value = False
    pp2, bp1 = monogenic_diabetes_predictor.predict_pp2bp1(seqvar, auto_acmg_data_gck)

    assert (
        pp2.prediction == AutoACMGPrediction.NotApplicable
    ), "PP2 should be NotMet for non-missense variants in GCK."
    assert (
        "not a missense variant" in pp2.summary
    ), "PP2 summary should confirm the non-missense nature."
    assert (
        bp1.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should still be NotApplicable for GCK."


# -------------- PP3 & BP4 --------------


@patch.object(MonogenicDiabetesPredictor, "_is_pathogenic_score")
@patch.object(MonogenicDiabetesPredictor, "_is_benign_score")
@patch.object(MonogenicDiabetesPredictor, "_affect_spliceAI")
@patch.object(MonogenicDiabetesPredictor, "_is_splice_variant")
@patch.object(MonogenicDiabetesPredictor, "_is_synonymous")
def test_verify_pp3bp4_prediction_logic(
    mock_is_synonymous,
    mock_is_splice_variant,
    mock_affect_spliceAI,
    mock_is_benign_score,
    mock_is_pathogenic_score,
    monogenic_diabetes_predictor,
    auto_acmg_data,
):
    """Test the prediction logic for PP3 and BP4."""
    mock_is_pathogenic_score.return_value = True
    mock_is_benign_score.return_value = False
    mock_affect_spliceAI.return_value = True
    mock_is_splice_variant.return_value = False
    mock_is_synonymous.return_value = False

    prediction, comment = monogenic_diabetes_predictor.verify_pp3bp4(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 is True
    assert prediction.BP4 is False


@pytest.mark.parametrize(
    "revel_score, spliceAI_scores, is_splice, is_synonymous, expected_pp3, expected_bp4",
    [
        (
            0.8,
            [0.3, 0.3, 0.3, 0.3],
            False,
            False,
            True,
            False,
        ),  # High REVEL score, high SpliceAI
        (
            0.1,
            [0.1, 0.1, 0.1, 0.1],
            False,
            False,
            False,
            True,
        ),  # Low REVEL score, low SpliceAI
        (
            0.5,
            [0.15, 0.15, 0.15, 0.15],
            False,
            False,
            False,
            False,
        ),  # Intermediate scores
        (
            0.8,
            [0.1, 0.1, 0.1, 0.1],
            False,
            False,
            True,
            False,
        ),  # High REVEL score, low SpliceAI
        # (0.1, [0.3, 0.3, 0.3, 0.3], False, False, True, False),  # Low REVEL score, high SpliceAI
        (
            0.1,
            [0.1, 0.1, 0.1, 0.1],
            True,
            False,
            False,
            True,
        ),  # Splice variant, low SpliceAI
        (
            0.1,
            [0.1, 0.1, 0.1, 0.1],
            False,
            True,
            False,
            True,
        ),  # Synonymous variant, low SpliceAI
    ],
)
def test_verify_pp3bp4_various_scenarios(
    monogenic_diabetes_predictor,
    auto_acmg_data,
    revel_score,
    spliceAI_scores,
    is_splice,
    is_synonymous,
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
        patch.object(MonogenicDiabetesPredictor, "_is_splice_variant", return_value=is_splice),
        patch.object(MonogenicDiabetesPredictor, "_is_synonymous", return_value=is_synonymous),
    ):
        prediction, _ = monogenic_diabetes_predictor.verify_pp3bp4(
            monogenic_diabetes_predictor.seqvar, auto_acmg_data
        )

        assert prediction.PP3 == expected_pp3
        assert prediction.BP4 == expected_bp4


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_missing_scores(monogenic_diabetes_predictor, auto_acmg_data):
    """Test behavior when scores are missing."""
    auto_acmg_data.scores.dbnsfp.revel = None

    prediction, comment = monogenic_diabetes_predictor.verify_pp3bp4(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
    )

    assert prediction is None
    assert "An error occurred during prediction" in comment


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_error_handling(monogenic_diabetes_predictor, auto_acmg_data):
    """Test error handling in verify_pp3bp4 method."""
    with patch.object(
        MonogenicDiabetesPredictor,
        "_is_pathogenic_score",
        side_effect=Exception("Test error"),
    ):
        prediction, comment = monogenic_diabetes_predictor.verify_pp3bp4(
            monogenic_diabetes_predictor.seqvar, auto_acmg_data
        )

        assert prediction is None
        assert "An error occurred during prediction" in comment
        assert "Test error" in comment


def test_verify_pp3bp4_spliceai_thresholds(monogenic_diabetes_predictor, auto_acmg_data):
    """Test that SpliceAI thresholds are correctly adjusted during PP3/BP4 prediction."""
    with (
        patch.object(MonogenicDiabetesPredictor, "_is_pathogenic_score", return_value=False),
        patch.object(MonogenicDiabetesPredictor, "_is_benign_score", return_value=False),
        patch.object(MonogenicDiabetesPredictor, "_affect_spliceAI", return_value=False),
        patch.object(MonogenicDiabetesPredictor, "_is_splice_variant", return_value=False),
        patch.object(MonogenicDiabetesPredictor, "_is_synonymous", return_value=False),
    ):
        monogenic_diabetes_predictor.verify_pp3bp4(
            monogenic_diabetes_predictor.seqvar, auto_acmg_data
        )

        # Check that thresholds were adjusted for PP3 and BP4
        assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
        assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
        assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
        assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2


# -------------- BP7 --------------


def test_predict_bp7_threshold_adjustment(monogenic_diabetes_predictor, auto_acmg_data):
    """Test that the BP7 thresholds are correctly adjusted for Monogenic Diabetes."""
    # Initial threshold values
    auto_acmg_data.thresholds.spliceAI_acceptor_gain = 0.1
    auto_acmg_data.thresholds.spliceAI_acceptor_loss = 0.1
    auto_acmg_data.thresholds.spliceAI_donor_gain = 0.1
    auto_acmg_data.thresholds.spliceAI_donor_loss = 0.1
    auto_acmg_data.thresholds.phyloP100 = 1.0

    # Call predict_bp7 method
    result = monogenic_diabetes_predictor.predict_bp7(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
    )

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
    ), "The spliceAI acceptor gain threshold should be adjusted to 0.2."
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
    ), "The spliceAI acceptor loss threshold should be adjusted to 0.2."
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
    ), "The spliceAI donor gain threshold should be adjusted to 0.2."
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2
    ), "The spliceAI donor loss threshold should be adjusted to 0.2."
    assert (
        auto_acmg_data.thresholds.phyloP100 == 2.0
    ), "The phyloP100 threshold should be adjusted to 2.0."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(MonogenicDiabetesPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default(
    mock_super_predict_bp7, monogenic_diabetes_predictor, auto_acmg_data
):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = monogenic_diabetes_predictor.predict_bp7(
        monogenic_diabetes_predictor.seqvar, auto_acmg_data
    )

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


def test_verify_pp3bp4_thresholds(monogenic_diabetes_predictor, auto_acmg_data):
    """Test that the thresholds for PP3/BP4 prediction are correctly set."""
    monogenic_diabetes_predictor.verify_pp3bp4(monogenic_diabetes_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.revel_pathogenic == 0.7
    assert auto_acmg_data.thresholds.revel_benign == 0.15
    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2
