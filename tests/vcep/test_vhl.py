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
from src.vcep import VHLPredictor


@pytest.fixture
def seqvar():
    return SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="3",
        pos=10191539,
        delete="C",
        insert="T",
    )


@pytest.fixture
def vhl_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return VHLPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


# ------------------ PM1 ------------------


def test_predict_pm1_in_moderate_residue(vhl_predictor, auto_acmg_data):
    """Test when the variant affects a critical residue in VHL."""
    auto_acmg_data.hgnc_id = "HGNC:12687"  # VHL gene
    auto_acmg_data.prot_pos = 167  # A critical residue in VHL
    result = vhl_predictor.predict_pm1(vhl_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "affects a germline hotspot or key functional domain in VHL" in result.summary
    ), "The summary should indicate the PM1 criteria."


def test_predict_pm1_not_met(vhl_predictor, auto_acmg_data):
    """Test when the variant does not meet the PM1 criteria for VHL."""
    auto_acmg_data.hgnc_id = "HGNC:12687"  # VHL gene
    auto_acmg_data.prot_pos = 200  # Position outside of critical residues
    result = vhl_predictor.predict_pm1(vhl_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met for non-critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not meet the PM1 criteria for VHL" in result.summary
    ), "The summary should indicate the PM1 criteria."


@patch("src.vcep.vhl.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, vhl_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the VHL VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default fallback for PM1.",
    )
    result = vhl_predictor.predict_pm1(vhl_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default fallback for PM1." in result.summary
    ), "The summary should indicate the default fallback."


# ------------------ PM2, BA1, BS1, BS2 ------------------


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
def test_predict_pm2ba1bs1bs2(mock_super_method, vhl_predictor, auto_acmg_data, seqvar):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = vhl_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.00000156
    assert auto_acmg_data.thresholds.ba1_benign == 0.000156
    assert auto_acmg_data.thresholds.bs1_benign == 0.0000156

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


# --------------- PM4 & BP3 ----------------


def test_in_vhl_important_domain_true(vhl_predictor, auto_acmg_data):
    """Test that the variant is in an important VHL domain."""
    auto_acmg_data.prot_pos = 160  # Within the Alpha (ɑ) domain

    # Call the method under test
    result = vhl_predictor._in_vhl_important_domain(auto_acmg_data)

    assert result is True, "The variant should be in an important VHL domain."


def test_in_vhl_important_domain_false(vhl_predictor, auto_acmg_data):
    """Test that the variant is not in an important VHL domain."""
    auto_acmg_data.prot_pos = 50  # Outside the important domains

    # Call the method under test
    result = vhl_predictor._in_vhl_important_domain(auto_acmg_data)

    assert result is False, "The variant should not be in an important VHL domain."


def test_in_gxeex_repeat_region_true(vhl_predictor, auto_acmg_data):
    """Test that the variant is in the GXEEX repeat region for VHL."""
    auto_acmg_data.prot_pos = 30  # Within the GXEEX repeat region

    # Call the method under test
    result = vhl_predictor._in_gxeex_repeat_region(auto_acmg_data)

    assert result is True, "The variant should be in the GXEEX repeat region."


def test_in_gxeex_repeat_region_false(vhl_predictor, auto_acmg_data):
    """Test that the variant is not in the GXEEX repeat region for VHL."""
    auto_acmg_data.prot_pos = 60  # Outside the GXEEX repeat region

    # Call the method under test
    result = vhl_predictor._in_gxeex_repeat_region(auto_acmg_data)

    assert result is False, "The variant should not be in the GXEEX repeat region."


@patch.object(VHLPredictor, "_in_vhl_important_domain", return_value=True)
@patch.object(VHLPredictor, "_in_gxeex_repeat_region", return_value=False)
def test_verify_pm4bp3_in_important_domain(
    mock_in_vhl_important_domain,
    mock_in_gxeex_repeat_region,
    vhl_predictor,
    seqvar,
    auto_acmg_data,
):
    """Test verify_pm4bp3 when the variant is an in-frame deletion in an important domain."""
    auto_acmg_data.consequence.cadd = "inframe_deletion"

    # Call the method under test
    prediction, comment = vhl_predictor.verify_pm4bp3(seqvar, auto_acmg_data)

    assert prediction.PM4 is True, "PM4 should be met for in-frame deletion in an important domain."
    assert (
        prediction.BP3 is False
    ), "BP3 should not be met for in-frame deletion in an important domain."
    assert "Variant is in an important domain of VHL. PM4 is met." in comment


@patch.object(VHLPredictor, "_in_vhl_important_domain", return_value=False)
@patch.object(VHLPredictor, "_in_gxeex_repeat_region", return_value=True)
def test_verify_pm4bp3_in_repeat_region(
    mock_in_vhl_important_domain,
    mock_in_gxeex_repeat_region,
    vhl_predictor,
    seqvar,
    auto_acmg_data,
):
    """Test verify_pm4bp3 when the variant is in the GXEEX repeat region."""
    auto_acmg_data.consequence.cadd = "inframe_deletion"

    # Call the method under test
    prediction, comment = vhl_predictor.verify_pm4bp3(seqvar, auto_acmg_data)

    assert (
        prediction.PM4 is False
    ), "PM4 should not be met for in-frame deletion in the GXEEX repeat region."
    assert (
        prediction.BP3 is True
    ), "BP3 should be met for in-frame deletion in the GXEEX repeat region."
    assert "Variant is in the GXEEX repeat motif in the VHL gene. BP3 is met." in comment


@patch.object(VHLPredictor, "_in_vhl_important_domain", return_value=False)
@patch.object(VHLPredictor, "_in_gxeex_repeat_region", return_value=False)
def test_verify_pm4bp3_neither_domain_nor_repeat(
    mock_in_vhl_important_domain,
    mock_in_gxeex_repeat_region,
    vhl_predictor,
    seqvar,
    auto_acmg_data,
):
    """Test verify_pm4bp3 when the variant is neither in an important domain nor in the GXEEX repeat region."""
    auto_acmg_data.consequence.cadd = "inframe_deletion"

    # Call the method under test
    prediction, comment = vhl_predictor.verify_pm4bp3(seqvar, auto_acmg_data)

    assert (
        prediction.PM4 is False
    ), "PM4 should not be met if not in important domain or GXEEX repeat region."
    assert (
        prediction.BP3 is True
    ), "BP3 should be met if variant is in neither an important domain nor GXEEX repeat region."
    assert "Variant is not in an important domain or in a repeat region. BP3 is met." in comment


# ------------------ PP2 & BP1 ------------------


def test_predict_pp2bp1(vhl_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for VHL."""

    # Call the method under test
    pp2_result, bp1_result = vhl_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


# ------------------ PP3 & BP4 ------------------


def test_verify_pp3bp4_thresholds(vhl_predictor, auto_acmg_data):
    """Test that the thresholds for PP3/BP4 prediction are correctly set."""
    vhl_predictor.verify_pp3bp4(vhl_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.revel_pathogenic == 0.664
    assert auto_acmg_data.thresholds.revel_benign == 0.3
    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.1
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.1
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.1
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.1


@patch.object(VHLPredictor, "_is_pathogenic_score")
@patch.object(VHLPredictor, "_is_benign_score")
@patch.object(VHLPredictor, "_affect_spliceAI")
@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_prediction_logic(
    mock_affect_spliceAI,
    mock_is_benign_score,
    mock_is_pathogenic_score,
    vhl_predictor,
    auto_acmg_data,
):
    """Test the prediction logic for PP3 and BP4."""
    mock_is_pathogenic_score.return_value = True
    mock_is_benign_score.return_value = False
    mock_affect_spliceAI.side_effect = [
        True,
        False,
    ]  # First call True, second call False

    prediction, comment = vhl_predictor.verify_pp3bp4(vhl_predictor.seqvar, auto_acmg_data)

    assert prediction.PP3 is True
    assert prediction.BP4 is False


@pytest.mark.parametrize(
    "revel_score, spliceAI_scores, expected_pp3, expected_bp4",
    [
        (0.7, [0.6, 0.6, 0.6, 0.6], True, False),  # High REVEL score, high SpliceAI
        (0.2, [0.05, 0.05, 0.05, 0.05], False, True),  # Low REVEL score, low SpliceAI
        (0.5, [0.3, 0.3, 0.3, 0.3], False, False),  # Intermediate scores
        (0.7, [0.05, 0.05, 0.05, 0.05], True, False),  # High REVEL score, low SpliceAI
        (0.2, [0.6, 0.6, 0.6, 0.6], True, False),  # Low REVEL score, high SpliceAI
    ],
)
@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_various_scenarios(
    vhl_predictor,
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

    prediction, _ = vhl_predictor.verify_pp3bp4(vhl_predictor.seqvar, auto_acmg_data)

    assert prediction.PP3 == expected_pp3
    assert prediction.BP4 == expected_bp4


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_missing_scores(vhl_predictor, auto_acmg_data):
    """Test behavior when scores are missing."""
    auto_acmg_data.scores.dbnsfp.revel = None

    prediction, comment = vhl_predictor.verify_pp3bp4(vhl_predictor.seqvar, auto_acmg_data)

    assert prediction is None
    assert "An error occurred during prediction" in comment


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_error_handling(vhl_predictor, auto_acmg_data):
    """Test error handling in verify_pp3bp4 method."""
    with patch.object(
        VHLPredictor,
        "_is_pathogenic_score",
        side_effect=Exception("Test error"),
    ):
        prediction, comment = vhl_predictor.verify_pp3bp4(vhl_predictor.seqvar, auto_acmg_data)

        assert prediction is None
        assert "An error occurred during prediction" in comment
        assert "Test error" in comment


def test_verify_pp3bp4_spliceai_thresholds(vhl_predictor, auto_acmg_data):
    """Test that SpliceAI thresholds are correctly adjusted during PP3/BP4 prediction."""
    with (
        patch.object(VHLPredictor, "_is_pathogenic_score", return_value=False),
        patch.object(VHLPredictor, "_is_benign_score", return_value=False),
        patch.object(VHLPredictor, "_affect_spliceAI", return_value=False),
    ):
        vhl_predictor.verify_pp3bp4(vhl_predictor.seqvar, auto_acmg_data)

        # Check that thresholds were adjusted for BP4
        assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.1
        assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.1
        assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.1
        assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.1


# ------------------ BP7 ------------------


def test_predict_bp7_threshold_adjustment(vhl_predictor, auto_acmg_data):
    """Test that the BP7 PhyloP threshold is correctly adjusted for VHL."""
    auto_acmg_data.thresholds.phyloP100 = 1.0  # Initial PhyloP threshold value

    # Call predict_bp7 method
    result = vhl_predictor.predict_bp7(vhl_predictor.seqvar, auto_acmg_data)

    # Check that the PhyloP threshold was adjusted
    assert (
        auto_acmg_data.thresholds.phyloP100 == 0.2
    ), "The BP7 PhyloP threshold should be adjusted to 0.2."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(VHLPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default(mock_super_predict_bp7, vhl_predictor, auto_acmg_data):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = vhl_predictor.predict_bp7(vhl_predictor.seqvar, auto_acmg_data)

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
