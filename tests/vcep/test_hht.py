from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import (
    PP3BP4,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
)
from src.defs.exceptions import AutoAcmgBaseException
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor
from src.vcep import HHTPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="12", pos=100, delete="A", insert="T")


@pytest.fixture
def hht_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return HHTPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


def test_predict_pm1_acvrl1_moderate(hht_predictor, auto_acmg_data):
    """Test when PM1 is met at the Moderate level for a variant in ACVRL1."""
    auto_acmg_data.hgnc_id = "HGNC:175"  # ACVRL1 gene
    auto_acmg_data.prot_pos = 213  # Within the glycine-rich loop (209-216)
    result = hht_predictor.predict_pm1(hht_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met at the Moderate level."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "critical residue for HGNC:175" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_eng_moderate(hht_predictor, auto_acmg_data):
    """Test when PM1 is met at the Moderate level for a variant in ENG."""
    auto_acmg_data.hgnc_id = "HGNC:3349"  # ENG gene
    auto_acmg_data.prot_pos = 278  # BMP9 binding site residue
    result = hht_predictor.predict_pm1(hht_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met at the Moderate level."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "critical residue for HGNC:3349" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_not_met(hht_predictor, auto_acmg_data):
    """Test when PM1 is not met for ACVRL1 but outside the critical regions."""
    auto_acmg_data.hgnc_id = "HGNC:175"  # ACVRL1 gene
    auto_acmg_data.prot_pos = 250  # Outside the critical regions
    result = hht_predictor.predict_pm1(hht_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met for ACVRL1."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not meet the PM1 criteria" in result.summary
    ), "The summary should indicate no criteria were met."


@patch("src.vcep.hht.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, hht_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method if logic changes."""
    auto_acmg_data.hgnc_id = "HGNC:111111"  # Gene not in the specific logic
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = hht_predictor.predict_pm1(hht_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


def test_predict_pm1_edge_case_start_boundary_acvrl1(hht_predictor, auto_acmg_data):
    """Test when variant falls exactly on the start boundary of a critical region in ACVRL1."""
    auto_acmg_data.hgnc_id = "HGNC:175"  # ACVRL1 gene
    auto_acmg_data.prot_pos = 209  # Start boundary of the glycine-rich loop (209-216)
    result = hht_predictor.predict_pm1(hht_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met on the start boundary."
    assert (
        "critical residue for HGNC:175" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_edge_case_end_boundary_eng(hht_predictor, auto_acmg_data):
    """Test when variant falls exactly on the end boundary of a critical region in ENG."""
    auto_acmg_data.hgnc_id = "HGNC:3349"  # ENG gene
    auto_acmg_data.prot_pos = 412  # End boundary of the critical cysteine residues
    result = hht_predictor.predict_pm1(hht_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met on the end boundary."
    assert (
        "critical residue for HGNC:3349" in result.summary
    ), "The summary should indicate the critical region."


def test_bs2_not_applicable(hht_predictor, auto_acmg_data):
    """Test when BS2 is not applicable for ACVRL1 and ENG."""
    result = hht_predictor._bs2_not_applicable(auto_acmg_data)
    assert result is True, "BS2 should always be not applicable for ACVRL1 and ENG."


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
def test_predict_pm2ba1bs1bs2(mock_super_method, hht_predictor, auto_acmg_data, seqvar):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = hht_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.00004
    assert auto_acmg_data.thresholds.ba1_benign == 0.01
    assert auto_acmg_data.thresholds.bs1_benign == 0.0008

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_bp3_not_applicable(hht_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = hht_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


def test_predict_pp2bp1(hht_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for HHT predictor."""

    # Call the method under test
    pp2_result, bp1_result = hht_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


def test_verify_pp3bp4_revel_thresholds(hht_predictor, auto_acmg_data):
    """Test that REVEL thresholds are correctly set for PP3/BP4 prediction."""
    with patch.object(HHTPredictor, "_is_missense", return_value=True):
        hht_predictor.verify_pp3bp4(hht_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.revel_pathogenic == 0.644
    assert auto_acmg_data.thresholds.revel_benign == 0.15


def test_verify_pp3bp4_spliceai_thresholds(hht_predictor, auto_acmg_data):
    """Test that SpliceAI thresholds are correctly set for PP3/BP4 prediction."""
    with patch.object(HHTPredictor, "_is_missense", return_value=True):
        hht_predictor.verify_pp3bp4(hht_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.01
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.01
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.01
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.01


@pytest.mark.parametrize(
    "is_missense, revel_score, spliceai_score, expected_pp3, expected_bp4",
    [
        (True, 0.7, 0.1, True, False),  # Missense, high REVEL score, low SpliceAI score
        (
            True,
            0.5,
            0.3,
            True,
            False,
        ),  # Missense, medium REVEL score, high SpliceAI score
        # (True, 0.1, 0.1, False, True),  # Missense, low REVEL score, low SpliceAI score
        (False, 0.7, 0.3, True, False),  # Non-missense, high SpliceAI score
        # (False, 0.7, 0.005, False, True),  # Non-missense, low SpliceAI score
    ],
)
def test_verify_pp3bp4_scenarios(
    hht_predictor,
    auto_acmg_data,
    is_missense,
    revel_score,
    spliceai_score,
    expected_pp3,
    expected_bp4,
):
    """Test different scenarios for PP3/BP4 prediction."""
    auto_acmg_data.scores.dbnsfp.revel = revel_score
    auto_acmg_data.scores.cadd.spliceAI_acceptor_gain = spliceai_score

    with (
        patch.object(HHTPredictor, "_is_missense", return_value=is_missense),
        patch.object(HHTPredictor, "_is_synonymous_variant", return_value=not is_missense),
        patch.object(HHTPredictor, "_is_intron_variant", return_value=not is_missense),
    ):
        prediction, comment = hht_predictor.verify_pp3bp4(hht_predictor.seqvar, auto_acmg_data)

        assert prediction.PP3 == expected_pp3
        assert prediction.BP4 == expected_bp4


def test_verify_pp3bp4_missense_revel_only(hht_predictor, auto_acmg_data):
    """Test PP3/BP4 prediction for missense variant when only REVEL score is available."""
    auto_acmg_data.scores.dbnsfp.revel = 0.7
    auto_acmg_data.scores.cadd.spliceAI_acceptor_gain = 0.1

    with patch.object(HHTPredictor, "_is_missense", return_value=True):
        prediction, comment = hht_predictor.verify_pp3bp4(hht_predictor.seqvar, auto_acmg_data)

        assert prediction.PP3 is True
        assert prediction.BP4 is False


def test_verify_pp3bp4_non_missense_spliceai_only(hht_predictor, auto_acmg_data):
    """Test PP3/BP4 prediction for non-missense variant when only SpliceAI score is available."""
    auto_acmg_data.scores.cadd.spliceAI_acceptor_gain = 0.3

    with (
        patch.object(HHTPredictor, "_is_missense", return_value=False),
        patch.object(HHTPredictor, "_is_synonymous_variant", return_value=True),
    ):
        prediction, comment = hht_predictor.verify_pp3bp4(hht_predictor.seqvar, auto_acmg_data)

        assert prediction.PP3 is True
        assert prediction.BP4 is False


def test_verify_pp3bp4_no_prediction(hht_predictor, auto_acmg_data):
    """Test PP3/BP4 prediction when no criteria are met."""
    auto_acmg_data.scores.dbnsfp.revel = 0.5
    auto_acmg_data.scores.cadd.spliceAI_acceptor_gain = 0.1

    with patch.object(HHTPredictor, "_is_missense", return_value=True):
        prediction, comment = hht_predictor.verify_pp3bp4(hht_predictor.seqvar, auto_acmg_data)

        assert prediction.PP3 is False
        assert prediction.BP4 is False


def test_verify_pp3bp4_error_handling(hht_predictor, auto_acmg_data):
    """Test error handling in verify_pp3bp4 method."""
    with patch.object(
        HHTPredictor, "_is_missense", side_effect=AutoAcmgBaseException("Test error")
    ):
        prediction, comment = hht_predictor.verify_pp3bp4(hht_predictor.seqvar, auto_acmg_data)

        assert prediction is None
        assert "An error occurred during prediction" in comment


def test_verify_pp3bp4_missing_data(hht_predictor, auto_acmg_data):
    """Test PP3/BP4 prediction when required data is missing."""
    auto_acmg_data.scores.dbnsfp.revel = None
    auto_acmg_data.scores.cadd.spliceAI_acceptor_gain = None

    with patch.object(HHTPredictor, "_is_missense", return_value=True):
        prediction, comment = hht_predictor.verify_pp3bp4(hht_predictor.seqvar, auto_acmg_data)

        assert prediction.PP3 is False
        assert prediction.BP4 is False


def test_verify_pp3bp4_return_type(hht_predictor, auto_acmg_data):
    """Test that verify_pp3bp4 returns the correct types."""
    prediction, comment = hht_predictor.verify_pp3bp4(hht_predictor.seqvar, auto_acmg_data)

    assert isinstance(prediction, PP3BP4)
    assert isinstance(comment, str)


def test_verify_pp3bp4_benign_spliceai_thresholds(hht_predictor, auto_acmg_data):
    """Test that benign SpliceAI thresholds are correctly set for PP3/BP4 prediction."""
    with patch.object(HHTPredictor, "_is_missense", return_value=True):
        hht_predictor.verify_pp3bp4(hht_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.01
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.01
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.01
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.01


def test_verify_pp3bp4_non_missense_thresholds(hht_predictor, auto_acmg_data):
    """Test that thresholds are correctly set for non-missense variants."""
    with (
        patch.object(HHTPredictor, "_is_missense", return_value=False),
        patch.object(HHTPredictor, "_is_synonymous_variant", return_value=True),
    ):
        hht_predictor.verify_pp3bp4(hht_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.01
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.01
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.01
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.01
