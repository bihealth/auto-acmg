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
from src.vcep import LysosomalDiseasesPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="17", pos=100, delete="A", insert="T")


@pytest.fixture
def lysosomal_diseases_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return LysosomalDiseasesPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


def test_predict_pm1_met_for_critical_residue(lysosomal_diseases_predictor, auto_acmg_data):
    """Test when the variant falls within the critical residues for GAA."""
    auto_acmg_data.prot_pos = 282  # Critical residue for GAA
    auto_acmg_data.hgnc_id = "HGNC:4065"  # GAA gene
    result = lysosomal_diseases_predictor.predict_pm1(
        lysosomal_diseases_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant affects a residue in GAA at position 282" in result.summary
    ), "The summary should indicate the affected residue."


def test_predict_pm1_not_met(lysosomal_diseases_predictor, auto_acmg_data):
    """Test when the variant does not fall within any critical residues for GAA."""
    auto_acmg_data.prot_pos = 700  # Not within any critical residues
    auto_acmg_data.hgnc_id = "HGNC:4065"  # GAA gene
    result = lysosomal_diseases_predictor.predict_pm1(
        lysosomal_diseases_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met for non-critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant does not meet the PM1 criteria" in result.summary
    ), "The summary should indicate the lack of criteria met."


@patch("src.vcep.lysosomal_diseases.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, lysosomal_diseases_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method if logic changes."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the specific logic
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = lysosomal_diseases_predictor.predict_pm1(
        lysosomal_diseases_predictor.seqvar, auto_acmg_data
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
    mock_super_method, lysosomal_diseases_predictor, auto_acmg_data, seqvar
):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = lysosomal_diseases_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.001
    assert auto_acmg_data.thresholds.ba1_benign == 0.01
    assert auto_acmg_data.thresholds.bs1_benign == 0.005

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_bp3_not_applicable(lysosomal_diseases_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = lysosomal_diseases_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


def test_predict_pp2bp1(lysosomal_diseases_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Lysosomal Diseases predictor."""

    # Call the method under test
    pp2_result, bp1_result = lysosomal_diseases_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


def test_verify_pp3bp4_thresholds(lysosomal_diseases_predictor, auto_acmg_data):
    """Test that the thresholds for PP3/BP4 prediction are correctly set."""
    lysosomal_diseases_predictor.verify_pp3bp4(lysosomal_diseases_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.revel_pathogenic == 0.7
    assert auto_acmg_data.thresholds.revel_benign == 0.5
    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2


@patch.object(LysosomalDiseasesPredictor, "_is_pathogenic_score")
@patch.object(LysosomalDiseasesPredictor, "_is_benign_score")
@patch.object(LysosomalDiseasesPredictor, "_affect_spliceAI")
def test_verify_pp3bp4_prediction_logic(
    mock_affect_spliceAI,
    mock_is_benign_score,
    mock_is_pathogenic_score,
    lysosomal_diseases_predictor,
    auto_acmg_data,
):
    """Test the prediction logic for PP3 and BP4."""
    mock_is_pathogenic_score.return_value = True
    mock_is_benign_score.return_value = False
    mock_affect_spliceAI.side_effect = [
        True,
        False,
    ]  # First call True, second call False

    prediction, comment = lysosomal_diseases_predictor.verify_pp3bp4(
        lysosomal_diseases_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 is True
    assert prediction.BP4 is False


@pytest.mark.parametrize(
    "revel_score, spliceAI_scores, expected_pp3, expected_bp4",
    [
        (0.8, [0.6, 0.6, 0.6, 0.6], True, False),  # High REVEL score, high SpliceAI
        (0.4, [0.1, 0.1, 0.1, 0.1], False, True),  # Low REVEL score, low SpliceAI
        (0.6, [0.3, 0.3, 0.3, 0.3], False, False),  # Intermediate scores
        (0.8, [0.1, 0.1, 0.1, 0.1], True, False),  # High REVEL score, low SpliceAI
        (0.4, [0.6, 0.6, 0.6, 0.6], True, False),  # Low REVEL score, high SpliceAI
    ],
)
def test_verify_pp3bp4_various_scenarios(
    lysosomal_diseases_predictor,
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

    prediction, _ = lysosomal_diseases_predictor.verify_pp3bp4(
        lysosomal_diseases_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 == expected_pp3
    assert prediction.BP4 == expected_bp4


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_missing_scores(lysosomal_diseases_predictor, auto_acmg_data):
    """Test behavior when scores are missing."""
    auto_acmg_data.scores.dbnsfp.revel = None

    prediction, comment = lysosomal_diseases_predictor.verify_pp3bp4(
        lysosomal_diseases_predictor.seqvar, auto_acmg_data
    )

    assert prediction is None
    assert "An error occurred during prediction" in comment


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_error_handling(lysosomal_diseases_predictor, auto_acmg_data):
    """Test error handling in verify_pp3bp4 method."""
    with patch.object(
        LysosomalDiseasesPredictor,
        "_is_pathogenic_score",
        side_effect=Exception("Test error"),
    ):
        prediction, comment = lysosomal_diseases_predictor.verify_pp3bp4(
            lysosomal_diseases_predictor.seqvar, auto_acmg_data
        )

        assert prediction is None
        assert "An error occurred during prediction" in comment
        assert "Test error" in comment


def test_verify_pp3bp4_spliceai_thresholds(lysosomal_diseases_predictor, auto_acmg_data):
    """Test that SpliceAI thresholds are correctly adjusted during PP3/BP4 prediction."""
    with (
        patch.object(LysosomalDiseasesPredictor, "_is_pathogenic_score", return_value=False),
        patch.object(LysosomalDiseasesPredictor, "_is_benign_score", return_value=False),
        patch.object(LysosomalDiseasesPredictor, "_affect_spliceAI", return_value=False),
    ):
        lysosomal_diseases_predictor.verify_pp3bp4(
            lysosomal_diseases_predictor.seqvar, auto_acmg_data
        )

        # Check that thresholds were adjusted for BP4
        assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
        assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
        assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
        assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2
