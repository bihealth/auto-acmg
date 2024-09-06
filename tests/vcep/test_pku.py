from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import (
    AlleleCondition,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    GenomicStrand,
)
from src.defs.exceptions import MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor
from src.vcep.pku import PKUPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="12", pos=100, delete="A", insert="T")


@pytest.fixture
def pku_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return PKUPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


def test_predict_pm1_residue_critical_for_pku(pku_predictor, auto_acmg_data):
    """Test when the variant affects a residue critical for PKU."""
    auto_acmg_data.hgnc_id = "HGNC:8582"  # PAH gene
    auto_acmg_data.prot_pos = 138  # Tyr138, an active site residue
    result = pku_predictor.predict_pm1(pku_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.Met, "PM1 should be met for critical residues."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert "critical for PKU" in result.summary, "The summary should indicate the critical residue."


def test_predict_pm1_residue_not_critical_for_pku(pku_predictor, auto_acmg_data):
    """Test when the variant does not affect a residue critical for PKU."""
    auto_acmg_data.hgnc_id = "HGNC:8582"  # PAH gene
    auto_acmg_data.prot_pos = 500  # Position outside the critical residues
    result = pku_predictor.predict_pm1(pku_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for non-critical residues."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not meet the PM1 criteria" in result.summary
    ), "The summary should indicate that criteria were not met."


@patch("src.vcep.pku.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, pku_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the PM1_CLUSTER_PKU
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = pku_predictor.predict_pm1(pku_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


@patch.object(
    PKUPredictor,
    "_get_allele_cond",
    return_value=AlleleCondition.Recessive,
)
@patch.object(
    PKUPredictor,
    "_get_control_af",
    return_value=MagicMock(bySex=MagicMock(overall=MagicMock(ac=10, nhomalt=6))),
)
@patch.object(PKUPredictor, "_get_any_af", return_value=None)
def test_check_zyg_homozygous_positive(
    mock_get_any_af,
    mock_get_control_af,
    mock_get_allele_cond,
    pku_predictor,
    seqvar,
    auto_acmg_data,
):
    pku_predictor.comment_pm2ba1bs1bs2 = ""
    assert pku_predictor._check_zyg(seqvar, auto_acmg_data) == True
    assert (
        "The variant is in a recessive (homozygous) disorder." in pku_predictor.comment_pm2ba1bs1bs2
    )


@patch.object(
    PKUPredictor,
    "_get_allele_cond",
    return_value=AlleleCondition.Recessive,
)
@patch.object(
    PKUPredictor,
    "_get_control_af",
    return_value=MagicMock(bySex=MagicMock(overall=MagicMock(ac=10, nhomalt=4))),
)
@patch.object(PKUPredictor, "_get_any_af", return_value=None)
def test_check_zyg_homozygous_negative(
    mock_get_any_af,
    mock_get_control_af,
    mock_get_allele_cond,
    pku_predictor,
    seqvar,
    auto_acmg_data,
):
    pku_predictor.comment_pm2ba1bs1bs2 = ""
    assert pku_predictor._check_zyg(seqvar, auto_acmg_data) == False


@patch.object(
    PKUPredictor,
    "_get_allele_cond",
    return_value=AlleleCondition.Recessive,
)
@patch.object(PKUPredictor, "_get_control_af", return_value=None)
@patch.object(PKUPredictor, "_get_any_af", return_value=None)
def test_check_zyg_missing_data_raises_error(
    mock_get_any_af,
    mock_get_control_af,
    mock_get_allele_cond,
    pku_predictor,
    seqvar,
    auto_acmg_data,
):
    pku_predictor.comment_pm2ba1bs1bs2 = ""
    with pytest.raises(MissingDataError):
        pku_predictor._check_zyg(seqvar, auto_acmg_data)


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
def test_predict_pm2ba1bs1bs2(mock_super_method, pku_predictor, auto_acmg_data, seqvar):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = pku_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.000002
    assert auto_acmg_data.thresholds.ba1_benign == 0.015
    assert auto_acmg_data.thresholds.bs1_benign == 0.002

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_predict_pp2bp1(pku_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for PKU."""

    # Call the method under test
    pp2_result, bp1_result = pku_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


def test_is_bp7_exception_first_base_of_exon(pku_predictor, auto_acmg_data):
    """Test that the BP7 exception is detected for a synonymous variant at the first base of an exon."""
    auto_acmg_data.exons = [
        MagicMock(altStartI=100, altEndI=200)
    ]  # Exon with start at position 100
    auto_acmg_data.strand = GenomicStrand.Plus

    assert pku_predictor._is_bp7_exception(
        pku_predictor.seqvar, auto_acmg_data
    ), "The variant should be detected as a BP7 exception at the first base of an exon."


def test_is_bp7_exception_last_three_bases_of_exon(pku_predictor, auto_acmg_data):
    """Test that the BP7 exception is detected for a synonymous variant in the last 3 bases of an exon."""
    auto_acmg_data.exons = [MagicMock(altStartI=100, altEndI=200)]  # Exon with end at position 200
    auto_acmg_data.strand = GenomicStrand.Plus
    pku_predictor.seqvar.pos = 198  # Position within the last 3 bases

    assert pku_predictor._is_bp7_exception(
        pku_predictor.seqvar, auto_acmg_data
    ), "The variant should be detected as a BP7 exception in the last 3 bases of an exon."


def test_is_bp7_exception_no_exception(pku_predictor, auto_acmg_data):
    """Test that no BP7 exception is detected when the variant is outside the first base or last 3 bases of an exon."""
    auto_acmg_data.exons = [MagicMock(altStartI=100, altEndI=200)]  # Exon with end at position 200
    auto_acmg_data.strand = GenomicStrand.Plus
    pku_predictor.seqvar.pos = 150  # Position not at the boundary

    assert not pku_predictor._is_bp7_exception(
        pku_predictor.seqvar, auto_acmg_data
    ), "The variant should not be detected as a BP7 exception."


def test_predict_bp7_threshold_adjustment(pku_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted."""
    auto_acmg_data.thresholds.bp7_donor = 1  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = pku_predictor.predict_bp7(pku_predictor.seqvar, auto_acmg_data)

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
def test_predict_bp7_fallback_to_default(mock_super_predict_bp7, pku_predictor, auto_acmg_data):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = pku_predictor.predict_bp7(pku_predictor.seqvar, auto_acmg_data)

    # Verify the result and ensure the superclass method was called
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.NotMet, "BP7 should return NotMet as mocked."
    assert (
        result.strength == AutoACMGStrength.BenignSupporting
    ), "The strength should be BenignSupporting."
    assert (
        "Default BP7 prediction fallback." in result.summary
    ), "The summary should indicate the fallback."
    assert mock_super_predict_bp7.called, "super().predict_bp7 should have been called."


def test_verify_pp3bp4_thresholds(pku_predictor, auto_acmg_data):
    """Test that the thresholds for PP3/BP4 prediction are correctly set."""
    pku_predictor.verify_pp3bp4(pku_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.revel_pathogenic == 0.644
    assert auto_acmg_data.thresholds.revel_benign == 0.29
    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.1
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.1
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.1
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.1


@patch.object(PKUPredictor, "_is_pathogenic_score")
@patch.object(PKUPredictor, "_is_benign_score")
@patch.object(PKUPredictor, "_affect_spliceAI")
def test_verify_pp3bp4_prediction_logic(
    mock_affect_spliceAI,
    mock_is_benign_score,
    mock_is_pathogenic_score,
    pku_predictor,
    auto_acmg_data,
):
    """Test the prediction logic for PP3 and BP4."""
    mock_is_pathogenic_score.return_value = True
    mock_is_benign_score.return_value = False
    mock_affect_spliceAI.side_effect = [True, False]  # First call True, second call False

    prediction, comment = pku_predictor.verify_pp3bp4(pku_predictor.seqvar, auto_acmg_data)

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
def test_verify_pp3bp4_various_scenarios(
    pku_predictor,
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

    prediction, _ = pku_predictor.verify_pp3bp4(pku_predictor.seqvar, auto_acmg_data)

    assert prediction.PP3 == expected_pp3
    assert prediction.BP4 == expected_bp4


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_missing_scores(pku_predictor, auto_acmg_data):
    """Test behavior when scores are missing."""
    auto_acmg_data.scores.dbnsfp.revel = None

    prediction, comment = pku_predictor.verify_pp3bp4(pku_predictor.seqvar, auto_acmg_data)

    assert prediction is None
    assert "An error occurred during prediction" in comment


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_error_handling(pku_predictor, auto_acmg_data):
    """Test error handling in verify_pp3bp4 method."""
    with patch.object(
        PKUPredictor,
        "_is_pathogenic_score",
        side_effect=Exception("Test error"),
    ):
        prediction, comment = pku_predictor.verify_pp3bp4(pku_predictor.seqvar, auto_acmg_data)

        assert prediction is None
        assert "An error occurred during prediction" in comment
        assert "Test error" in comment


def test_verify_pp3bp4_spliceai_thresholds(pku_predictor, auto_acmg_data):
    """Test that SpliceAI thresholds are correctly adjusted during PP3/BP4 prediction."""
    with (
        patch.object(PKUPredictor, "_is_pathogenic_score", return_value=False),
        patch.object(PKUPredictor, "_is_benign_score", return_value=False),
        patch.object(PKUPredictor, "_affect_spliceAI", return_value=False),
    ):

        pku_predictor.verify_pp3bp4(pku_predictor.seqvar, auto_acmg_data)

        # Check that thresholds were adjusted for BP4
        assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.1
        assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.1
        assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.1
        assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.1
