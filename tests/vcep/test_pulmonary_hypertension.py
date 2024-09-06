from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    GenomicStrand,
)
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor
from src.vcep.pulmonary_hypertension import PulmonaryHypertensionPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="2", pos=100, delete="A", insert="T")


@pytest.fixture
def pulmonary_hypertension_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return PulmonaryHypertensionPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


def test_predict_pm1_strong_criteria(pulmonary_hypertension_predictor, auto_acmg_data):
    """Test when the variant falls within strong criteria for BMPR2."""
    auto_acmg_data.hgnc_id = "HGNC:1078"  # BMPR2 gene
    auto_acmg_data.prot_pos = 210  # Within the strong criteria
    result = pulmonary_hypertension_predictor.predict_pm1(
        pulmonary_hypertension_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.Met, "PM1 should be met at the Strong level."
    assert (
        result.strength == AutoACMGStrength.PathogenicStrong
    ), "The strength should be PathogenicStrong."
    assert (
        "critical residue in BMPR2" in result.summary
    ), "The summary should indicate a critical residue."


def test_predict_pm1_moderate_criteria(pulmonary_hypertension_predictor, auto_acmg_data):
    """Test when the variant falls within moderate criteria for BMPR2."""
    auto_acmg_data.hgnc_id = "HGNC:1078"  # BMPR2 gene
    auto_acmg_data.prot_pos = 250  # Within the moderate criteria
    result = pulmonary_hypertension_predictor.predict_pm1(
        pulmonary_hypertension_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.Met, "PM1 should be met at the Moderate level."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "within the extracellular or kinase domain" in result.summary
    ), "The summary should indicate the domain."


def test_predict_pm1_non_critical_residue(pulmonary_hypertension_predictor, auto_acmg_data):
    """Test when the variant affects a non-critical residue in BMPR2."""
    auto_acmg_data.hgnc_id = "HGNC:1078"  # BMPR2 gene
    auto_acmg_data.prot_pos = 42  # Non-critical residue
    result = pulmonary_hypertension_predictor.predict_pm1(
        pulmonary_hypertension_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for non-critical residues."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "demonstrated to be non-critical for kinase activity" in result.summary
    ), "The summary should indicate non-critical residue."


@patch("src.vcep.pulmonary_hypertension.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, pulmonary_hypertension_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the BMPR2 VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = pulmonary_hypertension_predictor.predict_pm1(
        pulmonary_hypertension_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
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
    mock_super_method, pulmonary_hypertension_predictor, auto_acmg_data, seqvar
):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = pulmonary_hypertension_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.0001
    assert auto_acmg_data.thresholds.ba1_benign == 0.01
    assert auto_acmg_data.thresholds.bs1_benign == 0.001

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_predict_pp2bp1(pulmonary_hypertension_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Pulmonary Hypertension."""

    # Call the method under test
    pp2_result, bp1_result = pulmonary_hypertension_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


def test_is_bp7_exception_first_base_of_exon(pulmonary_hypertension_predictor, auto_acmg_data):
    """Test that the BP7 exception is detected for a synonymous variant at the first base of an exon."""
    auto_acmg_data.exons = [
        MagicMock(altStartI=100, altEndI=200)
    ]  # Exon with start at position 100
    auto_acmg_data.strand = GenomicStrand.Plus

    assert pulmonary_hypertension_predictor._is_bp7_exception(
        pulmonary_hypertension_predictor.seqvar, auto_acmg_data
    ), "The variant should be detected as a BP7 exception at the first base of an exon."


def test_is_bp7_exception_last_three_bases_of_exon(
    pulmonary_hypertension_predictor, auto_acmg_data
):
    """Test that the BP7 exception is detected for a synonymous variant in the last 3 bases of an exon."""
    auto_acmg_data.exons = [MagicMock(altStartI=100, altEndI=200)]  # Exon with end at position 200
    auto_acmg_data.strand = GenomicStrand.Plus
    pulmonary_hypertension_predictor.seqvar.pos = 198  # Position within the last 3 bases

    assert pulmonary_hypertension_predictor._is_bp7_exception(
        pulmonary_hypertension_predictor.seqvar, auto_acmg_data
    ), "The variant should be detected as a BP7 exception in the last 3 bases of an exon."


def test_is_bp7_exception_no_exception(pulmonary_hypertension_predictor, auto_acmg_data):
    """Test that no BP7 exception is detected when the variant is outside the first base or last 3 bases of an exon."""
    auto_acmg_data.exons = [MagicMock(altStartI=100, altEndI=200)]  # Exon with end at position 200
    auto_acmg_data.strand = GenomicStrand.Plus
    pulmonary_hypertension_predictor.seqvar.pos = 150  # Position not at the boundary

    assert not pulmonary_hypertension_predictor._is_bp7_exception(
        pulmonary_hypertension_predictor.seqvar, auto_acmg_data
    ), "The variant should not be detected as a BP7 exception."


def test_verify_pp3bp4_thresholds(pulmonary_hypertension_predictor, auto_acmg_data):
    """Test that the thresholds for PP3/BP4 prediction are correctly set."""
    pulmonary_hypertension_predictor.verify_pp3bp4(
        pulmonary_hypertension_predictor.seqvar, auto_acmg_data
    )

    assert auto_acmg_data.thresholds.revel_pathogenic == 0.75
    assert auto_acmg_data.thresholds.revel_benign == 0.25
    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.1
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.1
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.1
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.1


@patch.object(PulmonaryHypertensionPredictor, "_is_pathogenic_score")
@patch.object(PulmonaryHypertensionPredictor, "_is_benign_score")
@patch.object(PulmonaryHypertensionPredictor, "_affect_spliceAI")
def test_verify_pp3bp4_prediction_logic(
    mock_affect_spliceAI,
    mock_is_benign_score,
    mock_is_pathogenic_score,
    pulmonary_hypertension_predictor,
    auto_acmg_data,
):
    """Test the prediction logic for PP3 and BP4."""
    mock_is_pathogenic_score.return_value = True
    mock_is_benign_score.return_value = False
    mock_affect_spliceAI.side_effect = [True, False]  # First call True, second call False

    prediction, comment = pulmonary_hypertension_predictor.verify_pp3bp4(
        pulmonary_hypertension_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 is True
    assert prediction.BP4 is False


@pytest.mark.parametrize(
    "revel_score, cadd_score, spliceAI_scores, expected_pp3, expected_bp4",
    [
        (0.8, 25, [0.3, 0.3, 0.3, 0.3], True, False),  # High REVEL and CADD scores, high SpliceAI
        (0.2, 10, [0.05, 0.05, 0.05, 0.05], False, True),  # Low REVEL and CADD scores, low SpliceAI
        (0.5, 20, [0.15, 0.15, 0.15, 0.15], False, False),  # Intermediate scores
        (0.8, 15, [0.05, 0.05, 0.05, 0.05], True, False),  # High REVEL, low CADD and SpliceAI
        (0.2, 25, [0.3, 0.3, 0.3, 0.3], True, False),  # Low REVEL, high CADD and SpliceAI
    ],
)
def test_verify_pp3bp4_various_scenarios(
    pulmonary_hypertension_predictor,
    auto_acmg_data,
    revel_score,
    cadd_score,
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

    prediction, _ = pulmonary_hypertension_predictor.verify_pp3bp4(
        pulmonary_hypertension_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 == expected_pp3
    assert prediction.BP4 == expected_bp4


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_missing_scores(pulmonary_hypertension_predictor, auto_acmg_data):
    """Test behavior when scores are missing."""
    auto_acmg_data.scores.dbnsfp.revel = None
    auto_acmg_data.scores.cadd.phred = None

    prediction, comment = pulmonary_hypertension_predictor.verify_pp3bp4(
        pulmonary_hypertension_predictor.seqvar, auto_acmg_data
    )

    assert prediction is None
    assert "An error occurred during prediction" in comment


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_error_handling(pulmonary_hypertension_predictor, auto_acmg_data):
    """Test error handling in verify_pp3bp4 method."""
    with patch.object(
        PulmonaryHypertensionPredictor,
        "_is_pathogenic_score",
        side_effect=Exception("Test error"),
    ):
        prediction, comment = pulmonary_hypertension_predictor.verify_pp3bp4(
            pulmonary_hypertension_predictor.seqvar, auto_acmg_data
        )

        assert prediction is None
        assert "An error occurred during prediction" in comment
        assert "Test error" in comment


def test_verify_pp3bp4_spliceai_thresholds(pulmonary_hypertension_predictor, auto_acmg_data):
    """Test that SpliceAI thresholds are correctly adjusted during PP3/BP4 prediction."""
    with (
        patch.object(PulmonaryHypertensionPredictor, "_is_pathogenic_score", return_value=False),
        patch.object(PulmonaryHypertensionPredictor, "_is_benign_score", return_value=False),
        patch.object(PulmonaryHypertensionPredictor, "_affect_spliceAI", return_value=False),
    ):

        pulmonary_hypertension_predictor.verify_pp3bp4(
            pulmonary_hypertension_predictor.seqvar, auto_acmg_data
        )

        # Check that thresholds were adjusted for BP4
        assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.1
        assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.1
        assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.1
        assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.1
