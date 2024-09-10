from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import (
    PS1PM5,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    GenomicStrand,
)
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor
from src.vcep import LeberCongenitalAmaurosisPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def leber_congenital_amaurosis_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return LeberCongenitalAmaurosisPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


# ---------------- PS1 & PM5 ----------------


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_overrides(
    mock_super_verify, leber_congenital_amaurosis_predictor, seqvar, auto_acmg_data
):
    """Test that the overridden verify_ps1pm5 method in LeberCongenitalAmaurosisPredictor works correctly."""
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
    prediction, comment = leber_congenital_amaurosis_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

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
    mock_super_verify, leber_congenital_amaurosis_predictor, seqvar, auto_acmg_data
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
    prediction, comment = leber_congenital_amaurosis_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

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
def test_verify_ps1pm5_non_missense(
    mock_super_verify, leber_congenital_amaurosis_predictor, seqvar, auto_acmg_data
):
    """Test that PS1 and PM5 remain as per superclass for non-missense variants."""
    # Set up the mock to return PS1 and PM5 as not applicable initially
    mock_super_verify.return_value = (
        PS1PM5(PS1=False, PM5=False),
        "Not applicable for non-missense",
    )

    # Setup the data with a non-missense variant
    auto_acmg_data.consequence = MagicMock(mehari=["synonymous_variant"])

    # Run the method under test
    prediction, comment = leber_congenital_amaurosis_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

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
    mock_super_verify, leber_congenital_amaurosis_predictor, seqvar, auto_acmg_data
):
    """Test verify_ps1pm5 method exception handling."""
    mock_super_verify.side_effect = Exception("Test exception")

    with pytest.raises(Exception) as exc_info:
        leber_congenital_amaurosis_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    assert "Test exception" in str(exc_info.value), "Should raise the original exception"


# ---------------- PM1 ----------------


def test_predict_pm1_met_for_critical_residue(leber_congenital_amaurosis_predictor, auto_acmg_data):
    """Test when the variant falls within the critical residues for RPE65."""
    auto_acmg_data.prot_pos = 180  # Critical residue for RPE65
    auto_acmg_data.hgnc_id = "HGNC:10294"  # RPE65 gene
    result = leber_congenital_amaurosis_predictor.predict_pm1(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for critical residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant affects a residue in RPE65 at position 180" in result.summary
    ), "The summary should indicate the affected residue."


def test_predict_pm1_met_for_residues_range(leber_congenital_amaurosis_predictor, auto_acmg_data):
    """Test when the variant falls within the residue range for RPE65."""
    auto_acmg_data.prot_pos = 110  # Within the residue range 107-125 for RPE65
    auto_acmg_data.hgnc_id = "HGNC:10294"  # RPE65 gene
    result = leber_congenital_amaurosis_predictor.predict_pm1(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for residues within the specified range."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant affects a residue in RPE65 at position 110" in result.summary
    ), "The summary should indicate the affected residue."


def test_predict_pm1_not_met(leber_congenital_amaurosis_predictor, auto_acmg_data):
    """Test when the variant does not fall within any critical residues or ranges."""
    auto_acmg_data.prot_pos = 300  # Not within any critical residues or ranges
    auto_acmg_data.hgnc_id = "HGNC:10294"  # RPE65 gene
    result = leber_congenital_amaurosis_predictor.predict_pm1(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
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


@patch("src.vcep.leber_congenital_amaurosis.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, leber_congenital_amaurosis_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method if logic changes."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the specific logic
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = leber_congenital_amaurosis_predictor.predict_pm1(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


# ---------------- PM2, BA1, BS1, BS2 ----------------


@patch.object(LeberCongenitalAmaurosisPredictor, "_get_af", return_value=0.1)
@patch.object(LeberCongenitalAmaurosisPredictor, "_ba1_exception", return_value=False)
def test_verify_pm2ba1bs1bs2(
    mock_get_af,
    mock_ba1_exception,
    leber_congenital_amaurosis_predictor,
    auto_acmg_data,
    seqvar,
):
    # Setup: Adjusting the thresholds to test under different conditions
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.02
    auto_acmg_data.thresholds.pm2_pathogenic = 0.001

    # Call the method under test
    result, comment = leber_congenital_amaurosis_predictor.verify_pm2ba1bs1bs2(
        seqvar, auto_acmg_data
    )

    # Assertions to validate the expected behavior
    assert result.BA1 is True, "Expected PM2 to be True based on the mocked allele frequency"

    # Assert changed thresholds
    assert (
        auto_acmg_data.thresholds.ba1_benign == 0.008
    ), "BA1 threshold should be adjusted to 0.008"
    assert (
        auto_acmg_data.thresholds.bs1_benign == 0.0008
    ), "BS1 threshold should be adjusted to 0.0008"


# ------------- PM4 & BP3 -------------


def test_bp3_not_applicable(leber_congenital_amaurosis_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = leber_congenital_amaurosis_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


# ------------- PP2 & BP1 -------------


def test_predict_pp2bp1(leber_congenital_amaurosis_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Leber Congenital Amaurosis predictor."""

    # Call the method under test
    pp2_result, bp1_result = leber_congenital_amaurosis_predictor.predict_pp2bp1(
        seqvar, auto_acmg_data
    )

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


# ------------- PP3 & BP4 -------------


def test_verify_pp3bp4_thresholds(leber_congenital_amaurosis_predictor, auto_acmg_data):
    """Test that the thresholds for PP3/BP4 prediction are correctly set."""
    leber_congenital_amaurosis_predictor.verify_pp3bp4(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
    )

    assert auto_acmg_data.thresholds.revel_pathogenic == 0.644
    assert auto_acmg_data.thresholds.revel_benign == 0.29
    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.1
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.1
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.1
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.1


@patch.object(LeberCongenitalAmaurosisPredictor, "_is_pathogenic_score")
@patch.object(LeberCongenitalAmaurosisPredictor, "_is_benign_score")
@patch.object(LeberCongenitalAmaurosisPredictor, "_affect_spliceAI")
def test_verify_pp3bp4_prediction_logic(
    mock_affect_spliceAI,
    mock_is_benign_score,
    mock_is_pathogenic_score,
    leber_congenital_amaurosis_predictor,
    auto_acmg_data,
):
    """Test the prediction logic for PP3 and BP4."""
    mock_is_pathogenic_score.return_value = True
    mock_is_benign_score.return_value = False
    mock_affect_spliceAI.side_effect = [
        True,
        False,
    ]  # First call True, second call False

    prediction, comment = leber_congenital_amaurosis_predictor.verify_pp3bp4(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 is True
    assert prediction.BP4 is False


@pytest.mark.parametrize(
    "revel_score, spliceAI_scores, expected_pp3, expected_bp4",
    [
        (0.7, [0.3, 0.3, 0.3, 0.3], True, False),  # High REVEL score, high SpliceAI
        (0.2, [0.1, 0.1, 0.1, 0.1], False, True),  # Low REVEL score, low SpliceAI
        (0.5, [0.15, 0.15, 0.15, 0.15], False, False),  # Intermediate scores
        (0.7, [0.1, 0.1, 0.1, 0.1], True, False),  # High REVEL score, low SpliceAI
        (0.2, [0.3, 0.3, 0.3, 0.3], True, False),  # Low REVEL score, high SpliceAI
    ],
)
def test_verify_pp3bp4_various_scenarios(
    leber_congenital_amaurosis_predictor,
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

    prediction, _ = leber_congenital_amaurosis_predictor.verify_pp3bp4(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 == expected_pp3
    assert prediction.BP4 == expected_bp4


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_missing_scores(leber_congenital_amaurosis_predictor, auto_acmg_data):
    """Test behavior when scores are missing."""
    auto_acmg_data.scores.dbnsfp.revel = None

    prediction, comment = leber_congenital_amaurosis_predictor.verify_pp3bp4(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
    )

    assert prediction is None
    assert "An error occurred during prediction" in comment


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_error_handling(leber_congenital_amaurosis_predictor, auto_acmg_data):
    """Test error handling in verify_pp3bp4 method."""
    with patch.object(
        LeberCongenitalAmaurosisPredictor,
        "_is_pathogenic_score",
        side_effect=Exception("Test error"),
    ):
        prediction, comment = leber_congenital_amaurosis_predictor.verify_pp3bp4(
            leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
        )

        assert prediction is None
        assert "An error occurred during prediction" in comment
        assert "Test error" in comment


def test_verify_pp3bp4_spliceai_thresholds(leber_congenital_amaurosis_predictor, auto_acmg_data):
    """Test that SpliceAI thresholds are correctly adjusted during PP3/BP4 prediction."""
    with (
        patch.object(
            LeberCongenitalAmaurosisPredictor,
            "_is_pathogenic_score",
            return_value=False,
        ),
        patch.object(LeberCongenitalAmaurosisPredictor, "_is_benign_score", return_value=False),
        patch.object(LeberCongenitalAmaurosisPredictor, "_affect_spliceAI", return_value=False),
    ):
        leber_congenital_amaurosis_predictor.verify_pp3bp4(
            leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
        )

        # Check that thresholds were adjusted for BP4
        assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.1
        assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.1
        assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.1
        assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.1


# ------------- BP7 -------------


def test_is_bp7_exception_first_base_of_exon(leber_congenital_amaurosis_predictor, auto_acmg_data):
    """Test that the BP7 exception is detected for a synonymous variant at the first base of an exon."""
    auto_acmg_data.exons = [
        MagicMock(altStartI=100, altEndI=200)
    ]  # Exon with start at position 100
    auto_acmg_data.strand = GenomicStrand.Plus

    assert leber_congenital_amaurosis_predictor._is_bp7_exception(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
    ), "The variant should be detected as a BP7 exception at the first base of an exon."


def test_is_bp7_exception_last_three_bases_of_exon(
    leber_congenital_amaurosis_predictor, auto_acmg_data
):
    """Test that the BP7 exception is detected for a synonymous variant in the last 3 bases of an exon."""
    auto_acmg_data.exons = [MagicMock(altStartI=100, altEndI=200)]  # Exon with end at position 200
    auto_acmg_data.strand = GenomicStrand.Plus
    leber_congenital_amaurosis_predictor.seqvar.pos = 198  # Position within the last 3 bases

    assert leber_congenital_amaurosis_predictor._is_bp7_exception(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
    ), "The variant should be detected as a BP7 exception in the last 3 bases of an exon."


def test_is_bp7_exception_no_exception(leber_congenital_amaurosis_predictor, auto_acmg_data):
    """Test that no BP7 exception is detected when the variant is outside the first base or last 3 bases of an exon."""
    auto_acmg_data.exons = [MagicMock(altStartI=100, altEndI=200)]  # Exon with end at position 200
    auto_acmg_data.strand = GenomicStrand.Plus
    leber_congenital_amaurosis_predictor.seqvar.pos = 150  # Position not at the boundary

    assert not leber_congenital_amaurosis_predictor._is_bp7_exception(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
    ), "The variant should not be detected as a BP7 exception."


def test_predict_bp7_threshold_adjustment(leber_congenital_amaurosis_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted for RPE65."""
    auto_acmg_data.thresholds.bp7_donor = 1  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = leber_congenital_amaurosis_predictor.predict_bp7(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
    )

    # Check that the thresholds were adjusted for RPE65
    assert (
        auto_acmg_data.thresholds.bp7_donor == 7
    ), "The BP7 donor threshold should be adjusted to 7."
    assert (
        auto_acmg_data.thresholds.bp7_acceptor == 21
    ), "The BP7 acceptor threshold should be adjusted to 21."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(LeberCongenitalAmaurosisPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default(
    mock_super_predict_bp7, leber_congenital_amaurosis_predictor, auto_acmg_data
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
    result = leber_congenital_amaurosis_predictor.predict_bp7(
        leber_congenital_amaurosis_predictor.seqvar, auto_acmg_data
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
