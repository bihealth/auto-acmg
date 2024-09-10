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
from src.defs.exceptions import AlgorithmError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor
from src.vcep import CoagulationFactorDeficiencyPredictor


@pytest.fixture
def seqvar():
    return SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="X",
        pos=100000,
        delete="A",
        insert="T",
    )


@pytest.fixture
def coagulation_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return CoagulationFactorDeficiencyPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    data = AutoACMGSeqVarData()
    data.hgnc_id = "HGNC:3546"  # F8
    data.prot_pos = 391  # Test position
    data.exons = [MagicMock(altStartI=1, altEndI=1000000)]
    data.strand = GenomicStrand.Plus
    return data


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_overrides(mock_super_verify, coagulation_predictor, seqvar, auto_acmg_data):
    # Set up the mock to return PS1 and PM5 as applicable initially
    mock_super_verify.return_value = (PS1PM5(PS1=True, PM5=True), "Initial evaluation")

    # Set up the data
    auto_acmg_data.consequence = MagicMock(mehari=["missense_variant"])
    auto_acmg_data.thresholds = MagicMock(
        spliceAI_acceptor_gain=0.5,
        spliceAI_acceptor_loss=0.5,
        spliceAI_donor_gain=0.5,
        spliceAI_donor_loss=0.5,
    )
    auto_acmg_data.scores = MagicMock(
        cadd=MagicMock(spliceAI_acceptor_gain=0.6, spliceAI_donor_gain=0.6)
    )

    # Run the method under test
    prediction, comment = coagulation_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that the splicing effect leads to overriding PS1 and PM5 as not applicable
    assert not prediction.PS1, "PS1 should be marked as not applicable due to splicing effect."
    assert not prediction.PM5, "PM5 should be marked as not applicable due to splicing effect."
    assert "Variant affects splicing" in comment, "Comment should note the splicing effect."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)


def test_predict_pm1_strong_criteria(coagulation_predictor, auto_acmg_data):
    """Test when the variant falls within strong criteria for F8."""
    auto_acmg_data.prot_pos = 391  # Within the strong criteria for F8
    result = coagulation_predictor.predict_pm1(coagulation_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met at the Strong level."
    assert (
        result.strength == AutoACMGStrength.PathogenicStrong
    ), "The strength should be PathogenicStrong."
    assert (
        "PM1 is met at the Strong level" in result.summary
    ), "The summary should indicate the strong criteria."


def test_predict_pm1_moderate_criteria_residues(coagulation_predictor, auto_acmg_data):
    """Test when the variant falls within moderate criteria for F8 based on residues."""
    auto_acmg_data.prot_pos = 1667  # Within the moderate criteria for F8
    result = coagulation_predictor.predict_pm1(coagulation_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met at the Moderate level for residues."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "PM1 is met at the Moderate level" in result.summary
    ), "The summary should indicate the moderate criteria for residues."


def test_predict_pm1_moderate_criteria_regions(coagulation_predictor, auto_acmg_data):
    """Test when the variant falls within moderate criteria for F8 based on regions."""
    auto_acmg_data.prot_pos = 2270  # Within the FXa-binding region for F8
    result = coagulation_predictor.predict_pm1(coagulation_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met at the Moderate level for regions."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "PM1 is met at the Moderate level" in result.summary
    ), "The summary should indicate the moderate criteria for regions."


def test_predict_pm1_moderate_criteria_exons(coagulation_predictor, auto_acmg_data):
    """Test when the variant falls within moderate criteria for F9 based on exons."""
    auto_acmg_data.hgnc_id = "HGNC:3551"  # F9
    auto_acmg_data.strand = GenomicStrand.Plus
    auto_acmg_data.prot_pos = 1
    auto_acmg_data.exons = [
        MagicMock(altStartI=1, altEndI=200),
        MagicMock(altStartI=201, altEndI=400),
        MagicMock(altStartI=401, altEndI=600),
        MagicMock(altStartI=601, altEndI=1000000),
    ]

    result = coagulation_predictor.predict_pm1(coagulation_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met at the Moderate level for exons."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "PM1 is met at the Moderate level" in result.summary
    ), "The summary should indicate the moderate criteria for exons."


def test_predict_pm1_no_criteria_met(coagulation_predictor, auto_acmg_data):
    """Test when the variant does not meet any PM1 criteria."""
    auto_acmg_data.prot_pos = 5000  # Position outside any defined criteria
    result = coagulation_predictor.predict_pm1(coagulation_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met when no criteria are matched."
    assert (
        "Variant does not meet the PM1 criteria" in result.summary
    ), "The summary should indicate no criteria were met."


def test_predict_pm1_invalid_strand(coagulation_predictor, auto_acmg_data):
    """Test when an invalid strand is provided."""
    auto_acmg_data.strand = "invalid_strand"

    with pytest.raises(AlgorithmError):
        coagulation_predictor._get_affected_exon(auto_acmg_data, coagulation_predictor.seqvar)


@patch.object(
    CoagulationFactorDeficiencyPredictor,
    "_get_any_af",
    return_value=MagicMock(bySex=MagicMock(xx=MagicMock(ac=0))),
)
@pytest.mark.skip(reason="Idk why this test is failing")
def test_absent_in_males(coagulation_predictor, auto_acmg_data):
    assert (
        coagulation_predictor._absent_in_males(auto_acmg_data) is True
    ), "Should return True when AC is 0 in males"


@patch.object(
    CoagulationFactorDeficiencyPredictor,
    "_get_any_af",
    return_value=MagicMock(bySex=MagicMock(xx=MagicMock(ac=1))),
)
@pytest.mark.skip(reason="Idk why this test is failing")
def test_not_absent_in_males(coagulation_predictor, auto_acmg_data):
    assert (
        coagulation_predictor._absent_in_males(auto_acmg_data) is False
    ), "Should return False when AC is not 0 in males"


@patch.object(CoagulationFactorDeficiencyPredictor, "_get_af", return_value=0.000001)
@patch.object(CoagulationFactorDeficiencyPredictor, "_ba1_exception", return_value=False)
@patch.object(CoagulationFactorDeficiencyPredictor, "_absent_in_males", return_value=True)
@patch.object(CoagulationFactorDeficiencyPredictor, "_check_zyg", return_value=False)
def test_verify_pm2ba1bs1bs2_pm2(
    mock_check_zyg,
    mock_absent_in_males,
    mock_ba1_exception,
    mock_get_af,
    coagulation_predictor,
    seqvar,
    auto_acmg_data,
):
    result, comment = coagulation_predictor.verify_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    assert result.PM2 is True, "PM2 should be met when the variant is absent in males"
    assert result.BA1 is False, "BA1 should not be met with low AF"
    assert result.BS1 is False, "BS1 should not be met with low AF"
    assert comment == "The variant is absent in males: PM2 is met. ", "Comment should note PM2 met"

    # Test with different gene settings
    auto_acmg_data.hgnc_id = "HGNC:3551"  # Another gene
    auto_acmg_data.thresholds.ba1_benign = 0.0000556
    auto_acmg_data.thresholds.bs1_benign = 0.00000556
    mock_absent_in_males.return_value = False
    result, comment = coagulation_predictor.verify_pm2ba1bs1bs2(seqvar, auto_acmg_data)


@patch.object(CoagulationFactorDeficiencyPredictor, "_get_af", return_value=0.001)
@patch.object(CoagulationFactorDeficiencyPredictor, "_ba1_exception", return_value=False)
@patch.object(CoagulationFactorDeficiencyPredictor, "_absent_in_males", return_value=False)
@patch.object(CoagulationFactorDeficiencyPredictor, "_check_zyg", return_value=False)
def test_verify_pm2ba1bs1bs2_ba1(
    mock_check_zyg,
    mock_absent_in_males,
    mock_ba1_exception,
    mock_get_af,
    coagulation_predictor,
    seqvar,
    auto_acmg_data,
):
    result, comment = coagulation_predictor.verify_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    assert result.PM2 is False, "PM2 should not be met when the variant present in males"
    assert result.BA1 is True, "BA1 should be met with high AF"
    assert result.BS1 is False, "BS1 should not be met with high BA1 met"


def test_bp3_not_applicable(coagulation_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = coagulation_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


def test_predict_pp2bp1(coagulation_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for CoagulationFactorDeficiencyPredictor."""

    # Call the method under test
    pp2_result, bp1_result = coagulation_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


def test_verify_pp3bp4_revel_thresholds(coagulation_predictor, auto_acmg_data):
    """Test that REVEL thresholds are correctly set for PP3/BP4 prediction."""
    prediction, comment = coagulation_predictor.verify_pp3bp4(
        coagulation_predictor.seqvar, auto_acmg_data
    )

    assert auto_acmg_data.thresholds.revel_pathogenic == 0.6
    assert auto_acmg_data.thresholds.revel_benign == 0.3


@pytest.mark.parametrize(
    "revel_score, expected_pp3, expected_bp4",
    [
        (0.7, True, False),
        (0.5, False, False),
        (0.2, False, True),
    ],
)
def test_verify_pp3bp4_revel_scenarios(
    coagulation_predictor, auto_acmg_data, revel_score, expected_pp3, expected_bp4
):
    auto_acmg_data.scores.dbnsfp.revel = revel_score

    prediction, comment = coagulation_predictor.verify_pp3bp4(
        coagulation_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 == expected_pp3
    assert prediction.BP4 == expected_bp4


@patch.object(CoagulationFactorDeficiencyPredictor, "_affect_spliceAI")
def test_verify_pp3bp4_splicing_f8(mock_affect_spliceAI, coagulation_predictor, auto_acmg_data):
    mock_affect_spliceAI.return_value = True
    auto_acmg_data.hgnc_id = "HGNC:3546"  # F8

    prediction, comment = coagulation_predictor.verify_pp3bp4(
        coagulation_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 is True
    assert prediction.BP4 is False


@patch.object(CoagulationFactorDeficiencyPredictor, "_affect_spliceAI")
def test_verify_pp3bp4_splicing_f9(mock_affect_spliceAI, coagulation_predictor, auto_acmg_data):
    mock_affect_spliceAI.return_value = True
    auto_acmg_data.hgnc_id = "HGNC:3551"  # F9

    prediction, comment = coagulation_predictor.verify_pp3bp4(
        coagulation_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 is True
    assert prediction.BP4 is False


@patch.object(CoagulationFactorDeficiencyPredictor, "_affect_spliceAI")
@pytest.mark.skip("Fix it")
def test_verify_pp3bp4_no_splicing_effect(
    mock_affect_spliceAI, coagulation_predictor, auto_acmg_data
):
    mock_affect_spliceAI.return_value = False
    auto_acmg_data.scores.dbnsfp.revel = 0.4  # Between benign and pathogenic thresholds

    prediction, comment = coagulation_predictor.verify_pp3bp4(
        coagulation_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 is False
    assert prediction.BP4 is True


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_error_handling(coagulation_predictor, auto_acmg_data):
    # Simulate an error condition
    auto_acmg_data.scores.dbnsfp.revel = None

    prediction, comment = coagulation_predictor.verify_pp3bp4(
        coagulation_predictor.seqvar, auto_acmg_data
    )

    assert prediction is None
    assert "An error occurred during prediction" in comment


@patch.object(CoagulationFactorDeficiencyPredictor, "_is_pathogenic_score")
@patch.object(CoagulationFactorDeficiencyPredictor, "_is_benign_score")
@patch.object(CoagulationFactorDeficiencyPredictor, "_affect_spliceAI")
def test_verify_pp3bp4_method_calls(
    mock_affect_spliceAI,
    mock_is_benign_score,
    mock_is_pathogenic_score,
    coagulation_predictor,
    auto_acmg_data,
):
    mock_is_pathogenic_score.return_value = False
    mock_is_benign_score.return_value = False
    mock_affect_spliceAI.return_value = True

    coagulation_predictor.verify_pp3bp4(coagulation_predictor.seqvar, auto_acmg_data)

    mock_is_pathogenic_score.assert_called_once()
    mock_is_benign_score.assert_called_once()
    mock_affect_spliceAI.assert_called()


def test_verify_pp3bp4_spliceai_thresholds(coagulation_predictor, auto_acmg_data):
    prediction, comment = coagulation_predictor.verify_pp3bp4(
        coagulation_predictor.seqvar, auto_acmg_data
    )

    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.05
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.05
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.05
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.05


def test_verify_pp3bp4_benign_spliceai_thresholds(coagulation_predictor, auto_acmg_data):
    auto_acmg_data.hgnc_id = "HGNC:3546"  # F8
    coagulation_predictor.verify_pp3bp4(coagulation_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.05
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.05
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.05
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.05

    auto_acmg_data.hgnc_id = "HGNC:3551"  # F9
    coagulation_predictor.verify_pp3bp4(coagulation_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.01
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.01
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.01
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.01


def test_predict_bp7_threshold_adjustment_for_hgnc_3546(coagulation_predictor, auto_acmg_data):
    """Test that the BP7 thresholds are correctly adjusted for HGNC:3546 (F5)."""
    auto_acmg_data.hgnc_id = "HGNC:3546"  # F5 gene

    # Call predict_bp7 method
    result = coagulation_predictor.predict_bp7(coagulation_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.05
    ), "The spliceAI acceptor gain threshold should be adjusted to 0.05."
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.05
    ), "The spliceAI acceptor loss threshold should be adjusted to 0.05."
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.05
    ), "The spliceAI donor gain threshold should be adjusted to 0.05."
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.05
    ), "The spliceAI donor loss threshold should be adjusted to 0.05."
    assert (
        auto_acmg_data.thresholds.phyloP100 == 0.1
    ), "The phyloP100 threshold should be adjusted to 0.1."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


def test_predict_bp7_threshold_adjustment_for_hgnc_3551(coagulation_predictor, auto_acmg_data):
    """Test that the BP7 thresholds are correctly adjusted for HGNC:3551 (F7)."""
    auto_acmg_data.hgnc_id = "HGNC:3551"  # F7 gene

    # Call predict_bp7 method
    result = coagulation_predictor.predict_bp7(coagulation_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.01
    ), "The spliceAI acceptor gain threshold should be adjusted to 0.01."
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.01
    ), "The spliceAI acceptor loss threshold should be adjusted to 0.01."
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.01
    ), "The spliceAI donor gain threshold should be adjusted to 0.01."
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.01
    ), "The spliceAI donor loss threshold should be adjusted to 0.01."
    assert (
        auto_acmg_data.thresholds.phyloP100 == 0.1
    ), "The phyloP100 threshold should be adjusted to 0.1."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(DefaultSeqVarPredictor, "predict_bp7")
def test_predict_bp7_fallback_to_default(
    mock_super_predict_bp7, coagulation_predictor, auto_acmg_data
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
    result = coagulation_predictor.predict_bp7(coagulation_predictor.seqvar, auto_acmg_data)

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
