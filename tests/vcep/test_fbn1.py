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
from src.vcep import FBN1Predictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="15", pos=100, delete="A", insert="T")


@pytest.fixture
def fbn1_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return FBN1Predictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_overrides(mock_super_verify, fbn1_predictor, seqvar, auto_acmg_data):
    """Test that the overridden verify_ps1pm5 method in FBN1Predictor works correctly."""
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
        cadd=MagicMock(spliceAI_acceptor_gain=0.6, spliceAI_donor_gain=0.6)
    )

    # Run the method under test
    prediction, comment = fbn1_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that the splicing effect leads to overriding PS1 and PM5 as not applicable
    assert not prediction.PS1, "PS1 should be marked as not applicable due to splicing effect."
    assert not prediction.PM5, "PM5 should be marked as not applicable due to splicing effect."
    assert "Variant affects splicing" in comment, "Comment should note the splicing effect."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_no_splicing_effect(
    mock_super_verify, fbn1_predictor, seqvar, auto_acmg_data
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
            spliceAI_acceptor_gain=0.4,
            spliceAI_acceptor_loss=0.4,
            spliceAI_donor_gain=0.4,
            spliceAI_donor_loss=0.4,
        )
    )

    # Run the method under test
    prediction, comment = fbn1_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that PS1 and PM5 remain applicable
    assert prediction.PS1, "PS1 should remain applicable when there's no splicing effect."
    assert prediction.PM5, "PM5 should remain applicable when there's no splicing effect."
    assert "Initial evaluation" in comment, "Comment should reflect the initial evaluation."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_non_missense(mock_super_verify, fbn1_predictor, seqvar, auto_acmg_data):
    """Test that PS1 and PM5 remain as per superclass for non-missense variants."""
    # Set up the mock to return PS1 and PM5 as not applicable initially
    mock_super_verify.return_value = (
        PS1PM5(PS1=False, PM5=False),
        "Not applicable for non-missense",
    )

    # Setup the data with a non-missense variant
    auto_acmg_data.consequence = MagicMock(mehari=["synonymous_variant"])

    # Run the method under test
    prediction, comment = fbn1_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that PS1 and PM5 remain as per superclass evaluation
    assert not prediction.PS1, "PS1 should remain not applicable for non-missense variants."
    assert not prediction.PM5, "PM5 should remain not applicable for non-missense variants."
    assert (
        "Not applicable for non-missense" in comment
    ), "Comment should reflect superclass evaluation."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)


def test_predict_pm1_strong_residue(fbn1_predictor, auto_acmg_data):
    """Test when variant affects a strong critical residue in FBN1."""
    auto_acmg_data.hgnc_id = "HGNC:3603"  # FBN1 gene
    auto_acmg_data.prot_pos = 250  # Strong residue in FBN1
    result = fbn1_predictor.predict_pm1(fbn1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for a variant affecting a strong critical residue."
    assert (
        result.strength == AutoACMGStrength.PathogenicStrong
    ), "The strength should be PathogenicStrong for strong residues."
    assert (
        "affects a critical cysteine residue" in result.summary
    ), "The summary should indicate the strong critical residue."


def test_predict_pm1_moderate_residue(fbn1_predictor, auto_acmg_data):
    """Test when variant affects a moderate critical residue in FBN1."""
    auto_acmg_data.hgnc_id = "HGNC:3603"  # FBN1 gene
    auto_acmg_data.prot_pos = 100  # Moderate residue in FBN1
    result = fbn1_predictor.predict_pm1(fbn1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for a variant affecting a moderate critical residue."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate for moderate residues."
    assert (
        "Variant affects a residue in FBN1" in result.summary
    ), "The summary should indicate the moderate critical residue."


def test_predict_pm1_outside_critical_residues(fbn1_predictor, auto_acmg_data):
    """Test when variant does not fall within any critical residues in FBN1."""
    auto_acmg_data.hgnc_id = "HGNC:3603"  # FBN1 gene
    auto_acmg_data.prot_pos = 5000  # Position outside critical residues
    result = fbn1_predictor.predict_pm1(fbn1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for a variant outside critical residues."
    assert (
        "Variant does not meet the PM1 criteria for FBN1" in result.summary
    ), "The summary should indicate no criteria were met."


@patch("src.vcep.fbn1.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, fbn1_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for genes not in PM1_CLUSTER."""
    auto_acmg_data.hgnc_id = "HGNC:99999"  # Gene not in PM1_CLUSTER
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = fbn1_predictor.predict_pm1(fbn1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met if the gene is not in the PM1_CLUSTER."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should be from the default fallback."


def test_predict_pm1_edge_case_strong_boundary(fbn1_predictor, auto_acmg_data):
    """Test when variant falls exactly on the boundary of a strong critical residue."""
    auto_acmg_data.hgnc_id = "HGNC:3603"  # FBN1 gene
    auto_acmg_data.prot_pos = 2686  # Boundary of a strong critical residue in FBN1
    result = fbn1_predictor.predict_pm1(fbn1_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the boundary of a strong critical residue."
    assert (
        result.strength == AutoACMGStrength.PathogenicStrong
    ), "The strength should be PathogenicStrong for strong residues."
    assert (
        "affects a critical cysteine residue" in result.summary
    ), "The summary should indicate the strong critical residue."


def test_predict_pm1_edge_case_moderate_boundary(fbn1_predictor, auto_acmg_data):
    """Test when variant falls exactly on the boundary of a moderate critical residue."""
    auto_acmg_data.hgnc_id = "HGNC:3603"  # FBN1 gene
    auto_acmg_data.prot_pos = 2390  # Boundary of a moderate critical residue in FBN1
    result = fbn1_predictor.predict_pm1(fbn1_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the boundary of a moderate critical residue."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate for moderate residues."
    assert (
        "Variant affects a residue in FBN1" in result.summary
    ), "The summary should indicate the moderate critical residue."


def test_bs2_not_applicable(fbn1_predictor, auto_acmg_data):
    """Test BS2 is not applicable for FBN1 as overridden."""
    result = fbn1_predictor._bs2_not_applicable(auto_acmg_data)
    assert result is True, "BS2 should always be not applicable for FBN1."


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
def test_predict_pm2ba1bs1bs2(mock_super_method, fbn1_predictor, auto_acmg_data, seqvar):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = fbn1_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.000005
    assert auto_acmg_data.thresholds.ba1_benign == 0.001
    assert auto_acmg_data.thresholds.bs1_benign == 0.00005

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_bp3_not_applicable(fbn1_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = fbn1_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


@patch.object(FBN1Predictor, "verify_pp2bp1")
def test_predict_pp2bp1_applicable(mock_verify, fbn1_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 where variant prediction is applicable."""
    mock_verify.return_value = (
        MagicMock(PP2=True, PP2_strength=AutoACMGStrength.PathogenicModerate),
        "PP2 applicable.",
    )
    pp2, bp1 = fbn1_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert pp2.prediction == AutoACMGPrediction.Met, "PP2 should be Met."
    assert (
        pp2.strength == AutoACMGStrength.PathogenicModerate
    ), "PP2 strength should be PathogenicModerate."
    assert "PP2 applicable." in pp2.summary, "The summary should confirm PP2 applicability."
    assert bp1.prediction == AutoACMGPrediction.NotApplicable, "BP1 should be NotApplicable."


@patch.object(FBN1Predictor, "verify_pp2bp1")
def test_predict_pp2bp1_not_applicable(mock_verify, fbn1_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 where variant prediction is not applicable."""
    mock_verify.return_value = (
        MagicMock(PP2=False, PP2_strength=AutoACMGStrength.PathogenicSupporting),
        "PP2 not applicable.",
    )
    pp2, bp1 = fbn1_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert pp2.prediction == AutoACMGPrediction.NotMet, "PP2 should be NotMet."
    assert "PP2 not applicable." in pp2.summary, "The summary should confirm PP2 non-applicability."
    assert bp1.prediction == AutoACMGPrediction.NotApplicable, "BP1 should be NotApplicable."


@patch.object(FBN1Predictor, "verify_pp2bp1")
def test_predict_pp2bp1_failed(mock_verify, fbn1_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 when the prediction process fails."""
    mock_verify.return_value = (None, "Prediction process failed.")
    pp2, bp1 = fbn1_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert pp2.prediction == AutoACMGPrediction.Failed, "PP2 prediction should fail."
    assert "Prediction process failed." in pp2.summary, "The summary should report the failure."
    assert bp1.prediction == AutoACMGPrediction.NotApplicable, "BP1 should be NotApplicable."


@patch.object(FBN1Predictor, "verify_pp2bp1")
def test_predict_pp2bp1_exception_handling(mock_verify, fbn1_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 handling of unexpected exceptions."""
    mock_verify.side_effect = Exception("Internal error")
    with pytest.raises(Exception) as exc_info:
        pp2, bp1 = fbn1_predictor.predict_pp2bp1(seqvar, auto_acmg_data)
    assert "Internal error" in str(exc_info.value), "Should raise an internal error exception."
