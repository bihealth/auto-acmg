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
from src.vcep import HearingLossPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def hearing_loss_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return HearingLossPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_overrides(mock_super_verify, hearing_loss_predictor, seqvar, auto_acmg_data):
    """Test that the overridden verify_ps1pm5 method in HearingLossPredictor works correctly."""
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
    prediction, comment = hearing_loss_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that the splicing effect leads to overriding PS1 and PM5 as not applicable
    assert not prediction.PS1, "PS1 should be marked as not applicable due to splicing effect."
    assert not prediction.PM5, "PM5 should be marked as not applicable due to splicing effect."
    assert "Variant affects splicing" in comment, "Comment should note the splicing effect."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_no_splicing_effect(
    mock_super_verify, hearing_loss_predictor, seqvar, auto_acmg_data
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
    prediction, comment = hearing_loss_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that PS1 and PM5 remain applicable
    assert prediction.PS1, "PS1 should remain applicable when there's no splicing effect."
    assert prediction.PM5, "PM5 should remain applicable when there's no splicing effect."
    assert "Initial evaluation" in comment, "Comment should reflect the initial evaluation."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_non_missense(
    mock_super_verify, hearing_loss_predictor, seqvar, auto_acmg_data
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
    prediction, comment = hearing_loss_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that PS1 and PM5 remain as per superclass evaluation
    assert not prediction.PS1, "PS1 should remain not applicable for non-missense variants."
    assert not prediction.PM5, "PM5 should remain not applicable for non-missense variants."
    assert (
        "Not applicable for non-missense" in comment
    ), "Comment should reflect superclass evaluation."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_exception_handling(
    mock_super_verify, hearing_loss_predictor, seqvar, auto_acmg_data
):
    """Test verify_ps1pm5 method exception handling."""
    mock_super_verify.side_effect = Exception("Test exception")

    with pytest.raises(Exception) as exc_info:
        hearing_loss_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    assert "Test exception" in str(exc_info.value), "Should raise the original exception"


def test_predict_pm1_kcnq4_strong(hearing_loss_predictor, auto_acmg_data):
    """Test when PM1 is met at the Strong level for a variant in KCNQ4."""
    auto_acmg_data.hgnc_id = "HGNC:6298"  # KCNQ4 gene
    auto_acmg_data.prot_pos = 280  # Within the critical pore-forming intramembrane region
    result = hearing_loss_predictor.predict_pm1(hearing_loss_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met at the Strong level."
    assert (
        result.strength == AutoACMGStrength.PathogenicStrong
    ), "The strength should be PathogenicStrong."
    assert (
        "critical pore-forming intramembrane region" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_not_applicable(hearing_loss_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for other hearing loss genes."""
    auto_acmg_data.hgnc_id = "HGNC:4284"  # GJB2 gene
    result = hearing_loss_predictor.predict_pm1(hearing_loss_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for GJB2."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not applicable" in result.summary
    ), "The summary should indicate non-applicability."


def test_predict_pm1_not_met(hearing_loss_predictor, auto_acmg_data):
    """Test when PM1 is not met for KCNQ4 but outside the critical region."""
    auto_acmg_data.hgnc_id = "HGNC:6298"  # KCNQ4 gene
    auto_acmg_data.prot_pos = 300  # Outside the critical pore-forming intramembrane region
    result = hearing_loss_predictor.predict_pm1(hearing_loss_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.NotApplicable, "PM1 should not be met for KCNQ4."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant does not meet the PM1 criteria" in result.summary
    ), "The summary should indicate no criteria were met."


@patch("src.vcep.hearing_loss.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, hearing_loss_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method if logic changes."""
    auto_acmg_data.hgnc_id = "HGNC:111111"  # Gene not in the specific logic
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = hearing_loss_predictor.predict_pm1(hearing_loss_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


def test_predict_pm1_edge_case_start_boundary(hearing_loss_predictor, auto_acmg_data):
    """Test when variant falls exactly on the start boundary of the critical region."""
    auto_acmg_data.hgnc_id = "HGNC:6298"  # KCNQ4 gene
    auto_acmg_data.prot_pos = (
        271  # Start boundary of the critical pore-forming intramembrane region
    )
    result = hearing_loss_predictor.predict_pm1(hearing_loss_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met on the start boundary."
    assert (
        "critical pore-forming intramembrane region" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_edge_case_end_boundary(hearing_loss_predictor, auto_acmg_data):
    """Test when variant falls exactly on the end boundary of the critical region."""
    auto_acmg_data.hgnc_id = "HGNC:6298"  # KCNQ4 gene
    auto_acmg_data.prot_pos = 292  # End boundary of the critical pore-forming intramembrane region
    result = hearing_loss_predictor.predict_pm1(hearing_loss_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met on the end boundary."
    assert (
        "critical pore-forming intramembrane region" in result.summary
    ), "The summary should indicate the critical region."


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
def test_predict_pm2ba1bs1bs2(mock_super_method, hearing_loss_predictor, auto_acmg_data, seqvar):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = hearing_loss_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.00007
    assert auto_acmg_data.thresholds.ba1_benign == 0.001
    assert auto_acmg_data.thresholds.bs1_benign == 0.0007

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_predict_pp2bp1(hearing_loss_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Hearing Loss predictor."""

    # Call the method under test
    pp2_result, bp1_result = hearing_loss_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


def test_predict_pp3bp4_revel_strategy(hearing_loss_predictor, auto_acmg_data):
    """Test that REVEL is set as the strategy for PP3/BP4 prediction."""
    hearing_loss_predictor.predict_pp3bp4(hearing_loss_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.pp3bp4_strategy == "revel"


def test_predict_pp3bp4_revel_thresholds(hearing_loss_predictor, auto_acmg_data):
    """Test that REVEL thresholds are correctly set for PP3/BP4 prediction."""
    hearing_loss_predictor.predict_pp3bp4(hearing_loss_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.revel_pathogenic == 0.7
    assert auto_acmg_data.thresholds.revel_benign == 0.15


@patch("src.vcep.hearing_loss.DefaultSeqVarPredictor.predict_pp3bp4")
def test_predict_pp3bp4_calls_superclass(
    mock_super_predict_pp3bp4, hearing_loss_predictor, auto_acmg_data
):
    """Test that the superclass method is called with the correct parameters."""
    mock_super_predict_pp3bp4.return_value = (
        AutoACMGCriteria(
            name="PP3",
            prediction=AutoACMGPrediction.Applicable,
            strength=AutoACMGStrength.PathogenicSupporting,
        ),
        AutoACMGCriteria(
            name="BP4",
            prediction=AutoACMGPrediction.NotApplicable,
            strength=AutoACMGStrength.BenignSupporting,
        ),
    )

    pp3_result, bp4_result = hearing_loss_predictor.predict_pp3bp4(
        hearing_loss_predictor.seqvar, auto_acmg_data
    )

    mock_super_predict_pp3bp4.assert_called_once_with(hearing_loss_predictor.seqvar, auto_acmg_data)
    assert pp3_result.prediction == AutoACMGPrediction.Applicable
    assert bp4_result.prediction == AutoACMGPrediction.NotApplicable


@pytest.mark.parametrize(
    "revel_score, expected_pp3, expected_bp4",
    [
        (
            0.8,
            AutoACMGPrediction.Applicable,
            AutoACMGPrediction.NotApplicable,
        ),  # High REVEL score
        (
            0.5,
            AutoACMGPrediction.NotApplicable,
            AutoACMGPrediction.NotApplicable,
        ),  # Intermediate REVEL score
        (
            0.1,
            AutoACMGPrediction.NotApplicable,
            AutoACMGPrediction.Applicable,
        ),  # Low REVEL score
    ],
)
def test_predict_pp3bp4_revel_scenarios(
    hearing_loss_predictor, auto_acmg_data, revel_score, expected_pp3, expected_bp4
):
    """Test different REVEL score scenarios."""
    auto_acmg_data.scores.dbnsfp.revel = revel_score

    with patch(
        "src.vcep.hearing_loss.DefaultSeqVarPredictor.predict_pp3bp4"
    ) as mock_super_predict_pp3bp4:
        mock_super_predict_pp3bp4.return_value = (
            AutoACMGCriteria(
                name="PP3",
                prediction=expected_pp3,
                strength=AutoACMGStrength.PathogenicSupporting,
            ),
            AutoACMGCriteria(
                name="BP4",
                prediction=expected_bp4,
                strength=AutoACMGStrength.BenignSupporting,
            ),
        )

        pp3_result, bp4_result = hearing_loss_predictor.predict_pp3bp4(
            hearing_loss_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == expected_pp3
        assert bp4_result.prediction == expected_bp4


def test_predict_pp3bp4_strength(hearing_loss_predictor, auto_acmg_data):
    """Test that the strength of PP3 and BP4 is correctly set."""
    with patch(
        "src.vcep.hearing_loss.DefaultSeqVarPredictor.predict_pp3bp4"
    ) as mock_super_predict_pp3bp4:
        mock_super_predict_pp3bp4.return_value = (
            AutoACMGCriteria(
                name="PP3",
                prediction=AutoACMGPrediction.Applicable,
                strength=AutoACMGStrength.PathogenicSupporting,
            ),
            AutoACMGCriteria(
                name="BP4",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.BenignSupporting,
            ),
        )

        pp3_result, bp4_result = hearing_loss_predictor.predict_pp3bp4(
            hearing_loss_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.strength == AutoACMGStrength.PathogenicSupporting
        assert bp4_result.strength == AutoACMGStrength.BenignSupporting


def test_predict_pp3bp4_no_revel_score(hearing_loss_predictor, auto_acmg_data):
    """Test behavior when no REVEL score is available."""
    auto_acmg_data.scores.dbnsfp.revel = None

    with patch(
        "src.vcep.hearing_loss.DefaultSeqVarPredictor.predict_pp3bp4"
    ) as mock_super_predict_pp3bp4:
        mock_super_predict_pp3bp4.return_value = (
            AutoACMGCriteria(
                name="PP3",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
            ),
            AutoACMGCriteria(
                name="BP4",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.BenignSupporting,
            ),
        )

        pp3_result, bp4_result = hearing_loss_predictor.predict_pp3bp4(
            hearing_loss_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == AutoACMGPrediction.NotApplicable
        assert bp4_result.prediction == AutoACMGPrediction.NotApplicable


def test_predict_pp3bp4_error_handling(hearing_loss_predictor, auto_acmg_data):
    """Test error handling in predict_pp3bp4 method."""
    with patch(
        "src.vcep.hearing_loss.DefaultSeqVarPredictor.predict_pp3bp4",
        side_effect=Exception("Test error"),
    ):
        with pytest.raises(Exception) as exc_info:
            hearing_loss_predictor.predict_pp3bp4(hearing_loss_predictor.seqvar, auto_acmg_data)

        assert str(exc_info.value) == "Test error"


def test_predict_pp3bp4_summary(hearing_loss_predictor, auto_acmg_data):
    """Test that the summary of PP3 and BP4 is correctly set."""
    with patch(
        "src.vcep.hearing_loss.DefaultSeqVarPredictor.predict_pp3bp4"
    ) as mock_super_predict_pp3bp4:
        mock_super_predict_pp3bp4.return_value = (
            AutoACMGCriteria(
                name="PP3",
                prediction=AutoACMGPrediction.Applicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary="REVEL score indicates pathogenicity",
            ),
            AutoACMGCriteria(
                name="BP4",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.BenignSupporting,
                summary="REVEL score does not indicate benign",
            ),
        )

        pp3_result, bp4_result = hearing_loss_predictor.predict_pp3bp4(
            hearing_loss_predictor.seqvar, auto_acmg_data
        )

        assert "REVEL score indicates pathogenicity" in pp3_result.summary
        assert "REVEL score does not indicate benign" in bp4_result.summary


def test_predict_pp3bp4_edge_cases(hearing_loss_predictor, auto_acmg_data):
    """Test edge cases for REVEL scores."""
    edge_cases = [0.7, 0.15]  # Exactly at the thresholds
    for score in edge_cases:
        auto_acmg_data.scores.dbnsfp.revel = score
        with patch(
            "src.vcep.hearing_loss.DefaultSeqVarPredictor.predict_pp3bp4"
        ) as mock_super_predict_pp3bp4:
            mock_super_predict_pp3bp4.return_value = (
                AutoACMGCriteria(
                    name="PP3",
                    prediction=(
                        AutoACMGPrediction.Applicable
                        if score >= 0.7
                        else AutoACMGPrediction.NotApplicable
                    ),
                    strength=AutoACMGStrength.PathogenicSupporting,
                ),
                AutoACMGCriteria(
                    name="BP4",
                    prediction=(
                        AutoACMGPrediction.Applicable
                        if score <= 0.15
                        else AutoACMGPrediction.NotApplicable
                    ),
                    strength=AutoACMGStrength.BenignSupporting,
                ),
            )

            pp3_result, bp4_result = hearing_loss_predictor.predict_pp3bp4(
                hearing_loss_predictor.seqvar, auto_acmg_data
            )

            if score == 0.7:
                assert pp3_result.prediction == AutoACMGPrediction.Applicable
                assert bp4_result.prediction == AutoACMGPrediction.NotApplicable
            elif score == 0.15:
                assert pp3_result.prediction == AutoACMGPrediction.NotApplicable
                assert bp4_result.prediction == AutoACMGPrediction.Applicable
