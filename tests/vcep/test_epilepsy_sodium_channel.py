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
from src.vcep import EpilepsySodiumChannelPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="2", pos=100, delete="A", insert="T")


@pytest.fixture
def epilepsy_sodium_channel_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return EpilepsySodiumChannelPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


def test_predict_pm1_not_applicable_scn1b(epilepsy_sodium_channel_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for SCN1B."""
    auto_acmg_data.hgnc_id = "HGNC:10586"  # SCN1B gene
    result = epilepsy_sodium_channel_predictor.predict_pm1(
        epilepsy_sodium_channel_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for SCN1B."
    assert (
        result.summary == "PM1 is not applicable for SCN1B."
    ), "The summary should indicate that PM1 is not applicable for SCN1B."


def test_predict_pm1_in_critical_region(epilepsy_sodium_channel_predictor, auto_acmg_data):
    """Test when variant falls within a critical region for SCN1A."""
    auto_acmg_data.hgnc_id = "HGNC:10585"  # SCN1A gene
    auto_acmg_data.prot_pos = 240  # Within the critical region (226-246)
    result = epilepsy_sodium_channel_predictor.predict_pm1(
        epilepsy_sodium_channel_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for a variant in a critical region."
    assert (
        result.summary
        == "Variant falls within a critical residue region for HGNC:10585 between positions 226-246. PM1 is met."
    ), "The summary should indicate the critical region."


def test_predict_pm1_outside_critical_region(epilepsy_sodium_channel_predictor, auto_acmg_data):
    """Test when variant does not fall within any critical region for SCN1A."""
    auto_acmg_data.hgnc_id = "HGNC:10585"  # SCN1A gene
    auto_acmg_data.prot_pos = 300  # Outside the critical region
    result = epilepsy_sodium_channel_predictor.predict_pm1(
        epilepsy_sodium_channel_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met for a variant outside any critical region."
    assert (
        result.summary == "Variant does not meet the PM1 criteria for HGNC:10585."
    ), "The summary should indicate no critical region was met."


@patch("src.vcep.epilepsy_sodium_channel.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, epilepsy_sodium_channel_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method for genes not in PM1_CLUSTER."""
    auto_acmg_data.hgnc_id = "HGNC:99999"  # Gene not in PM1_CLUSTER
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = epilepsy_sodium_channel_predictor.predict_pm1(
        epilepsy_sodium_channel_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met if the gene is not in the PM1_CLUSTER."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should be from the default fallback."


def test_predict_pm1_edge_case_start_boundary(epilepsy_sodium_channel_predictor, auto_acmg_data):
    """Test when variant falls exactly on the start boundary of a critical region."""
    auto_acmg_data.hgnc_id = "HGNC:10585"  # SCN1A gene
    auto_acmg_data.prot_pos = 226  # Start boundary of the critical region (226-246)
    result = epilepsy_sodium_channel_predictor.predict_pm1(
        epilepsy_sodium_channel_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met when on the start boundary of a critical region."
    assert (
        result.summary
        == "Variant falls within a critical residue region for HGNC:10585 between positions 226-246. PM1 is met."
    ), "The summary should indicate the critical region."


def test_predict_pm1_edge_case_end_boundary(epilepsy_sodium_channel_predictor, auto_acmg_data):
    """Test when variant falls exactly on the end boundary of a critical region."""
    auto_acmg_data.hgnc_id = "HGNC:10585"  # SCN1A gene
    auto_acmg_data.prot_pos = 246  # End boundary of the critical region (226-246)
    result = epilepsy_sodium_channel_predictor.predict_pm1(
        epilepsy_sodium_channel_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met when on the end boundary of a critical region."
    assert (
        result.summary
        == "Variant falls within a critical residue region for HGNC:10585 between positions 226-246. PM1 is met."
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
@pytest.mark.parametrize(
    "hgnc_id,expected_ba1,expected_bs1",
    [
        ("HGNC:10585", 0.0002, 0.000004),
        ("HGNC:10588", 0.0001, 0.000002),
        ("HGNC:10590", 0.0001, 0.000002),
        ("HGNC:10596", 0.0001, 0.000002),
        ("HGNC:10586", 0.003, 0.0001),
    ],
)
def test_predict_pm2ba1bs1bs2_with_varied_thresholds(
    mock_super_method,
    epilepsy_sodium_channel_predictor,
    auto_acmg_data,
    seqvar,
    hgnc_id,
    expected_ba1,
    expected_bs1,
):
    # Setup
    auto_acmg_data.hgnc_id = hgnc_id

    # Method call
    epilepsy_sodium_channel_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Validate thresholds are set correctly
    assert (
        auto_acmg_data.thresholds.pm2_pathogenic == 0.000001
    )  # Ensure PM2 threshold is set to a very rare occurrence
    assert auto_acmg_data.thresholds.ba1_benign == expected_ba1
    assert auto_acmg_data.thresholds.bs1_benign == expected_bs1

    # Check that the superclass method was called with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Reset mock for the next iteration
    mock_super_method.reset_mock()


def test_predict_pp2bp1(epilepsy_sodium_channel_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Epilepsy Sodium Channel predictor."""
    # Call the method under test
    pp2_result, bp1_result = epilepsy_sodium_channel_predictor.predict_pp2bp1(
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


def test_predict_pp3bp4_revel_strategy(epilepsy_sodium_channel_predictor, auto_acmg_data):
    """Test that REVEL is set as the strategy for PP3/BP4 prediction."""
    pp3_result, bp4_result = epilepsy_sodium_channel_predictor.predict_pp3bp4(
        epilepsy_sodium_channel_predictor.seqvar, auto_acmg_data
    )

    assert auto_acmg_data.thresholds.pp3bp4_strategy == "revel"


@patch("src.vcep.epilepsy_sodium_channel.DefaultSeqVarPredictor.predict_pp3bp4")
def test_predict_pp3bp4_calls_superclass(
    mock_super_predict_pp3bp4, epilepsy_sodium_channel_predictor, auto_acmg_data
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

    pp3_result, bp4_result = epilepsy_sodium_channel_predictor.predict_pp3bp4(
        epilepsy_sodium_channel_predictor.seqvar, auto_acmg_data
    )

    mock_super_predict_pp3bp4.assert_called_once_with(
        epilepsy_sodium_channel_predictor.seqvar, auto_acmg_data
    )
    assert pp3_result.prediction == AutoACMGPrediction.Applicable
    assert bp4_result.prediction == AutoACMGPrediction.NotApplicable


@pytest.mark.parametrize(
    "revel_score, expected_pp3, expected_bp4",
    [
        (
            0.9,
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
    epilepsy_sodium_channel_predictor,
    auto_acmg_data,
    revel_score,
    expected_pp3,
    expected_bp4,
):
    """Test different REVEL score scenarios."""
    auto_acmg_data.scores.dbnsfp.revel = revel_score

    with patch(
        "src.vcep.epilepsy_sodium_channel.DefaultSeqVarPredictor.predict_pp3bp4"
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

        pp3_result, bp4_result = epilepsy_sodium_channel_predictor.predict_pp3bp4(
            epilepsy_sodium_channel_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == expected_pp3
        assert bp4_result.prediction == expected_bp4


def test_predict_pp3bp4_strength(epilepsy_sodium_channel_predictor, auto_acmg_data):
    """Test that the strength of PP3 and BP4 is correctly set."""
    with patch(
        "src.vcep.epilepsy_sodium_channel.DefaultSeqVarPredictor.predict_pp3bp4"
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

        pp3_result, bp4_result = epilepsy_sodium_channel_predictor.predict_pp3bp4(
            epilepsy_sodium_channel_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.strength == AutoACMGStrength.PathogenicSupporting
        assert bp4_result.strength == AutoACMGStrength.BenignSupporting


def test_predict_pp3bp4_no_revel_score(epilepsy_sodium_channel_predictor, auto_acmg_data):
    """Test behavior when no REVEL score is available."""
    auto_acmg_data.scores.dbnsfp.revel = None

    with patch(
        "src.vcep.epilepsy_sodium_channel.DefaultSeqVarPredictor.predict_pp3bp4"
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

        pp3_result, bp4_result = epilepsy_sodium_channel_predictor.predict_pp3bp4(
            epilepsy_sodium_channel_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == AutoACMGPrediction.NotApplicable
        assert bp4_result.prediction == AutoACMGPrediction.NotApplicable


def test_predict_pp3bp4_error_handling(epilepsy_sodium_channel_predictor, auto_acmg_data):
    """Test error handling in predict_pp3bp4 method."""
    with patch(
        "src.vcep.epilepsy_sodium_channel.DefaultSeqVarPredictor.predict_pp3bp4",
        side_effect=Exception("Test error"),
    ):
        with pytest.raises(Exception) as exc_info:
            epilepsy_sodium_channel_predictor.predict_pp3bp4(
                epilepsy_sodium_channel_predictor.seqvar, auto_acmg_data
            )

        assert str(exc_info.value) == "Test error"
