from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
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
    return AutoACMGData()


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
        result.prediction == AutoACMGPrediction.Met
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
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for a variant outside any critical region."
    assert (
        result.summary == "Variant does not meet the PM1 criteria for HGNC:10585."
    ), "The summary should indicate no critical region was met."


@patch("src.vcep.epilepsy_sodium_channel.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, epilepsy_sodium_channel_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method for genes not in PM1_CLUSTER."""
    auto_acmg_data.hgnc_id = "HGNC:99999"  # Gene not in PM1_CLUSTER
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = epilepsy_sodium_channel_predictor.predict_pm1(
        epilepsy_sodium_channel_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
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
        result.prediction == AutoACMGPrediction.Met
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
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the end boundary of a critical region."
    assert (
        result.summary
        == "Variant falls within a critical residue region for HGNC:10585 between positions 226-246. PM1 is met."
    ), "The summary should indicate the critical region."


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
