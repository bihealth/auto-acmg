from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
)
from src.defs.exceptions import AutoAcmgBaseException
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import ENIGMAPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="17", pos=100, delete="A", insert="T")


@pytest.fixture
def enigma_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return ENIGMAPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


def test_predict_pm1_not_applicable_brca1(enigma_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for BRCA1."""
    auto_acmg_data.hgnc_id = "HGNC:1100"  # BRCA1 gene
    result = enigma_predictor.predict_pm1(enigma_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for BRCA1."
    assert (
        result.summary == "PM1 is not applicable for BRCA1 and BRCA2."
    ), "The summary should indicate that PM1 is not applicable for BRCA1."


def test_predict_pm1_not_applicable_brca2(enigma_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for BRCA2."""
    auto_acmg_data.hgnc_id = "HGNC:1101"  # BRCA2 gene
    result = enigma_predictor.predict_pm1(enigma_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for BRCA2."
    assert (
        result.summary == "PM1 is not applicable for BRCA1 and BRCA2."
    ), "The summary should indicate that PM1 is not applicable for BRCA2."


@patch("src.vcep.enigma.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, enigma_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for genes other than BRCA1 and BRCA2."""
    auto_acmg_data.hgnc_id = "HGNC:99999"  # Gene not BRCA1 or BRCA2
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = enigma_predictor.predict_pm1(enigma_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for genes other than BRCA1 and BRCA2."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should be from the default fallback."


def test_predict_pm1_name(enigma_predictor, auto_acmg_data):
    """Test the name of the criteria returned by the ENIGMA predictor."""
    auto_acmg_data.hgnc_id = "HGNC:1100"  # BRCA1 gene
    result = enigma_predictor.predict_pm1(enigma_predictor.seqvar, auto_acmg_data)

    assert result.name == "PM1", "The name of the criteria should be 'PM1'."


def test_predict_pm1_strength(enigma_predictor, auto_acmg_data):
    """Test the strength level returned by the ENIGMA predictor."""
    auto_acmg_data.hgnc_id = "HGNC:1101"  # BRCA2 gene
    result = enigma_predictor.predict_pm1(enigma_predictor.seqvar, auto_acmg_data)

    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate for ENIGMA."


@patch.object(ENIGMAPredictor, "_get_af", return_value=0.1)
@patch.object(ENIGMAPredictor, "_ba1_exception", return_value=False)
def test_verify_pm2ba1bs1bs2(
    mock_get_af, mock_ba1_exception, enigma_predictor, auto_acmg_data, seqvar
):
    # Setup: Adjusting the thresholds to test under different conditions
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.02
    auto_acmg_data.thresholds.pm2_pathogenic = 0.001

    # Call the method under test
    result, comment = enigma_predictor.verify_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assertions to validate the expected behavior
    assert result.BA1 is True, "Expected PM2 to be True based on the mocked allele frequency"

    # Assert changed thresholds
    assert (
        auto_acmg_data.thresholds.ba1_benign == 0.001
    ), "BA1 threshold should be adjusted to 0.001"
    assert (
        auto_acmg_data.thresholds.bs1_benign == 0.00002
    ), "BS1 threshold should be adjusted to 0.00002"


def test_predict_pm4bp3_not_applicable(enigma_predictor, seqvar, auto_acmg_data):
    """Test that PM4 and BP3 are marked as Not Applicable for the brain malformations VCEP."""
    pm4_result, bp3_result = enigma_predictor.predict_pm4bp3(seqvar, auto_acmg_data)

    # Check PM4 result
    assert isinstance(
        pm4_result, AutoACMGCriteria
    ), "The PM4 result should be of type AutoACMGCriteria."
    assert pm4_result.prediction == AutoACMGPrediction.NotApplicable, "PM4 should be NotApplicable."
    assert (
        pm4_result.strength == AutoACMGStrength.PathogenicModerate
    ), "PM4 strength should be PathogenicModerate."
    assert (
        "PM4 is not applicable" in pm4_result.summary
    ), "The summary should indicate PM4 is not applicable."

    # Check BP3 result
    assert isinstance(
        bp3_result, AutoACMGCriteria
    ), "The BP3 result should be of type AutoACMGCriteria."
    assert bp3_result.prediction == AutoACMGPrediction.NotApplicable, "BP3 should be NotApplicable."
    assert (
        bp3_result.strength == AutoACMGStrength.BenignSupporting
    ), "BP3 strength should be BenignSupporting."
    assert (
        "BP3 is not applicable" in bp3_result.summary
    ), "The summary should indicate BP3 is not applicable."


def test_in_important_domain(enigma_predictor, auto_acmg_data):
    """Test if a variant is correctly identified as being in an important domain."""
    auto_acmg_data.hgnc_id = "HGNC:1100"  # BRCA1 gene
    auto_acmg_data.prot_pos = 1695  # Position in an important domain for BRCA1

    important_domain = enigma_predictor._in_important_domain(auto_acmg_data)

    assert important_domain, "The variant should be in an important domain."

    auto_acmg_data.prot_pos = 2000  # Position outside any important domain

    important_domain = enigma_predictor._in_important_domain(auto_acmg_data)

    assert not important_domain, "The variant should not be in an important domain."


def test_predict_pp2bp1_missense_not_in_domain_not_splice_affecting(
    enigma_predictor, seqvar, auto_acmg_data
):
    """Test PP2 and BP1 where the variant is missense, not in an important domain, and does not affect splicing."""
    auto_acmg_data.consequence = MagicMock(mehari=["missense_variant"], cadd=None)
    auto_acmg_data.prot_pos = 2000  # Outside any important domain

    pp2, bp1 = enigma_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert pp2.prediction == AutoACMGPrediction.NotApplicable, "PP2 should be NotApplicable."
    assert (
        bp1.prediction == AutoACMGPrediction.Met
    ), "BP1 should be Met as the criteria are satisfied."
    assert (
        bp1.strength == AutoACMGStrength.PathogenicStrong
    ), "BP1 strength should be PathogenicStrong."
    assert (
        "not in an important domain" in bp1.summary
    ), "The summary should indicate that the variant is not in an important domain."


def test_predict_pp2bp1_missense_in_important_domain(enigma_predictor, seqvar, auto_acmg_data):
    """Test PP2 and BP1 where the variant is missense and in an important domain."""
    auto_acmg_data.consequence = MagicMock(mehari=["missense_variant"], cadd=None)
    auto_acmg_data.prot_pos = 1700  # Within an important domain

    with patch.object(ENIGMAPredictor, "_in_important_domain", return_value=True):
        pp2, bp1 = enigma_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert pp2.prediction == AutoACMGPrediction.NotApplicable, "PP2 should be NotApplicable."
    assert (
        bp1.prediction == AutoACMGPrediction.NotMet
    ), "BP1 should not be Met as the variant is in an important domain."


def test_predict_pp2bp1_synonymous_splice_affecting(enigma_predictor, seqvar, auto_acmg_data):
    """Test PP2 and BP1 where the variant is synonymous and predicted to affect splicing."""
    auto_acmg_data.consequence = MagicMock(mehari=["synonymous_variant"], cadd=None)
    auto_acmg_data.prot_pos = 2000  # Outside any important domain

    with patch.object(ENIGMAPredictor, "_spliceai_impact", return_value=True):
        pp2, bp1 = enigma_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert pp2.prediction == AutoACMGPrediction.NotApplicable, "PP2 should be NotApplicable."
    assert (
        bp1.prediction == AutoACMGPrediction.NotMet
    ), "BP1 should not be Met as the variant is predicted to affect splicing."


def test_predict_pp2bp1_inframe_deletion_not_in_important_domain_not_splice_affecting(
    enigma_predictor, seqvar, auto_acmg_data
):
    """Test PP2 and BP1 for an inframe deletion that is not in an important domain and does not affect splicing."""
    auto_acmg_data.consequence = MagicMock(mehari=["inframe_deletion"], cadd=None)
    auto_acmg_data.prot_pos = 2500  # Outside any important domain

    pp2, bp1 = enigma_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert pp2.prediction == AutoACMGPrediction.NotApplicable, "PP2 should be NotApplicable."
    assert (
        bp1.prediction == AutoACMGPrediction.Met
    ), "BP1 should be Met as the criteria are satisfied."
    assert (
        "not in an important domain" in bp1.summary
    ), "The summary should reflect the non-impact on important domains."


def test_predict_pp2bp1_non_relevant_variant_type(enigma_predictor, seqvar, auto_acmg_data):
    """Test PP2 and BP1 for a variant type that is not missense, synonymous, or inframe (e.g., frameshift)."""
    auto_acmg_data.consequence = MagicMock(mehari=["frameshift_variant"], cadd=None)
    auto_acmg_data.prot_pos = 2000  # Irrelevant due to variant type

    pp2, bp1 = enigma_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert pp2.prediction == AutoACMGPrediction.NotApplicable, "PP2 should be NotApplicable."
    assert (
        bp1.prediction == AutoACMGPrediction.NotMet
    ), "BP1 should not be Met due to variant type being not relevant."
    assert (
        "not synonymous, missense or inframe indel" in bp1.summary
    ), "The summary should indicate the variant type."


@patch.object(ENIGMAPredictor, "_in_important_domain", return_value=True)
@patch.object(ENIGMAPredictor, "_is_synonymous", return_value=True)
@patch.object(ENIGMAPredictor, "_is_intronic", return_value=False)
@patch.object(ENIGMAPredictor, "_affect_canonical_ss", return_value=False)
@patch.object(ENIGMAPredictor, "_is_conserved", return_value=False)
def test_verify_bp7_met(
    mock_is_conserved,
    mock_affect_canonical_ss,
    mock_is_intronic,
    mock_is_synonymous,
    mock_in_important_domain,
    enigma_predictor,
    auto_acmg_data,
):
    """Test that BP7 is met when conditions are satisfied."""
    auto_acmg_data.hgnc_id = "HGNC:1100"  # BRCA1 gene

    prediction_bp7, comment_bp7 = enigma_predictor.verify_bp7(
        enigma_predictor.seqvar, auto_acmg_data
    )

    assert prediction_bp7.BP7, "BP7 should be met when the conditions are satisfied."
    assert "BP7 is met" in comment_bp7, "The comment should indicate that BP7 is met."


@patch.object(ENIGMAPredictor, "_in_important_domain", return_value=False)
@patch.object(ENIGMAPredictor, "_is_synonymous", return_value=True)
@patch.object(ENIGMAPredictor, "_is_intronic", return_value=False)
@patch.object(ENIGMAPredictor, "_affect_canonical_ss", return_value=False)
@patch.object(ENIGMAPredictor, "_is_conserved", return_value=False)
def test_verify_bp7_not_met(
    mock_is_conserved,
    mock_affect_canonical_ss,
    mock_is_intronic,
    mock_is_synonymous,
    mock_in_important_domain,
    enigma_predictor,
    auto_acmg_data,
):
    """Test that BP7 is not met when the conditions are not satisfied."""
    auto_acmg_data.hgnc_id = "HGNC:1100"  # BRCA1 gene

    prediction_bp7, comment_bp7 = enigma_predictor.verify_bp7(
        enigma_predictor.seqvar, auto_acmg_data
    )

    assert not prediction_bp7.BP7, "BP7 should not be met when the conditions are not satisfied."
    assert "BP7 is not met" in comment_bp7, "The comment should indicate that BP7 is not met."


@patch.object(ENIGMAPredictor, "_in_important_domain", side_effect=AutoAcmgBaseException("Error"))
@patch.object(ENIGMAPredictor, "_is_synonymous", return_value=True)
def test_verify_bp7_exception_handling(
    mock_in_important_domain, mock_is_synonymous, enigma_predictor, auto_acmg_data
):
    """Test that BP7 prediction handles exceptions properly."""
    auto_acmg_data.hgnc_id = "HGNC:1100"  # BRCA1 gene

    prediction_bp7, comment_bp7 = enigma_predictor.verify_bp7(
        enigma_predictor.seqvar, auto_acmg_data
    )

    assert prediction_bp7 is None, "BP7 prediction should be None when an exception occurs."
    assert (
        "Failed to predict BP7 criterion" in comment_bp7
    ), "The comment should indicate a failure."


def test_verify_bp7_threshold_adjustment(enigma_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted for BRCA1 and BRCA2."""
    auto_acmg_data.thresholds.bp7_donor = 1  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    auto_acmg_data.hgnc_id = "HGNC:1100"  # BRCA1 gene

    enigma_predictor.verify_bp7(enigma_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.bp7_donor == 7
    ), "The BP7 donor threshold should be adjusted to 7."
    assert (
        auto_acmg_data.thresholds.bp7_acceptor == 21
    ), "The BP7 acceptor threshold should be adjusted to 21."
