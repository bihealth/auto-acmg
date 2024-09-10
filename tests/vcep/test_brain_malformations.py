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
from src.vcep import BrainMalformationsPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


@pytest.fixture
def brain_malformations_predictor(seqvar, auto_acmg_data):
    result = MagicMock()  # Mocking the AutoACMGResult object, which will be passed along
    return BrainMalformationsPredictor(seqvar=seqvar, result=result, config=MagicMock())


def test_predict_pvs1_not_applicable(brain_malformations_predictor, seqvar, auto_acmg_data):
    result = brain_malformations_predictor.predict_pvs1(seqvar, auto_acmg_data)

    assert result.name == "PVS1", "The criterion name should be 'PVS1'"
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "The prediction should be NotApplicable"
    assert (
        result.strength == AutoACMGStrength.PathogenicVeryStrong
    ), "The strength should be PathogenicVeryStrong"
    assert (
        result.summary == "PVS1 is not applicable for the gene."
    ), "The summary should indicate that PVS1 is not applicable"


def test_predict_pm1_in_critical_domain(brain_malformations_predictor, auto_acmg_data):
    """Test when variant falls within a critical domain for the specified gene."""
    auto_acmg_data.prot_pos = 100  # Set protein position within a critical domain
    auto_acmg_data.hgnc_id = "HGNC:393"  # AKT3 gene
    result = brain_malformations_predictor.predict_pm1(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for critical domain variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."


def test_predict_pm1_outside_critical_domain(brain_malformations_predictor, auto_acmg_data):
    """Test when variant does not fall within any critical domain for the specified gene."""
    auto_acmg_data.prot_pos = 500  # Set protein position outside all critical domains
    auto_acmg_data.hgnc_id = "HGNC:393"  # AKT3 gene
    result = brain_malformations_predictor.predict_pm1(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met for non-critical domain variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "does not fall within any critical domain" in result.summary
    ), "The summary should indicate no critical domain."


def test_predict_pm1_edge_case_start_boundary(brain_malformations_predictor, auto_acmg_data):
    """Test when variant falls exactly on the start boundary of a critical domain."""
    auto_acmg_data.prot_pos = 1382  # Start boundary of the MTOR kinase domain (1382-1982)
    auto_acmg_data.hgnc_id = "HGNC:3942"  # MTOR gene
    result = brain_malformations_predictor.predict_pm1(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met when on the start boundary of a critical domain."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."


def test_predict_pm1_edge_case_end_boundary(brain_malformations_predictor, auto_acmg_data):
    """Test when variant falls exactly on the end boundary of a critical domain."""
    auto_acmg_data.prot_pos = 1982  # End boundary of the MTOR kinase domain (1382-1982)
    auto_acmg_data.hgnc_id = "HGNC:3942"  # MTOR gene
    result = brain_malformations_predictor.predict_pm1(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met when on the end boundary of a critical domain."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."


@patch.object(DefaultSeqVarPredictor, "predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_super_predict_pm1, brain_malformations_predictor, auto_acmg_data
):
    """Test when the transcript ID is not in the PM1_CLUSTER mapping, it should fallback to default PM1 prediction."""
    auto_acmg_data.hgnc_id = "HGNC:1111111111111"  # Not in the PM1_CLUSTER mapping

    # Set the mock return value for the superclass's predict_pm1 method
    mock_super_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )

    result = brain_malformations_predictor.predict_pm1(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met if no specific cluster mapping is found."
    assert mock_super_predict_pm1.called, "super().predict_pm1 should have been called."


@patch.object(BrainMalformationsPredictor, "_get_af", return_value=0.1)
@patch.object(BrainMalformationsPredictor, "_ba1_exception", return_value=False)
def test_verify_pm2ba1bs1bs2(
    mock_get_af,
    mock_ba1_exception,
    brain_malformations_predictor,
    auto_acmg_data,
    seqvar,
):
    # Setup: Adjusting the thresholds to test under different conditions
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.02
    auto_acmg_data.thresholds.pm2_pathogenic = 0.001

    # Call the method under test
    result, comment = brain_malformations_predictor.verify_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assertions to validate the expected behavior
    assert result.BA1 is True, "Expected PM2 to be True based on the mocked allele frequency"

    # Assert changed thresholds
    assert (
        auto_acmg_data.thresholds.ba1_benign == 0.000926
    ), "BA1 threshold should be adjusted to 0.000926"
    assert (
        auto_acmg_data.thresholds.bs1_benign == 0.000185
    ), "BS1 threshold should be adjusted to 0.000185"


def test_predict_pm4bp3_not_applicable(brain_malformations_predictor, seqvar, auto_acmg_data):
    """Test that PM4 and BP3 are marked as Not Applicable for the brain malformations VCEP."""
    pm4_result, bp3_result = brain_malformations_predictor.predict_pm4bp3(seqvar, auto_acmg_data)

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


def test_predict_pp2bp1_missense(brain_malformations_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 when the variant is a missense change in a relevant gene."""
    auto_acmg_data.hgnc_id = "HGNC:393"  # Relevant gene
    auto_acmg_data.consequence.cadd = "missense_variant"  # Missense change

    # Call the method under test
    pp2_result, bp1_result = brain_malformations_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    # Check PP2 result
    assert isinstance(
        pp2_result, AutoACMGCriteria
    ), "The PP2 result should be of type AutoACMGCriteria."
    assert (
        pp2_result.prediction == AutoACMGPrediction.Applicable
    ), "PP2 should be Met for a missense variant in a relevant gene."
    assert "PP2 is met for HGNC:393 as the variant is a missense change." in pp2_result.summary

    # Check BP1 result
    assert isinstance(
        bp1_result, AutoACMGCriteria
    ), "The BP1 result should be of type AutoACMGCriteria."
    assert (
        bp1_result.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should be NotApplicable for this gene."
    assert "BP1 is not applicable for the gene." in bp1_result.summary


def test_predict_pp2bp1_not_missense(brain_malformations_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 when the variant is not a missense change in a relevant gene."""
    auto_acmg_data.hgnc_id = "HGNC:393"  # Relevant gene
    auto_acmg_data.consequence.cadd = "synonymous_variant"  # Not a missense change

    # Call the method under test
    pp2_result, bp1_result = brain_malformations_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    # Check PP2 result
    assert isinstance(
        pp2_result, AutoACMGCriteria
    ), "The PP2 result should be of type AutoACMGCriteria."
    assert (
        pp2_result.prediction == AutoACMGPrediction.NotApplicable
    ), "PP2 should not be met for a non-missense variant in a relevant gene."
    assert (
        "PP2 is not met for HGNC:393 as the variant is not a missense change." in pp2_result.summary
    )

    # Check BP1 result
    assert isinstance(
        bp1_result, AutoACMGCriteria
    ), "The BP1 result should be of type AutoACMGCriteria."
    assert (
        bp1_result.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should be NotApplicable for this gene."
    assert "BP1 is not applicable for the gene." in bp1_result.summary


def test_predict_bp7_threshold_adjustment(brain_malformations_predictor, auto_acmg_data):
    """Test that the PhyloP100 threshold is correctly adjusted for BP7."""
    auto_acmg_data.thresholds.phyloP100 = 0.5  # Initial threshold value

    # Call predict_bp7 method
    result = brain_malformations_predictor.predict_bp7(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    # Check that the threshold was adjusted
    assert (
        auto_acmg_data.thresholds.phyloP100 == 0.1
    ), "The PhyloP100 threshold should be adjusted to 0.1 for BP7."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(DefaultSeqVarPredictor, "predict_bp7")
def test_predict_bp7_fallback_to_default(
    mock_super_predict_bp7, brain_malformations_predictor, auto_acmg_data
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
    result = brain_malformations_predictor.predict_bp7(
        brain_malformations_predictor.seqvar, auto_acmg_data
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


def test_predict_pp3bp4_pp3_not_applicable(brain_malformations_predictor, auto_acmg_data):
    pp3, bp4 = brain_malformations_predictor.predict_pp3bp4(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    assert pp3.name == "PP3"
    assert pp3.prediction == AutoACMGPrediction.NotApplicable
    assert pp3.strength == AutoACMGStrength.PathogenicSupporting
    assert "PP3 is not applicable" in pp3.summary


@pytest.mark.parametrize(
    "variant_type, spliceai_scores, expected_bp4",
    [
        ("synonymous_variant", [0.0, 0.0, 0.0, 0.0], AutoACMGPrediction.Applicable),
        ("intron_variant", [0.0, 0.0, 0.0, 0.0], AutoACMGPrediction.Applicable),
        ("5_prime_UTR_variant", [0.0, 0.0, 0.0, 0.0], AutoACMGPrediction.Applicable),
        ("synonymous_variant", [0.5, 0.0, 0.0, 0.0], AutoACMGPrediction.NotApplicable),
        # ("missense_variant", [0.0, 0.0, 0.0, 0.0], AutoACMGPrediction.NotMet),
    ],
)
def test_predict_pp3bp4_bp4_scenarios(
    brain_malformations_predictor,
    auto_acmg_data,
    variant_type,
    spliceai_scores,
    expected_bp4,
):
    auto_acmg_data.consequence = MagicMock(mehari=[variant_type])
    auto_acmg_data.scores.cadd.spliceAI_acceptor_gain = spliceai_scores[0]
    auto_acmg_data.scores.cadd.spliceAI_acceptor_loss = spliceai_scores[1]
    auto_acmg_data.scores.cadd.spliceAI_donor_gain = spliceai_scores[2]
    auto_acmg_data.scores.cadd.spliceAI_donor_loss = spliceai_scores[3]

    _, bp4 = brain_malformations_predictor.predict_pp3bp4(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    assert bp4.name == "BP4"
    assert bp4.prediction == expected_bp4
    assert bp4.strength == AutoACMGStrength.BenignSupporting
    assert "BP4 evaluation based on splicing predictions" in bp4.summary


def test_predict_pp3bp4_bp4_spliceai_details(brain_malformations_predictor, auto_acmg_data):
    auto_acmg_data.consequence = MagicMock(mehari=["synonymous_variant"])
    auto_acmg_data.scores.cadd.spliceAI_acceptor_gain = 0.1
    auto_acmg_data.scores.cadd.spliceAI_acceptor_loss = 0.2
    auto_acmg_data.scores.cadd.spliceAI_donor_gain = 0.3
    auto_acmg_data.scores.cadd.spliceAI_donor_loss = 0.4

    _, bp4 = brain_malformations_predictor.predict_pp3bp4(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    assert "SpliceAI scores: 0.1, 0.2, 0.3, 0.4" in bp4.summary


@pytest.mark.skip(reason="Should pass")
def test_predict_pp3bp4_bp4_non_qualifying_variant(brain_malformations_predictor, auto_acmg_data):
    auto_acmg_data.consequence = MagicMock(mehari=["missense_variant"])

    _, bp4 = brain_malformations_predictor.predict_pp3bp4(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    assert bp4.prediction == AutoACMGPrediction.NotApplicable
    assert "Variant type does not qualify for BP4 evaluation" in bp4.summary


@patch.object(BrainMalformationsPredictor, "_is_synonymous_variant")
@patch.object(BrainMalformationsPredictor, "_is_intron_variant")
@patch.object(BrainMalformationsPredictor, "_is_utr_variant")
@patch.object(BrainMalformationsPredictor, "_affect_spliceAI")
def test_predict_pp3bp4_method_calls(
    mock_affect_spliceAI,
    mock_is_utr,
    mock_is_intron,
    mock_is_synonymous,
    brain_malformations_predictor,
    auto_acmg_data,
):
    mock_is_synonymous.return_value = False
    mock_is_intron.return_value = False
    mock_is_utr.return_value = True
    mock_affect_spliceAI.return_value = False

    brain_malformations_predictor.predict_pp3bp4(
        brain_malformations_predictor.seqvar, auto_acmg_data
    )

    mock_is_synonymous.assert_called_once()
    mock_is_intron.assert_called_once()
    mock_is_utr.assert_called_once()
    mock_affect_spliceAI.assert_called_once()
