from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import (
    PS1PM5,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
)
from src.defs.exceptions import AutoAcmgBaseException
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor
from src.vcep import ENIGMAPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="17", pos=100, delete="A", insert="T")


@pytest.fixture
def enigma_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return ENIGMAPredictor(seqvar=seqvar, result=result)


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


# --------------- PS1 & PM5 ---------------


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_overrides(mock_super_verify, enigma_predictor, seqvar, auto_acmg_data):
    """Test that the overridden verify_ps1pm5 method in ENIGMAPredictor works correctly."""
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
    prediction, comment = enigma_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that the splicing effect leads to overriding PS1 and PM5 as not applicable
    assert not prediction.PS1, "PS1 should be marked as not applicable due to splicing effect."
    assert not prediction.PM5, "PM5 should be marked as not applicable due to splicing effect."
    assert "Variant affects splicing" in comment, "Comment should note the splicing effect."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_no_splicing_effect(
    mock_super_verify, enigma_predictor, seqvar, auto_acmg_data
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
    prediction, comment = enigma_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that PS1 and PM5 remain applicable
    assert prediction.PS1, "PS1 should remain applicable when there's no splicing effect."
    assert prediction.PM5, "PM5 should remain applicable when there's no splicing effect."
    assert "Initial evaluation" in comment, "Comment should reflect the initial evaluation."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_non_missense(mock_super_verify, enigma_predictor, seqvar, auto_acmg_data):
    """Test that PS1 and PM5 remain as per superclass for non-missense variants."""
    # Set up the mock to return PS1 and PM5 as not applicable initially
    mock_super_verify.return_value = (
        PS1PM5(PS1=False, PM5=False),
        "Not applicable for non-missense",
    )

    # Setup the data with a non-missense variant
    auto_acmg_data.consequence = MagicMock(mehari=["synonymous_variant"])

    # Run the method under test
    prediction, comment = enigma_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that PS1 and PM5 remain as per superclass evaluation
    assert not prediction.PS1, "PS1 should remain not applicable for non-missense variants."
    assert not prediction.PM5, "PM5 should remain not applicable for non-missense variants."
    assert (
        "Not applicable for non-missense" in comment
    ), "Comment should reflect superclass evaluation."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)


# --------------- PM1 ---------------


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
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = enigma_predictor.predict_pm1(enigma_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
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


# --------------- PM2, BA1, BS1, BS2 ---------------


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


# --------------- PM4 & BP3 ---------------


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


# --------------- PP2 & BP1 ---------------


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
        bp1.prediction == AutoACMGPrediction.Applicable
    ), "BP1 should be Met as the criteria are satisfied."
    assert (
        bp1.strength == AutoACMGStrength.BenignSupporting
    ), "BP1 strength should be BenignSupporting."
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
        bp1.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should not be Met as the variant is in an important domain."


def test_predict_pp2bp1_synonymous_splice_affecting(enigma_predictor, seqvar, auto_acmg_data):
    """Test PP2 and BP1 where the variant is synonymous and predicted to affect splicing."""
    auto_acmg_data.consequence = MagicMock(mehari=["synonymous_variant"], cadd=None)
    auto_acmg_data.prot_pos = 2000  # Outside any important domain

    with patch.object(ENIGMAPredictor, "_spliceai_impact", return_value=True):
        pp2, bp1 = enigma_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert pp2.prediction == AutoACMGPrediction.NotApplicable, "PP2 should be NotApplicable."
    assert (
        bp1.prediction == AutoACMGPrediction.NotApplicable
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
        bp1.prediction == AutoACMGPrediction.Applicable
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
        bp1.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should not be Met due to variant type being not relevant."
    assert (
        "not synonymous, missense or inframe indel" in bp1.summary
    ), "The summary should indicate the variant type."


# --------------- PP3 & BP4 ---------------


def test_predict_pp3bp4_thresholds(enigma_predictor, auto_acmg_data):
    """Test that thresholds are correctly set for PP3/BP4 prediction."""
    enigma_predictor.predict_pp3bp4(enigma_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.bayesDel_noAF_pathogenic == 0.521
    assert auto_acmg_data.thresholds.bayesDel_noAF_benign == -0.476
    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.1
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.1
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.1
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.1


@pytest.mark.parametrize(
    "is_missense, is_in_domain, bayesdel_score, spliceai_score, expected_pp3, expected_bp4",
    [
        (
            True,
            True,
            0.3,
            0.1,
            True,
            False,
        ),  # Missense in domain, high BayesDel, low SpliceAI
        (
            True,
            True,
            0.1,
            0.1,
            False,
            True,
        ),  # Missense in domain, low BayesDel, low SpliceAI
        (True, False, 0.3, 0.1, False, False),  # Missense not in domain, high BayesDel
        (False, False, 0.1, 0.3, True, False),  # Non-missense, high SpliceAI
        (False, False, 0.1, 0.05, False, True),  # Non-missense, low SpliceAI
    ],
)
def test_predict_pp3bp4_scenarios(
    enigma_predictor,
    auto_acmg_data,
    is_missense,
    is_in_domain,
    bayesdel_score,
    spliceai_score,
    expected_pp3,
    expected_bp4,
):
    with (
        patch.object(ENIGMAPredictor, "_is_missense_variant", return_value=is_missense),
        patch.object(ENIGMAPredictor, "_is_inframe_indel", return_value=False),
        patch.object(ENIGMAPredictor, "_in_important_domain", return_value=is_in_domain),
        patch.object(ENIGMAPredictor, "_is_intronic", return_value=not is_missense),
        patch.object(ENIGMAPredictor, "_is_synonymous_variant", return_value=not is_missense),
    ):
        auto_acmg_data.scores.dbnsfp.bayesDel_noAF = bayesdel_score
        auto_acmg_data.scores.cadd.spliceAI_acceptor_gain = spliceai_score

        pp3_result, bp4_result = enigma_predictor.predict_pp3bp4(
            enigma_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == (
            AutoACMGPrediction.Applicable if expected_pp3 else AutoACMGPrediction.NotApplicable
        )
        assert bp4_result.prediction == (
            AutoACMGPrediction.Applicable if expected_bp4 else AutoACMGPrediction.NotApplicable
        )


def test_predict_pp3bp4_missense_in_domain_high_bayesdel(enigma_predictor, auto_acmg_data):
    with (
        patch.object(ENIGMAPredictor, "_is_missense_variant", return_value=True),
        patch.object(ENIGMAPredictor, "_in_important_domain", return_value=True),
        patch.object(ENIGMAPredictor, "_is_pathogenic_score", return_value=True),
    ):
        pp3_result, bp4_result = enigma_predictor.predict_pp3bp4(
            enigma_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == AutoACMGPrediction.Applicable
        assert "BayesDel_noAF score" in pp3_result.summary
        assert bp4_result.prediction == AutoACMGPrediction.NotApplicable


def test_predict_pp3bp4_missense_in_domain_low_bayesdel_no_splice(enigma_predictor, auto_acmg_data):
    with (
        patch.object(ENIGMAPredictor, "_is_missense_variant", return_value=True),
        patch.object(ENIGMAPredictor, "_in_important_domain", return_value=True),
        patch.object(ENIGMAPredictor, "_is_benign_score", return_value=True),
        patch.object(ENIGMAPredictor, "_affect_spliceAI", return_value=False),
    ):
        pp3_result, bp4_result = enigma_predictor.predict_pp3bp4(
            enigma_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == AutoACMGPrediction.NotApplicable
        assert bp4_result.prediction == AutoACMGPrediction.Applicable
        assert "BayesDel_noAF score" in bp4_result.summary


def test_predict_pp3bp4_splice_effect(enigma_predictor, auto_acmg_data):
    with (
        patch.object(ENIGMAPredictor, "_is_missense_variant", return_value=True),
        patch.object(ENIGMAPredictor, "_affect_spliceAI", return_value=True),
    ):
        pp3_result, bp4_result = enigma_predictor.predict_pp3bp4(
            enigma_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == AutoACMGPrediction.Applicable
        assert "SpliceAI ≥0.2" in pp3_result.summary
        assert bp4_result.prediction == AutoACMGPrediction.NotApplicable


def test_predict_pp3bp4_intronic_no_splice_effect(enigma_predictor, auto_acmg_data):
    with (
        patch.object(ENIGMAPredictor, "_is_intronic", return_value=True),
        patch.object(ENIGMAPredictor, "_affect_spliceAI", return_value=False),
    ):
        pp3_result, bp4_result = enigma_predictor.predict_pp3bp4(
            enigma_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == AutoACMGPrediction.NotApplicable
        assert bp4_result.prediction == AutoACMGPrediction.Applicable
        assert "SpliceAI ≥0.2" in bp4_result.summary


def test_predict_pp3bp4_strength(enigma_predictor, auto_acmg_data):
    pp3_result, bp4_result = enigma_predictor.predict_pp3bp4(
        enigma_predictor.seqvar, auto_acmg_data
    )

    assert pp3_result.strength == AutoACMGStrength.PathogenicSupporting
    assert bp4_result.strength == AutoACMGStrength.BenignSupporting


def test_predict_pp3bp4_no_criteria_met(enigma_predictor, auto_acmg_data):
    with (
        patch.object(ENIGMAPredictor, "_is_missense_variant", return_value=False),
        patch.object(ENIGMAPredictor, "_is_intronic", return_value=False),
        patch.object(ENIGMAPredictor, "_is_synonymous_variant", return_value=False),
        patch.object(ENIGMAPredictor, "_affect_spliceAI", return_value=False),
    ):
        pp3_result, bp4_result = enigma_predictor.predict_pp3bp4(
            enigma_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == AutoACMGPrediction.NotApplicable
        assert bp4_result.prediction == AutoACMGPrediction.NotApplicable
        assert "PP3 criteria not met." in pp3_result.summary
        assert "BP4 criteria not met." in bp4_result.summary


# -------------- BP7 --------------


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
