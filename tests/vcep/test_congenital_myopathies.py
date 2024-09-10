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
from src.vcep import CongenitalMyopathiesPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="19", pos=100, delete="A", insert="T")


@pytest.fixture
def congenital_myopathies_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return CongenitalMyopathiesPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


@patch.object(
    DefaultSeqVarPredictor,
    "predict_pvs1",
    return_value=AutoACMGCriteria(
        name="PVS1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicVeryStrong,
        summary="Default behavior.",
    ),
)
def test_predict_pvs1_not_applicable_for_specific_gene(
    mock_super_predict_pvs1, congenital_myopathies_predictor, seqvar, auto_acmg_data
):
    # Set the HGNC ID to one that makes PVS1 not applicable
    auto_acmg_data.hgnc_id = "HGNC:2974"
    result = congenital_myopathies_predictor.predict_pvs1(seqvar, auto_acmg_data)

    # Verify the outcome is as expected for the specific HGNC ID
    assert result.name == "PVS1"
    assert result.prediction == AutoACMGPrediction.NotApplicable
    assert result.strength == AutoACMGStrength.PathogenicVeryStrong
    assert result.summary == "PVS1 is not applicable for the gene."
    mock_super_predict_pvs1.assert_not_called()  # Ensure superclass method is not called


@patch.object(
    DefaultSeqVarPredictor,
    "predict_pvs1",
    return_value=AutoACMGCriteria(
        name="PVS1",
        prediction=AutoACMGPrediction.Applicable,
        strength=AutoACMGStrength.PathogenicVeryStrong,
        summary="Superclass default behavior.",
    ),
)
def test_predict_pvs1_calls_superclass_when_not_specific_gene(
    mock_super_predict_pvs1, congenital_myopathies_predictor, seqvar, auto_acmg_data
):
    # Set the HGNC ID to one not affecting the PVS1 applicability
    auto_acmg_data.hgnc_id = "HGNC:Random"
    result = congenital_myopathies_predictor.predict_pvs1(seqvar, auto_acmg_data)

    # Verify that the superclass method is called
    mock_super_predict_pvs1.assert_called_once_with(seqvar, auto_acmg_data)
    # Check the response from the superclass method
    assert result.name == "PVS1"
    assert result.prediction == AutoACMGPrediction.Applicable
    assert result.strength == AutoACMGStrength.PathogenicVeryStrong
    assert result.summary == "Superclass default behavior."


def test_predict_pm1_in_critical_domain(congenital_myopathies_predictor, auto_acmg_data):
    """Test when the variant falls within a critical domain for RYR1."""
    auto_acmg_data.prot_pos = 4850  # Within the critical domain (4800-4950) for RYR1
    auto_acmg_data.hgnc_id = "HGNC:10483"  # RYR1 gene
    result = congenital_myopathies_predictor.predict_pm1(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for a variant in a critical domain."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."


def test_predict_pm1_outside_critical_domain(congenital_myopathies_predictor, auto_acmg_data):
    """Test when the variant does not fall within any critical domain for RYR1."""
    auto_acmg_data.prot_pos = 5000  # Outside the critical domain (4800-4950) for RYR1
    auto_acmg_data.hgnc_id = "HGNC:10483"  # RYR1 gene
    result = congenital_myopathies_predictor.predict_pm1(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met for a variant outside any critical domain."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "Variant does not fall within a critical domain." in result.summary
    ), "The summary should indicate no critical domain."


def test_predict_pm1_not_applicable(congenital_myopathies_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for NEB, ACTA1, DNM2, MTM1."""
    auto_acmg_data.hgnc_id = "HGNC:7720"  # NEB gene
    result = congenital_myopathies_predictor.predict_pm1(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for NEB."
    assert (
        "Not applicable for NEB" in result.summary
    ), "The summary should indicate that PM1 is not applicable for NEB."


@patch("src.vcep.congenital_myopathies.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, congenital_myopathies_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method."""
    auto_acmg_data.hgnc_id = "HGNC:99999"  # Not in the PM1_CLUSTER mapping
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = congenital_myopathies_predictor.predict_pm1(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should be from the default fallback."


def test_predict_pm1_edge_case_start_boundary(congenital_myopathies_predictor, auto_acmg_data):
    """Test when variant falls exactly on the start boundary of a critical domain."""
    auto_acmg_data.prot_pos = 4800  # Start boundary of the RYR1 critical domain (4800-4950)
    auto_acmg_data.hgnc_id = "HGNC:10483"  # RYR1 gene
    result = congenital_myopathies_predictor.predict_pm1(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met when on the start boundary of a critical domain."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."


def test_predict_pm1_edge_case_end_boundary(congenital_myopathies_predictor, auto_acmg_data):
    """Test when variant falls exactly on the end boundary of a critical domain."""
    auto_acmg_data.prot_pos = 4950  # End boundary of the RYR1 critical domain (4800-4950)
    auto_acmg_data.hgnc_id = "HGNC:10483"  # RYR1 gene
    result = congenital_myopathies_predictor.predict_pm1(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met when on the end boundary of a critical domain."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."


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
    "hgnc_id,expected_pm2,expected_ba1,expected_bs1",
    [
        ("HGNC:7720", 0.0000001, 0.00559, 0.000237),
        ("HGNC:129", 0.0000001, 0.0000781, 0.00000781),
        ("HGNC:2974", 0.0000001, 0.0000015, 0.00000015),
        ("HGNC:7448", 0.0000001, 0.000016, 0.0000016),
        ("HGNC:10483", 0.0000001, 0.0000486, 0.00000486),
    ],
)
def test_predict_pm2ba1bs1bs2_with_varied_thresholds(
    mock_super_method,
    congenital_myopathies_predictor,
    auto_acmg_data,
    seqvar,
    hgnc_id,
    expected_pm2,
    expected_ba1,
    expected_bs1,
):
    # Setup
    auto_acmg_data.hgnc_id = hgnc_id

    # Method call
    congenital_myopathies_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Validate thresholds are set correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == expected_pm2
    assert auto_acmg_data.thresholds.ba1_benign == expected_ba1
    assert auto_acmg_data.thresholds.bs1_benign == expected_bs1

    # Check that the superclass method was called with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Reset mock for the next iteration
    mock_super_method.reset_mock()


def test_bp3_not_applicable(congenital_myopathies_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = congenital_myopathies_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


def test_predict_pp2bp1_missense(congenital_myopathies_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 when the variant is a missense change in a relevant gene."""
    auto_acmg_data.hgnc_id = "HGNC:129"  # Relevant gene
    auto_acmg_data.consequence.cadd = "missense_variant"  # Missense change

    # Call the method under test
    pp2_result, bp1_result = congenital_myopathies_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    # Check PP2 result
    assert isinstance(
        pp2_result, AutoACMGCriteria
    ), "The PP2 result should be of type AutoACMGCriteria."
    assert (
        pp2_result.prediction == AutoACMGPrediction.Applicable
    ), "PP2 should be Met for a missense variant in a relevant gene."
    assert "PP2 is met for HGNC:129 as the variant is a missense change." in pp2_result.summary

    # Check BP1 result
    assert isinstance(
        bp1_result, AutoACMGCriteria
    ), "The BP1 result should be of type AutoACMGCriteria."
    assert (
        bp1_result.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should be NotApplicable for this gene."
    assert "BP1 is not applicable for the gene." in bp1_result.summary


def test_predict_pp2bp1_not_missense(congenital_myopathies_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 when the variant is not a missense change in a relevant gene."""
    auto_acmg_data.hgnc_id = "HGNC:129"  # Relevant gene
    auto_acmg_data.consequence.cadd = "synonymous_variant"  # Not a missense change

    # Call the method under test
    pp2_result, bp1_result = congenital_myopathies_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    # Check PP2 result
    assert isinstance(
        pp2_result, AutoACMGCriteria
    ), "The PP2 result should be of type AutoACMGCriteria."
    assert (
        pp2_result.prediction == AutoACMGPrediction.NotApplicable
    ), "PP2 should not be met for a non-missense variant in a relevant gene."
    assert (
        "PP2 is not met for HGNC:129 as the variant is not a missense change." in pp2_result.summary
    )

    # Check BP1 result
    assert isinstance(
        bp1_result, AutoACMGCriteria
    ), "The BP1 result should be of type AutoACMGCriteria."
    assert (
        bp1_result.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should be NotApplicable for this gene."
    assert "BP1 is not applicable for the gene." in bp1_result.summary


def test_verify_pp3bp4_revel_thresholds(congenital_myopathies_predictor, auto_acmg_data):
    """Test that REVEL thresholds are correctly set for PP3/BP4 prediction."""
    prediction, comment = congenital_myopathies_predictor.verify_pp3bp4(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    assert auto_acmg_data.thresholds.revel_pathogenic == 0.7
    assert auto_acmg_data.thresholds.revel_benign == 0.15


@pytest.mark.parametrize(
    "revel_score, expected_pp3, expected_bp4",
    [
        (0.8, True, False),
        (0.5, False, False),
        (0.1, False, True),
    ],
)
def test_verify_pp3bp4_revel_scenarios(
    congenital_myopathies_predictor,
    auto_acmg_data,
    revel_score,
    expected_pp3,
    expected_bp4,
):
    auto_acmg_data.scores.dbnsfp.revel = revel_score

    prediction, comment = congenital_myopathies_predictor.verify_pp3bp4(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 == expected_pp3
    assert prediction.BP4 == expected_bp4


@patch.object(CongenitalMyopathiesPredictor, "_affect_spliceAI")
def test_verify_pp3bp4_splicing(
    mock_affect_spliceAI, congenital_myopathies_predictor, auto_acmg_data
):
    mock_affect_spliceAI.return_value = True
    auto_acmg_data.scores.dbnsfp.revel = 0.5  # Between benign and pathogenic thresholds

    prediction, comment = congenital_myopathies_predictor.verify_pp3bp4(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 is True
    assert prediction.BP4 is False


@patch.object(CongenitalMyopathiesPredictor, "_affect_spliceAI")
@pytest.mark.skip("Fix it")
def test_verify_pp3bp4_no_splicing_effect(
    mock_affect_spliceAI, congenital_myopathies_predictor, auto_acmg_data
):
    mock_affect_spliceAI.return_value = False
    auto_acmg_data.scores.dbnsfp.revel = 0.5  # Between benign and pathogenic thresholds

    prediction, comment = congenital_myopathies_predictor.verify_pp3bp4(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 is False
    assert prediction.BP4 is True


@pytest.mark.skip("Fix it")
def test_verify_pp3bp4_error_handling(congenital_myopathies_predictor, auto_acmg_data):
    # Simulate an error condition
    auto_acmg_data.scores.dbnsfp.revel = None

    prediction, comment = congenital_myopathies_predictor.verify_pp3bp4(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    assert prediction is None
    assert "An error occurred during prediction" in comment


@patch.object(CongenitalMyopathiesPredictor, "_is_pathogenic_score")
@patch.object(CongenitalMyopathiesPredictor, "_is_benign_score")
@patch.object(CongenitalMyopathiesPredictor, "_affect_spliceAI")
def test_verify_pp3bp4_method_calls(
    mock_affect_spliceAI,
    mock_is_benign_score,
    mock_is_pathogenic_score,
    congenital_myopathies_predictor,
    auto_acmg_data,
):
    mock_is_pathogenic_score.return_value = False
    mock_is_benign_score.return_value = False
    mock_affect_spliceAI.return_value = True

    congenital_myopathies_predictor.verify_pp3bp4(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    mock_is_pathogenic_score.assert_called_once()
    mock_is_benign_score.assert_called_once()
    mock_affect_spliceAI.assert_called()


def test_verify_pp3bp4_spliceai_thresholds(congenital_myopathies_predictor, auto_acmg_data):
    prediction, comment = congenital_myopathies_predictor.verify_pp3bp4(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.05
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.05
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.05
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.05


def test_verify_pp3bp4_benign_spliceai_thresholds(congenital_myopathies_predictor, auto_acmg_data):
    prediction, comment = congenital_myopathies_predictor.verify_pp3bp4(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.05
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.05
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.05
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.05
