from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
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
    return AutoACMGData()


def test_predict_pm1_in_critical_domain(congenital_myopathies_predictor, auto_acmg_data):
    """Test when the variant falls within a critical domain for RYR1."""
    auto_acmg_data.prot_pos = 4850  # Within the critical domain (4800-4950) for RYR1
    auto_acmg_data.hgnc_id = "HGNC:10483"  # RYR1 gene
    result = congenital_myopathies_predictor.predict_pm1(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
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
        result.prediction == AutoACMGPrediction.NotMet
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


@patch("src.vcep.congenital_myopathies.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, congenital_myopathies_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method."""
    auto_acmg_data.hgnc_id = "HGNC:99999"  # Not in the PM1_CLUSTER mapping
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = congenital_myopathies_predictor.predict_pm1(
        congenital_myopathies_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
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
        result.prediction == AutoACMGPrediction.Met
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
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met when on the end boundary of a critical domain."
    assert (
        "falls within a critical domain" in result.summary
    ), "The summary should indicate the critical domain."


@patch.object(
    DefaultPredictor,
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
        pp2_result.prediction == AutoACMGPrediction.Met
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
        pp2_result.prediction == AutoACMGPrediction.NotMet
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
