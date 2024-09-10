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
from src.vcep import ACADVLPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="17", pos=100, delete="A", insert="T")


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


@pytest.fixture
def acadvl_predictor(seqvar, auto_acmg_data):
    result = MagicMock()  # Mocking the AutoACMGResult object, which will be passed along
    return ACADVLPredictor(seqvar=seqvar, result=result, config=MagicMock())


def test_predict_pm1_in_critical_region(acadvl_predictor, auto_acmg_data):
    """Test when variant falls within a critical region for ACADVL."""
    auto_acmg_data.prot_pos = 220  # Set protein position within a critical region (214-223)
    result = acadvl_predictor.predict_pm1(acadvl_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for critical region variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "falls within a critical region" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_outside_critical_region(acadvl_predictor, auto_acmg_data):
    """Test when variant does not fall within any critical region for ACADVL."""
    auto_acmg_data.prot_pos = 300  # Set protein position outside all critical regions
    result = acadvl_predictor.predict_pm1(acadvl_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met for non-critical region variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not fall within any critical region" in result.summary
    ), "The summary should indicate no critical region."


def test_predict_pm1_edge_case_start_boundary(acadvl_predictor, auto_acmg_data):
    """Test when variant falls exactly on the start boundary of a critical region."""
    auto_acmg_data.prot_pos = 214  # Start boundary of the first critical region (214-223)
    result = acadvl_predictor.predict_pm1(acadvl_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met when on the start boundary of a critical region."
    assert (
        "falls within a critical region" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_edge_case_end_boundary(acadvl_predictor, auto_acmg_data):
    """Test when variant falls exactly on the end boundary of a critical region."""
    auto_acmg_data.prot_pos = 223  # End boundary of the first critical region (214-223)
    result = acadvl_predictor.predict_pm1(acadvl_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met when on the end boundary of a critical region."
    assert (
        "falls within a critical region" in result.summary
    ), "The summary should indicate the critical region."


def test_bs2_not_applicable_acadvl(acadvl_predictor, auto_acmg_data):
    """Test BS2 is not applicable for ACADVL as overridden."""
    result = acadvl_predictor._bs2_not_applicable(auto_acmg_data)
    assert result is True, "BS2 should always be not applicable for ACADVL."


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
def test_predict_pm2ba1bs1bs2(mock_super_method, acadvl_predictor, auto_acmg_data, seqvar):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = acadvl_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.001
    assert auto_acmg_data.thresholds.ba1_benign == 0.007
    assert auto_acmg_data.thresholds.bs1_benign == 0.0035

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_bp3_not_applicable_acadvl(acadvl_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = acadvl_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable for ACADVL."


def test_predict_pp2bp1(acadvl_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for ACADVL predictor."""

    # Call the method under test
    pp2_result, bp1_result = acadvl_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


def test_predict_pp3bp4_missense(acadvl_predictor, seqvar, auto_acmg_data):
    """Test predict_pp3bp4 for a missense variant with high REVEL score."""
    auto_acmg_data.consequence = MagicMock(cadd={"missense": True}, mehari=["missense_variant"])
    auto_acmg_data.scores.dbnsfp.revel = 0.8

    pp3, bp4 = acadvl_predictor.predict_pp3bp4(seqvar, auto_acmg_data)

    assert pp3.prediction == AutoACMGPrediction.Applicable
    assert bp4.prediction == AutoACMGPrediction.NotApplicable
    assert "REVEL score 0.8 > 0.75, PP3 met" in pp3.summary


def test_predict_pp3bp4_missense_low_revel(acadvl_predictor, seqvar, auto_acmg_data):
    """Test predict_pp3bp4 for a missense variant with low REVEL score."""
    auto_acmg_data.consequence = MagicMock(cadd={"missense": True}, mehari=["missense_variant"])
    auto_acmg_data.scores.dbnsfp.revel = 0.4

    pp3, bp4 = acadvl_predictor.predict_pp3bp4(seqvar, auto_acmg_data)

    assert pp3.prediction == AutoACMGPrediction.NotApplicable
    assert bp4.prediction == AutoACMGPrediction.Applicable
    assert "REVEL score 0.4 < 0.5, BP4 met" in bp4.summary


def test_predict_pp3bp4_inframe_indel(acadvl_predictor, seqvar, auto_acmg_data):
    """Test predict_pp3bp4 for an in-frame indel."""
    auto_acmg_data.consequence = MagicMock(cadd={"inframe": True}, mehari=["inframe_deletion"])
    auto_acmg_data.scores.dbnsfp.provean = -3.0
    auto_acmg_data.scores.dbnsfp.mutationTaster = 0.6

    pp3, bp4 = acadvl_predictor.predict_pp3bp4(seqvar, auto_acmg_data)

    assert pp3.prediction == AutoACMGPrediction.Applicable
    assert bp4.prediction == AutoACMGPrediction.NotApplicable


@pytest.mark.skip(reason="Need to fix")
def test_predict_pp3bp4_splice(acadvl_predictor, seqvar, auto_acmg_data):
    """Test predict_pp3bp4 for a splice variant."""
    auto_acmg_data.consequence = MagicMock(cadd="splice", mehari=["splice_variant"])
    auto_acmg_data.scores.dbscsnv.ada = 0.8
    auto_acmg_data.scores.dbscsnv.rf = 0.9
    auto_acmg_data.scores.cadd.ada = 0.8
    auto_acmg_data.scores.cadd.rf = 0.9
    auto_acmg_data.thresholds.ada = 0.5
    auto_acmg_data.thresholds.rf = 0.5

    pp3, bp4 = acadvl_predictor.predict_pp3bp4(seqvar, auto_acmg_data)

    assert pp3.prediction == AutoACMGPrediction.Applicable
    assert bp4.prediction == AutoACMGPrediction.NotApplicable


def test_predict_pp3bp4_other(acadvl_predictor, seqvar, auto_acmg_data):
    """Test predict_pp3bp4 for a variant that doesn't meet any criteria."""
    auto_acmg_data.consequence = MagicMock(cadd={}, mehari=["other_variant"])
    pp3, bp4 = acadvl_predictor.predict_pp3bp4(seqvar, auto_acmg_data)

    assert pp3.prediction == AutoACMGPrediction.NotApplicable
    assert bp4.prediction == AutoACMGPrediction.NotApplicable


@patch.object(ACADVLPredictor, "_is_pathogenic_splicing")
@patch.object(ACADVLPredictor, "_is_benign_splicing")
@pytest.mark.skip(reason="Need to fix")
def test_predict_pp3bp4_splice_methods(
    mock_benign, mock_pathogenic, acadvl_predictor, seqvar, auto_acmg_data
):
    """Test that _is_pathogenic_splicing and _is_benign_splicing are called for splice variants."""
    auto_acmg_data.consequence = MagicMock(cadd="splice", mehari=["splice_variant"])
    auto_acmg_data.scores.dbscsnv.ada = 0.8
    auto_acmg_data.scores.dbscsnv.rf = 0.9
    auto_acmg_data.scores.cadd.ada = 0.8
    auto_acmg_data.scores.cadd.rf = 0.9
    auto_acmg_data.thresholds.ada = 0.6
    auto_acmg_data.thresholds.rf = 0.7

    mock_pathogenic.return_value = True
    mock_benign.return_value = False

    pp3, bp4 = acadvl_predictor.predict_pp3bp4(seqvar, auto_acmg_data)

    mock_pathogenic.assert_called_once()
    mock_benign.assert_called_once()
    assert pp3.prediction == AutoACMGPrediction.Applicable
    assert bp4.prediction == AutoACMGPrediction.NotApplicable
