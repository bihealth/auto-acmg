from unittest.mock import MagicMock, patch

import pytest

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import (
    AlleleCondition,
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
)
from src.defs.exceptions import MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep.scid import SCIDPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="10", pos=100, delete="A", insert="T")


@pytest.fixture
def scid_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return SCIDPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_in_moderate_domain(scid_predictor, auto_acmg_data):
    """Test when the variant falls within a moderate domain for a SCID gene."""
    auto_acmg_data.hgnc_id = "HGNC:12765"  # FOXN1 gene
    auto_acmg_data.prot_pos = 300  # Within the DNA binding forkhead domain (270-367)
    result = scid_predictor.predict_pm1(scid_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for critical domain variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "critical residue in HGNC:12765" in result.summary
    ), "The summary should indicate the critical domain."


def test_predict_pm1_in_supporting_domain(scid_predictor, auto_acmg_data):
    """Test when the variant falls within a supporting domain for a SCID gene."""
    auto_acmg_data.hgnc_id = "HGNC:9831"  # RAG1 gene
    auto_acmg_data.prot_pos = 600  # Within the core domain (387-1011)
    result = scid_predictor.predict_pm1(scid_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for supporting domain variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "supporting region of HGNC:9831" in result.summary
    ), "The summary should indicate the supporting domain."


def test_predict_pm1_strong_residue(scid_predictor, auto_acmg_data):
    """Test when the variant affects a strong residue for a SCID gene."""
    auto_acmg_data.hgnc_id = "HGNC:6010"  # IL2RG gene
    auto_acmg_data.prot_pos = 62  # A conserved cysteine residue
    result = scid_predictor.predict_pm1(scid_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for strong residue variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicStrong
    ), "The strength should be PathogenicStrong."
    assert (
        "strong residue in HGNC:6010" in result.summary
    ), "The summary should indicate the strong residue."


def test_predict_pm1_not_applicable(scid_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for ADA, DCLRE1C, or IL7R."""
    auto_acmg_data.hgnc_id = "HGNC:186"  # ADA gene
    result = scid_predictor.predict_pm1(scid_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be applicable for ADA."
    assert "not applicable" in result.summary, "The summary should indicate non-applicability."


@patch("src.vcep.scid.DefaultPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, scid_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the SCID VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default fallback for PM1.",
    )

    result = scid_predictor.predict_pm1(scid_predictor.seqvar, auto_acmg_data)
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default fallback for PM1." in result.summary
    ), "The summary should indicate the default fallback."


@patch.object(
    SCIDPredictor,
    "_get_allele_cond",
    return_value=AlleleCondition.Recessive,
)
@patch.object(
    SCIDPredictor,
    "_get_control_af",
    return_value=MagicMock(bySex=MagicMock(overall=MagicMock(ac=10, nhomalt=6))),
)
@patch.object(SCIDPredictor, "_get_any_af", return_value=None)
def test_check_zyg_homozygous_positive(
    mock_get_any_af,
    mock_get_control_af,
    mock_get_allele_cond,
    scid_predictor,
    seqvar,
    auto_acmg_data,
):
    scid_predictor.comment_pm2ba1bs1bs2 = ""
    assert scid_predictor._check_zyg(seqvar, auto_acmg_data) == True
    assert (
        "The variant is in a recessive (homozygous) disorder."
        in scid_predictor.comment_pm2ba1bs1bs2
    )


@patch.object(
    SCIDPredictor,
    "_get_allele_cond",
    return_value=AlleleCondition.Recessive,
)
@patch.object(
    SCIDPredictor,
    "_get_control_af",
    return_value=MagicMock(bySex=MagicMock(overall=MagicMock(ac=10, nhomalt=4))),
)
@patch.object(SCIDPredictor, "_get_any_af", return_value=None)
def test_check_zyg_homozygous_negative(
    mock_get_any_af,
    mock_get_control_af,
    mock_get_allele_cond,
    scid_predictor,
    seqvar,
    auto_acmg_data,
):
    scid_predictor.comment_pm2ba1bs1bs2 = ""
    assert scid_predictor._check_zyg(seqvar, auto_acmg_data) == False


@patch.object(
    SCIDPredictor,
    "_get_allele_cond",
    return_value=AlleleCondition.Recessive,
)
@patch.object(SCIDPredictor, "_get_control_af", return_value=None)
@patch.object(SCIDPredictor, "_get_any_af", return_value=None)
def test_check_zyg_missing_data_raises_error(
    mock_get_any_af,
    mock_get_control_af,
    mock_get_allele_cond,
    scid_predictor,
    seqvar,
    auto_acmg_data,
):
    scid_predictor.comment_pm2ba1bs1bs2 = ""
    with pytest.raises(MissingDataError):
        scid_predictor._check_zyg(seqvar, auto_acmg_data)


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
        ("HGNC:12765", 0.00002412, 0.00447, 0.00141),
        ("HGNC:186", 0.0001742, 0.00721, 0.00161),
        ("HGNC:17642", 0.00003266, 0.00346, 0.00078),
        ("HGNC:6024", 0.00004129, 0.00566, 0.00126),
        ("HGNC:6193", 0.000115, 0.00447, 0.001),
        ("HGNC:9831", 0.000102, 0.00872, 0.00195),
        ("HGNC:9832", 0.0000588, 0.00872, 0.00195),
        ("HGNC:6010", 0.000124, 0.01110, 0.00249),
    ],
)
def test_predict_pm2ba1bs1bs2_gene_specific(
    mock_super_method,
    scid_predictor,
    auto_acmg_data,
    seqvar,
    hgnc_id,
    expected_pm2,
    expected_ba1,
    expected_bs1,
):
    # Set gene ID
    auto_acmg_data.hgnc_id = hgnc_id

    # Call the method under test
    scid_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds are set correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == expected_pm2
    assert auto_acmg_data.thresholds.ba1_benign == expected_ba1
    assert auto_acmg_data.thresholds.bs1_benign == expected_bs1

    # Validate that the superclass method was called correctly with modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Reset mock for the next iteration
    mock_super_method.reset_mock()


def test_is_conserved_scid_gene(scid_predictor, auto_acmg_data):
    """Test the _is_conserved method for SCID genes."""
    auto_acmg_data.hgnc_id = "HGNC:6024"  # IL7R gene

    # Call _is_conserved method
    result = scid_predictor._is_conserved(auto_acmg_data)

    # Check that the conservation is considered for FOXN1
    assert result is False, "The conservation check should be False for SCID genes."


def test_is_conserved_non_scid_gene(scid_predictor, auto_acmg_data):
    """Test the _is_conserved method for non-SCID genes."""
    auto_acmg_data.hgnc_id = "HGNC:6010"  # IL2RG gene

    # Call _is_conserved method
    result = scid_predictor._is_conserved(auto_acmg_data)

    # Check that conservation is ignored for non-SCID genes
    assert result is False, "The conservation check should be ignored for non-SCID genes."


def test_predict_pp2bp1(scid_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for SCID."""

    # Call the method under test
    pp2_result, bp1_result = scid_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


def test_predict_bp7_threshold_adjustment(scid_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted for SCID."""
    auto_acmg_data.thresholds.bp7_donor = 1  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = scid_predictor.predict_bp7(scid_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.bp7_donor == 7
    ), "The BP7 donor threshold should be adjusted to 7."
    assert (
        auto_acmg_data.thresholds.bp7_acceptor == 21
    ), "The BP7 acceptor threshold should be adjusted to 21."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(SCIDPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default(mock_super_predict_bp7, scid_predictor, auto_acmg_data):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = scid_predictor.predict_bp7(scid_predictor.seqvar, auto_acmg_data)

    # Verify the result and ensure the superclass method was called
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert result.prediction == AutoACMGPrediction.NotMet, "BP7 should return NotMet as mocked."
    assert (
        result.strength == AutoACMGStrength.BenignSupporting
    ), "The strength should be BenignSupporting."
    assert (
        "Default BP7 prediction fallback." in result.summary
    ), "The summary should indicate the fallback."
    assert mock_super_predict_bp7.called, "super().predict_bp7 should have been called."
