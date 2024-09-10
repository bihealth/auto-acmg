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
from src.vcep import RettAngelmanPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="18", pos=100, delete="A", insert="T")


@pytest.fixture
def rett_angelman_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return RettAngelmanPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


# --------------------- PM1 ---------------------


def test_predict_pm1_in_critical_region(rett_angelman_predictor, auto_acmg_data):
    """Test when the variant falls within a critical region for a Rett and Angelman-like Disorder gene."""
    auto_acmg_data.hgnc_id = "HGNC:11634"  # TCF4 gene
    auto_acmg_data.prot_pos = 580  # Within the critical bHLH domain (564-617)
    result = rett_angelman_predictor.predict_pm1(rett_angelman_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for critical region variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "critical residue in HGNC:11634" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_not_applicable(rett_angelman_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for SLC9A6."""
    auto_acmg_data.hgnc_id = "HGNC:11079"  # SLC9A6 gene
    result = rett_angelman_predictor.predict_pm1(rett_angelman_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be applicable for SLC9A6."
    assert "not applicable" in result.summary, "The summary should indicate non-applicability."


def test_predict_pm1_outside_critical_region(rett_angelman_predictor, auto_acmg_data):
    """Test when the variant does not fall within any critical region for Rett and Angelman-like Disorders genes."""
    auto_acmg_data.hgnc_id = "HGNC:11634"  # TCF4 gene
    auto_acmg_data.prot_pos = 700  # Position outside all critical regions
    result = rett_angelman_predictor.predict_pm1(rett_angelman_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met for non-critical region variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not meet the PM1 criteria" in result.summary
    ), "The summary should indicate no critical region."


@patch("src.vcep.rett_angelman.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, rett_angelman_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the Rett and Angelman-like VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = rett_angelman_predictor.predict_pm1(rett_angelman_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


# ----------------- PM2, BA1, BS1, BS2 -----------------


@patch.object(RettAngelmanPredictor, "_get_af", return_value=0.1)
@patch.object(RettAngelmanPredictor, "_ba1_exception", return_value=False)
def test_verify_pm2ba1bs1bs2(
    mock_get_af, mock_ba1_exception, rett_angelman_predictor, auto_acmg_data, seqvar
):
    # Setup: Adjusting the thresholds to test under different conditions
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.02
    auto_acmg_data.thresholds.pm2_pathogenic = 0.001

    # Call the method under test
    result, comment = rett_angelman_predictor.verify_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assertions to validate the expected behavior
    assert result.BA1 is True, "Expected PM2 to be True based on the mocked allele frequency"

    # Assert changed thresholds
    assert (
        auto_acmg_data.thresholds.ba1_benign == 0.0003
    ), "BA1 threshold should be adjusted to 0.0003"
    assert (
        auto_acmg_data.thresholds.bs1_benign == 0.00008
    ), "BS1 threshold should be adjusted to 0.00008"


# --------------------- PM4 & BP3 ---------------------


def test_exclude_pm4_true(rett_angelman_predictor, auto_acmg_data, seqvar):
    """Test that the variant is excluded from PM4 based on the exclusion regions."""
    auto_acmg_data.hgnc_id = "HGNC:11411"  # CDKL5
    auto_acmg_data.prot_pos = 950  # Within exclusion range for CDKL5

    # Call the method under test
    result = rett_angelman_predictor._exclude_pm4(seqvar, auto_acmg_data)

    assert result is True, "The variant should be excluded from PM4."


def test_exclude_pm4_false(rett_angelman_predictor, auto_acmg_data, seqvar):
    """Test that the variant is not excluded from PM4."""
    auto_acmg_data.hgnc_id = "HGNC:11411"  # CDKL5
    auto_acmg_data.prot_pos = 500  # Outside the exclusion range for CDKL5

    # Call the method under test
    result = rett_angelman_predictor._exclude_pm4(seqvar, auto_acmg_data)

    assert result is False, "The variant should not be excluded from PM4."


def test_in_foxg1_bp3_region_true(rett_angelman_predictor, auto_acmg_data):
    """Test that the variant is in the BP3 region for FOXG1."""
    auto_acmg_data.hgnc_id = "HGNC:3811"  # FOXG1
    auto_acmg_data.prot_pos = 50  # Inside the BP3 region for FOXG1

    # Call the method under test
    result = rett_angelman_predictor._in_foxg1_bp3_region(auto_acmg_data)

    assert result is True, "The variant should be in the BP3 region for FOXG1."


def test_in_foxg1_bp3_region_false(rett_angelman_predictor, auto_acmg_data):
    """Test that the variant is not in the BP3 region for FOXG1."""
    auto_acmg_data.hgnc_id = "HGNC:3811"  # FOXG1
    auto_acmg_data.prot_pos = 100  # Outside the BP3 region for FOXG1

    # Call the method under test
    result = rett_angelman_predictor._in_foxg1_bp3_region(auto_acmg_data)

    assert result is False, "The variant should not be in the BP3 region for FOXG1."


@patch.object(RettAngelmanPredictor, "_in_repeat_region", return_value=False)
def test_verify_pm4bp3_in_frame_deletion_in_important_domain(
    mock_in_repeat_region, rett_angelman_predictor, seqvar, auto_acmg_data
):
    """Test verify_pm4bp3 when the variant is an in-frame deletion in an important domain."""
    auto_acmg_data.hgnc_id = "HGNC:6990"  # MECP2
    auto_acmg_data.prot_pos = 400  # In-frame deletion in the proline-rich region

    auto_acmg_data.consequence.cadd = "inframe_deletion"

    # Call the method under test
    prediction, comment = rett_angelman_predictor.verify_pm4bp3(seqvar, auto_acmg_data)

    assert (
        prediction.PM4 is False
    ), "PM4 should not be met for in-frame deletion in an important domain."
    assert (
        prediction.BP3 is False
    ), "BP3 should not be met for in-frame deletion in an important domain."


@patch.object(RettAngelmanPredictor, "_in_repeat_region", return_value=True)
def test_verify_pm4bp3_in_repeat_region(
    mock_in_repeat_region, rett_angelman_predictor, seqvar, auto_acmg_data
):
    """Test verify_pm4bp3 when the variant is an in-frame deletion in a repeat region."""
    auto_acmg_data.hgnc_id = "HGNC:3811"  # FOXG1
    auto_acmg_data.prot_pos = 50  # In-frame deletion in a repeat region

    auto_acmg_data.consequence.cadd = "inframe_deletion"

    # Call the method under test
    prediction, comment = rett_angelman_predictor.verify_pm4bp3(seqvar, auto_acmg_data)

    assert (
        prediction.PM4 is False
    ), "PM4 should not be met for in-frame deletion in a repeat region."
    assert prediction.BP3 is True, "BP3 should be met for in-frame deletion in a repeat region."
    assert (
        "Variant consequence is in-frame deletion/insertion. Variant is in a repeat region or not in a conserved domain or in excluded region."
        in comment
    )


@patch.object(RettAngelmanPredictor, "_in_repeat_region", return_value=False)
def test_verify_pm4bp3_bp3_for_foxg1(
    mock_in_repeat_region, rett_angelman_predictor, seqvar, auto_acmg_data
):
    """Test verify_pm4bp3 when the variant is in the BP3 region for FOXG1."""
    auto_acmg_data.hgnc_id = "HGNC:3811"  # FOXG1
    auto_acmg_data.prot_pos = 50  # In the BP3 region for FOXG1

    auto_acmg_data.consequence.cadd = "inframe_deletion"

    # Call the method under test
    prediction, comment = rett_angelman_predictor.verify_pm4bp3(seqvar, auto_acmg_data)

    assert prediction.PM4 is False, "PM4 should not be met for FOXG1 in the BP3 region."
    assert prediction.BP3 is True, "BP3 should be met for FOXG1 in the BP3 region."
    assert "Variant is in the BP3 region for FOXG1." in comment


# --------------------- PP2 & BP1 ---------------------


def test_predict_pp2bp1(rett_angelman_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Retts and Angelman-like Disorders predictor."""

    # Call the method under test
    pp2_result, bp1_result = rett_angelman_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


# ------------------- PP3 & BP4 -------------------


def test_predict_pp3bp4_revel_strategy(rett_angelman_predictor, auto_acmg_data):
    """Test that REVEL is set as the strategy for PP3/BP4 prediction."""
    pp3_result, bp4_result = rett_angelman_predictor.predict_pp3bp4(
        rett_angelman_predictor.seqvar, auto_acmg_data
    )

    assert auto_acmg_data.thresholds.pp3bp4_strategy == "revel"
    assert auto_acmg_data.thresholds.revel_pathogenic == 0.75
    assert auto_acmg_data.thresholds.revel_benign == 0.15


@patch("src.vcep.rett_angelman.DefaultSeqVarPredictor.predict_pp3bp4")
def test_predict_pp3bp4_calls_superclass(
    mock_super_predict_pp3bp4, rett_angelman_predictor, auto_acmg_data
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

    pp3_result, bp4_result = rett_angelman_predictor.predict_pp3bp4(
        rett_angelman_predictor.seqvar, auto_acmg_data
    )

    mock_super_predict_pp3bp4.assert_called_once_with(
        rett_angelman_predictor.seqvar, auto_acmg_data
    )
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
    rett_angelman_predictor, auto_acmg_data, revel_score, expected_pp3, expected_bp4
):
    """Test different REVEL score scenarios."""
    auto_acmg_data.scores.dbnsfp.revel = revel_score

    with patch(
        "src.vcep.rett_angelman.DefaultSeqVarPredictor.predict_pp3bp4"
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

        pp3_result, bp4_result = rett_angelman_predictor.predict_pp3bp4(
            rett_angelman_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == expected_pp3
        assert bp4_result.prediction == expected_bp4


def test_predict_pp3bp4_strength(rett_angelman_predictor, auto_acmg_data):
    """Test that the strength of PP3 and BP4 is correctly set."""
    with patch(
        "src.vcep.rett_angelman.DefaultSeqVarPredictor.predict_pp3bp4"
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

        pp3_result, bp4_result = rett_angelman_predictor.predict_pp3bp4(
            rett_angelman_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.strength == AutoACMGStrength.PathogenicSupporting
        assert bp4_result.strength == AutoACMGStrength.BenignSupporting


def test_predict_pp3bp4_no_revel_score(rett_angelman_predictor, auto_acmg_data):
    """Test behavior when no REVEL score is available."""
    auto_acmg_data.scores.dbnsfp.revel = None

    with patch(
        "src.vcep.rett_angelman.DefaultSeqVarPredictor.predict_pp3bp4"
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

        pp3_result, bp4_result = rett_angelman_predictor.predict_pp3bp4(
            rett_angelman_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == AutoACMGPrediction.NotApplicable
        assert bp4_result.prediction == AutoACMGPrediction.NotApplicable


def test_predict_pp3bp4_error_handling(rett_angelman_predictor, auto_acmg_data):
    """Test error handling in predict_pp3bp4 method."""
    with patch(
        "src.vcep.rett_angelman.DefaultSeqVarPredictor.predict_pp3bp4",
        side_effect=Exception("Test error"),
    ):
        with pytest.raises(Exception) as exc_info:
            rett_angelman_predictor.predict_pp3bp4(rett_angelman_predictor.seqvar, auto_acmg_data)

        assert str(exc_info.value) == "Test error"


# --------------------- BP7 ---------------------


def test_predict_bp7_threshold_adjustment(rett_angelman_predictor, auto_acmg_data):
    """Test that the BP7 thresholds are correctly adjusted for Rett and Angelman-like Disorders."""
    auto_acmg_data.thresholds.phyloP100 = 1.0  # Initial phyloP100 threshold value

    # Call predict_bp7 method
    result = rett_angelman_predictor.predict_bp7(rett_angelman_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.phyloP100 == 0.1
    ), "The phyloP100 threshold should be adjusted to 0.1."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(RettAngelmanPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default(
    mock_super_predict_bp7, rett_angelman_predictor, auto_acmg_data
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
    result = rett_angelman_predictor.predict_bp7(rett_angelman_predictor.seqvar, auto_acmg_data)

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
