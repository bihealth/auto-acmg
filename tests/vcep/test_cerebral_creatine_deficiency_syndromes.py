from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import (
    PS1PM5,
    AlleleCondition,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
)
from src.defs.exceptions import MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor
from src.vcep import CerebralCreatineDeficiencySyndromesPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def cerebral_creatine_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return CerebralCreatineDeficiencySyndromesPredictor(
        seqvar=seqvar, result=result, config=MagicMock()
    )


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_overrides(
    mock_verify_ps1pm5, cerebral_creatine_predictor, seqvar, auto_acmg_data
):
    """Test that the overridden verify_ps1pm5 method works correctly."""
    mock_verify_ps1pm5.return_value = (
        PS1PM5(PS1=True, PM5=True),
        "Variant affects splicing. PS1 is not applicable.",
    )
    # Setup the data
    auto_acmg_data.hgnc_id = "HGNC:4175"
    auto_acmg_data.consequence = MagicMock(mehari=["missense_variant"])
    auto_acmg_data.thresholds.spliceAI_acceptor_gain = 0.5
    auto_acmg_data.thresholds.spliceAI_acceptor_loss = 0.5
    auto_acmg_data.thresholds.spliceAI_donor_gain = 0.5
    auto_acmg_data.thresholds.spliceAI_donor_loss = 0.5
    auto_acmg_data.scores.cadd.spliceAI_acceptor_gain = 0.6
    auto_acmg_data.scores.cadd.spliceAI_donor_gain = 0.6

    prediction, comment = cerebral_creatine_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    assert not prediction.PS1, "PS1 should be marked as not applicable due to splicing effect."
    assert not prediction.PM5, "PM5 should be marked as not applicable due to splicing effect."


def test_predict_pm1_not_applicable(cerebral_creatine_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for Cerebral Creatine Deficiency Syndromes."""
    result = cerebral_creatine_predictor.predict_pm1(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for Cerebral Creatine Deficiency Syndromes."
    assert (
        result.summary == "PM1 is not applicable for Cerebral Creatine Deficiency Syndromes."
    ), "The summary should indicate that PM1 is not applicable for Cerebral Creatine Deficiency Syndromes."


def test_predict_pm1_strength(cerebral_creatine_predictor, auto_acmg_data):
    """Test the strength level returned by the Cerebral Creatine Deficiency Syndromes predictor."""
    result = cerebral_creatine_predictor.predict_pm1(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate for Cerebral Creatine Deficiency Syndromes."


def test_predict_pm1_name(cerebral_creatine_predictor, auto_acmg_data):
    """Test the name of the criteria returned by the Cerebral Creatine Deficiency Syndromes predictor."""
    result = cerebral_creatine_predictor.predict_pm1(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    assert result.name == "PM1", "The name of the criteria should be 'PM1'."


@patch("src.vcep.cerebral_creatine_deficiency_syndromes.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, cerebral_creatine_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method (if implemented)."""
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = cerebral_creatine_predictor.predict_pm1(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    # In this specific case, the fallback should never happen since PM1 is always not applicable,
    # but this test ensures that if something changes, the fallback works correctly.
    assert result.prediction == AutoACMGPrediction.NotApplicable, "PM1 should remain NotApplicable."


@patch.object(
    CerebralCreatineDeficiencySyndromesPredictor,
    "_get_allele_cond",
    return_value=AlleleCondition.Recessive,
)
@patch.object(
    CerebralCreatineDeficiencySyndromesPredictor,
    "_get_control_af",
    return_value=MagicMock(bySex=MagicMock(overall=MagicMock(ac=10, nhomalt=6))),
)
@patch.object(CerebralCreatineDeficiencySyndromesPredictor, "_get_any_af", return_value=None)
def test_check_zyg_homozygous_positive(
    mock_get_any_af,
    mock_get_control_af,
    mock_get_allele_cond,
    cerebral_creatine_predictor,
    seqvar,
    auto_acmg_data,
):
    cerebral_creatine_predictor.comment_pm2ba1bs1bs2 = ""
    assert cerebral_creatine_predictor._check_zyg(seqvar, auto_acmg_data) == True
    assert (
        "The variant is in a recessive (homozygous) disorder."
        in cerebral_creatine_predictor.comment_pm2ba1bs1bs2
    )


@patch.object(
    CerebralCreatineDeficiencySyndromesPredictor,
    "_get_allele_cond",
    return_value=AlleleCondition.Recessive,
)
@patch.object(
    CerebralCreatineDeficiencySyndromesPredictor,
    "_get_control_af",
    return_value=MagicMock(bySex=MagicMock(overall=MagicMock(ac=10, nhomalt=4))),
)
@patch.object(CerebralCreatineDeficiencySyndromesPredictor, "_get_any_af", return_value=None)
def test_check_zyg_homozygous_negative(
    mock_get_any_af,
    mock_get_control_af,
    mock_get_allele_cond,
    cerebral_creatine_predictor,
    seqvar,
    auto_acmg_data,
):
    cerebral_creatine_predictor.comment_pm2ba1bs1bs2 = ""
    assert cerebral_creatine_predictor._check_zyg(seqvar, auto_acmg_data) == False


@patch.object(
    CerebralCreatineDeficiencySyndromesPredictor,
    "_get_allele_cond",
    return_value=AlleleCondition.Recessive,
)
@patch.object(CerebralCreatineDeficiencySyndromesPredictor, "_get_control_af", return_value=None)
@patch.object(CerebralCreatineDeficiencySyndromesPredictor, "_get_any_af", return_value=None)
def test_check_zyg_missing_data_raises_error(
    mock_get_any_af,
    mock_get_control_af,
    mock_get_allele_cond,
    cerebral_creatine_predictor,
    seqvar,
    auto_acmg_data,
):
    cerebral_creatine_predictor.comment_pm2ba1bs1bs2 = ""
    with pytest.raises(MissingDataError):
        cerebral_creatine_predictor._check_zyg(seqvar, auto_acmg_data)


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
def test_predict_pm2ba1bs1bs2_gatm(
    mock_super_method, cerebral_creatine_predictor, auto_acmg_data, seqvar
):
    # Default thresholds
    auto_acmg_data.hgnc_id = "HGNC:4175"
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = cerebral_creatine_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.000055
    assert auto_acmg_data.thresholds.ba1_benign == 0.0005
    assert auto_acmg_data.thresholds.bs1_benign == 0.0001

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


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
def test_predict_pm2ba1bs1bs2_gamt(
    mock_super_method, cerebral_creatine_predictor, auto_acmg_data, seqvar
):
    # Default thresholds
    auto_acmg_data.hgnc_id = "HGNC:4136"  # GAMT
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = cerebral_creatine_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.0004
    assert auto_acmg_data.thresholds.ba1_benign == 0.003
    assert auto_acmg_data.thresholds.bs1_benign == 0.001

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


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
def test_predict_pm2ba1bs1bs2_slc6a8(
    mock_super_method, cerebral_creatine_predictor, auto_acmg_data, seqvar
):
    # Default thresholds
    auto_acmg_data.hgnc_id = "HGNC:11055"  # SLC6A8
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = cerebral_creatine_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.00002
    assert auto_acmg_data.thresholds.ba1_benign == 0.002
    assert auto_acmg_data.thresholds.bs1_benign == 0.0002

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_bp3_not_applicable(cerebral_creatine_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = cerebral_creatine_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


def test_predict_pp2bp1(cerebral_creatine_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Cerebral Creatine Deficiency Syndromes."""

    # Call the method under test
    pp2_result, bp1_result = cerebral_creatine_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

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


@pytest.mark.parametrize(
    "revel_score, expected_pp3, expected_bp4",
    [
        # (0.8, AutoACMGPrediction.Met, AutoACMGPrediction.NotMet),
        # (0.5, AutoACMGPrediction.NotMet, AutoACMGPrediction.NotMet),
        (0.1, AutoACMGPrediction.NotMet, AutoACMGPrediction.Met),
    ],
)
def test_predict_pp3bp4_revel_scenarios(
    cerebral_creatine_predictor, auto_acmg_data, revel_score, expected_pp3, expected_bp4
):
    auto_acmg_data.scores.dbnsfp.revel = revel_score

    pp3, bp4 = cerebral_creatine_predictor.predict_pp3bp4(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    assert pp3.prediction == expected_pp3
    assert bp4.prediction == expected_bp4
    if revel_score >= 0.75:
        assert f"REVEL score {revel_score} >= 0.75, meeting PP3." in pp3.summary
    elif revel_score <= 0.15:
        assert f"REVEL score {revel_score} <= 0.15, meeting BP4." in bp4.summary


@patch.object(CerebralCreatineDeficiencySyndromesPredictor, "_is_inframe_indel")
@patch.object(CerebralCreatineDeficiencySyndromesPredictor, "_affect_spliceAI")
def test_predict_pp3bp4_inframe_indel(
    mock_affect_spliceAI, mock_is_inframe_indel, cerebral_creatine_predictor, auto_acmg_data
):
    mock_is_inframe_indel.return_value = True
    mock_affect_spliceAI.return_value = False
    auto_acmg_data.scores.dbnsfp.provean = -3.0
    auto_acmg_data.scores.dbnsfp.mutationTaster = 0.6

    pp3, bp4 = cerebral_creatine_predictor.predict_pp3bp4(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    assert pp3.prediction == AutoACMGPrediction.Met
    assert "In-frame indel predicted deleterious by PROVEAN and MutationTaster." in pp3.summary


@patch.object(CerebralCreatineDeficiencySyndromesPredictor, "_affect_spliceAI")
def test_predict_pp3bp4_splicing(mock_affect_spliceAI, cerebral_creatine_predictor, auto_acmg_data):
    mock_affect_spliceAI.return_value = True

    pp3, bp4 = cerebral_creatine_predictor.predict_pp3bp4(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    assert pp3.prediction == AutoACMGPrediction.Met
    assert "Splicing predictions indicate an impact, meeting PP3." in pp3.summary


def test_predict_pp3bp4_no_criteria_met(cerebral_creatine_predictor, auto_acmg_data):
    auto_acmg_data.scores.dbnsfp.revel = 0.5
    with (
        patch.object(
            CerebralCreatineDeficiencySyndromesPredictor, "_is_inframe_indel", return_value=False
        ),
        patch.object(
            CerebralCreatineDeficiencySyndromesPredictor, "_affect_spliceAI", return_value=False
        ),
    ):
        pp3, bp4 = cerebral_creatine_predictor.predict_pp3bp4(
            cerebral_creatine_predictor.seqvar, auto_acmg_data
        )

    assert pp3.prediction == AutoACMGPrediction.NotMet
    assert bp4.prediction == AutoACMGPrediction.Met
    assert "No significant splicing impact predicted, meeting BP4." in bp4.summary


def test_predict_pp3bp4_strength(cerebral_creatine_predictor, auto_acmg_data):
    pp3, bp4 = cerebral_creatine_predictor.predict_pp3bp4(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    assert pp3.strength == AutoACMGStrength.PathogenicSupporting
    assert bp4.strength == AutoACMGStrength.BenignSupporting


@patch.object(CerebralCreatineDeficiencySyndromesPredictor, "_is_inframe_indel")
@patch.object(CerebralCreatineDeficiencySyndromesPredictor, "_affect_spliceAI")
def test_predict_pp3bp4_method_calls(
    mock_affect_spliceAI, mock_is_inframe_indel, cerebral_creatine_predictor, auto_acmg_data
):
    mock_is_inframe_indel.return_value = False
    mock_affect_spliceAI.return_value = True

    cerebral_creatine_predictor.predict_pp3bp4(cerebral_creatine_predictor.seqvar, auto_acmg_data)

    mock_is_inframe_indel.assert_called_once()
    mock_affect_spliceAI.assert_called_once()


def test_predict_pp3bp4_multiple_criteria(cerebral_creatine_predictor, auto_acmg_data):
    auto_acmg_data.scores.dbnsfp.revel = 0.8
    with (
        patch.object(
            CerebralCreatineDeficiencySyndromesPredictor, "_is_inframe_indel", return_value=True
        ),
        patch.object(
            CerebralCreatineDeficiencySyndromesPredictor, "_affect_spliceAI", return_value=True
        ),
    ):
        pp3, bp4 = cerebral_creatine_predictor.predict_pp3bp4(
            cerebral_creatine_predictor.seqvar, auto_acmg_data
        )

    assert pp3.prediction == AutoACMGPrediction.Met
    assert "REVEL score 0.8 >= 0.75, meeting PP3." in pp3.summary
    assert "Splicing predictions indicate an impact, meeting PP3." in pp3.summary
    assert bp4.prediction == AutoACMGPrediction.NotMet


def test_predict_bp7_threshold_adjustment(cerebral_creatine_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted."""
    auto_acmg_data.thresholds.bp7_donor = 1  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value
    # Set the hgnc_id to GAMT (HGNC:4136) to trigger the threshold adjustment
    auto_acmg_data.hgnc_id = "HGNC:4136"

    # Call predict_bp7 method
    result = cerebral_creatine_predictor.predict_bp7(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.bp7_donor == 7
    ), "The BP7 donor threshold should be adjusted to 7."
    assert (
        auto_acmg_data.thresholds.bp7_acceptor == 21
    ), "The BP7 acceptor threshold should be adjusted to 4."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(DefaultSeqVarPredictor, "predict_bp7")
def test_predict_bp7_fallback_to_default(
    mock_super_predict_bp7, cerebral_creatine_predictor, auto_acmg_data
):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotMet,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = cerebral_creatine_predictor.predict_bp7(
        cerebral_creatine_predictor.seqvar, auto_acmg_data
    )

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
