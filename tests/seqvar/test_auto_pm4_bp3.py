from unittest.mock import MagicMock, patch

import pytest
import tabix

from src.defs.auto_acmg import AutoACMGPrediction, AutoACMGSeqVarData, AutoACMGStrength
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.auto_pm4_bp3 import AutoPM4BP3


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def auto_pm4bp3():
    return AutoPM4BP3()


@pytest.fixture
def var_data():
    data = MagicMock(spec=AutoACMGSeqVarData)
    data.consequence = MagicMock()
    return data


# ============== _in_repeat_region =================


@patch("src.seqvar.auto_pm4_bp3.tabix.open")
def test_in_repeat_region_true(mock_tabix_open, auto_pm4bp3, seqvar):
    mock_tabix = MagicMock()
    mock_tabix.query.return_value = iter(
        [["chr1", "99", "101"]]
    )  # Simulating a record that matches the position
    mock_tabix_open.return_value = mock_tabix

    result = auto_pm4bp3._in_repeat_region(seqvar)
    assert result is True


@patch("src.seqvar.auto_pm4_bp3.tabix.open")
def test_in_repeat_region_false(mock_tabix_open, auto_pm4bp3, seqvar):
    mock_tabix = MagicMock()
    mock_tabix.query.return_value = iter([])  # No records returned
    mock_tabix_open.return_value = mock_tabix

    result = auto_pm4bp3._in_repeat_region(seqvar)
    assert result is False


@patch("src.seqvar.auto_pm4_bp3.tabix.open")
def test_in_repeat_region_error(mock_tabix_open, auto_pm4bp3, seqvar):
    mock_tabix_open.side_effect = tabix.TabixError("Failed to open file")

    with pytest.raises(AlgorithmError) as excinfo:
        auto_pm4bp3._in_repeat_region(seqvar)
    assert "Failed to check if the variant is in a repeat region." in str(excinfo.value)


# ============== _is_stop_loss =================


@pytest.fixture
def var_data_is_stop_loss():
    data = MagicMock()
    data.consequence = MagicMock()
    return data


def test_is_stop_loss_cadd_true(var_data_is_stop_loss):
    var_data_is_stop_loss.consequence.cadd = "some_effect stop_loss other_effect"
    assert AutoPM4BP3()._is_stop_loss(var_data_is_stop_loss) is True


def test_is_stop_loss_mehari_true(var_data_is_stop_loss):
    var_data_is_stop_loss.consequence.cadd = "some_effect other_effect"
    var_data_is_stop_loss.consequence.mehari = [
        "some_effect",
        "stop_loss",
        "other_effect",
    ]
    assert AutoPM4BP3()._is_stop_loss(var_data_is_stop_loss) is True


def test_is_stop_loss_false(var_data_is_stop_loss):
    var_data_is_stop_loss.consequence.cadd = "some_effect other_effect"
    var_data_is_stop_loss.consequence.mehari = ["some_effect", "other_effect"]
    assert AutoPM4BP3()._is_stop_loss(var_data_is_stop_loss) is False


# ============== _is_inframe_delins =================


@pytest.fixture
def var_data_is_inframe_delins():
    data = MagicMock()
    data.consequence = MagicMock()
    return data


def test_inframe_deletion_cadd_true(var_data_is_inframe_delins):
    var_data_is_inframe_delins.consequence.cadd = "some_effect inframe_deletion other_effect"
    assert AutoPM4BP3().is_inframe_delins(var_data_is_inframe_delins) is True


def test_inframe_insertion_cadd_true(var_data_is_inframe_delins):
    var_data_is_inframe_delins.consequence.cadd = "inframe_insertion some_effect other_effect"
    assert AutoPM4BP3().is_inframe_delins(var_data_is_inframe_delins) is True


def test_inframe_deletion_mehari_true(var_data_is_inframe_delins):
    var_data_is_inframe_delins.consequence.cadd = "some_effect other_effect"
    var_data_is_inframe_delins.consequence.mehari = [
        "some_effect",
        "inframe_deletion",
        "other_effect",
    ]
    assert AutoPM4BP3().is_inframe_delins(var_data_is_inframe_delins) is True


def test_inframe_insertion_mehari_true(var_data_is_inframe_delins):
    var_data_is_inframe_delins.consequence.cadd = "some_effect other_effect"
    var_data_is_inframe_delins.consequence.mehari = [
        "some_effect",
        "inframe_insertion",
        "other_effect",
    ]
    assert AutoPM4BP3().is_inframe_delins(var_data_is_inframe_delins) is True


def test_inframe_delins_false(var_data_is_inframe_delins):
    var_data_is_inframe_delins.consequence.cadd = "some_effect other_effect"
    var_data_is_inframe_delins.consequence.mehari = ["some_effect", "other_effect"]
    assert AutoPM4BP3().is_inframe_delins(var_data_is_inframe_delins) is False


# ================== _bp3_not_applicable ==================


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


def test_bp3_not_applicable_mitochondrial(auto_pm4bp3, seqvar, auto_acmg_data):
    """Test BP3 is not applicable when the variant is in mitochondrial DNA."""
    seqvar.chrom = "MT"  # Mitochondrial chromosome
    result = auto_pm4bp3._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should not be applicable for mitochondrial DNA variants."


def test_bp3_applicable_nuclear_dna(auto_pm4bp3, seqvar, auto_acmg_data):
    """Test BP3 is applicable when the variant is in nuclear DNA."""
    seqvar.chrom = "1"  # Nuclear chromosome
    result = auto_pm4bp3._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is False, "BP3 should be applicable for nuclear DNA variants."


# ============== verify_pm4bp3 =================


@patch.object(AutoPM4BP3, "_is_stop_loss", return_value=True)
def test_stop_loss_variant(mock_is_stop_loss, auto_pm4bp3, seqvar, var_data):
    prediction, comment = auto_pm4bp3.verify_pm4bp3(seqvar, var_data)
    assert prediction.PM4 is True
    assert prediction.BP3 is False
    assert "stop-loss" in comment


@patch.object(AutoPM4BP3, "_is_stop_loss", return_value=False)
@patch.object(AutoPM4BP3, "is_inframe_delins", return_value=True)
@patch.object(AutoPM4BP3, "_in_repeat_region", return_value=False)
def test_inframe_delins_outside_repeat(
    mock_in_repeat_region,
    mock_is_inframe_delins,
    mock_is_stop_loss,
    auto_pm4bp3,
    seqvar,
    var_data,
):
    prediction, comment = auto_pm4bp3.verify_pm4bp3(seqvar, var_data)
    assert prediction.PM4 is True
    assert prediction.BP3 is False
    assert "not in a repeat region" in comment


@patch.object(AutoPM4BP3, "_is_stop_loss", return_value=False)
@patch.object(AutoPM4BP3, "is_inframe_delins", return_value=True)
@patch.object(AutoPM4BP3, "_in_repeat_region", return_value=True)
def test_inframe_delins_inside_repeat(
    mock_in_repeat_region,
    mock_is_inframe_delins,
    mock_is_stop_loss,
    auto_pm4bp3,
    seqvar,
    var_data,
):
    prediction, comment = auto_pm4bp3.verify_pm4bp3(seqvar, var_data)
    assert prediction.PM4 is False
    assert prediction.BP3 is True
    assert "in a repeat region" in comment


@patch.object(AutoPM4BP3, "_is_stop_loss", return_value=False)
@patch.object(AutoPM4BP3, "is_inframe_delins", return_value=False)
def test_variant_not_meeting_criteria(
    mock_is_inframe_delins, mock_is_stop_loss, auto_pm4bp3, seqvar, var_data
):
    prediction, comment = auto_pm4bp3.verify_pm4bp3(seqvar, var_data)
    assert prediction.PM4 is False
    assert prediction.BP3 is False
    assert "not met" in comment


@patch.object(AutoPM4BP3, "_is_stop_loss", side_effect=AutoAcmgBaseException("Test Error"))
def test_exception_handling(mock_is_stop_loss, auto_pm4bp3, seqvar, var_data):
    prediction, comment = auto_pm4bp3.verify_pm4bp3(seqvar, var_data)
    assert prediction is None
    assert "Test Error" in comment


# ============== predict_pm4bp3 ==============


@patch.object(AutoPM4BP3, "verify_pm4bp3")
def test_predict_pm4bp3_success(mock_verify, auto_pm4bp3, seqvar, var_data):
    # Setup mock to return a successful prediction for both PM4 and BP3
    mock_verify.return_value = (
        MagicMock(
            PM4=True,
            BP3=False,
            PM4_strength=AutoACMGStrength.PathogenicModerate,
            BP3_strength=AutoACMGStrength.BenignSupporting,
        ),
        "Verification successful.",
    )
    pm4_criteria, bp3_criteria = auto_pm4bp3.predict_pm4bp3(seqvar, var_data)

    assert pm4_criteria.name == "PM4"
    assert pm4_criteria.prediction == AutoACMGPrediction.Applicable
    assert pm4_criteria.strength == AutoACMGStrength.PathogenicModerate
    assert "Verification successful." in pm4_criteria.summary

    assert bp3_criteria.name == "BP3"
    assert bp3_criteria.prediction == AutoACMGPrediction.NotApplicable
    assert bp3_criteria.strength == AutoACMGStrength.BenignSupporting
    assert "Verification successful." in bp3_criteria.summary


@patch.object(AutoPM4BP3, "verify_pm4bp3")
def test_predict_pm4bp3_failure(mock_verify, auto_pm4bp3, seqvar, var_data):
    # Setup mock to return None indicating a failed verification process
    mock_verify.return_value = (None, "Verification failed due to an error.")
    pm4_criteria, bp3_criteria = auto_pm4bp3.predict_pm4bp3(seqvar, var_data)

    assert pm4_criteria.name == "PM4"
    assert pm4_criteria.prediction == AutoACMGPrediction.Failed
    assert pm4_criteria.strength == AutoACMGStrength.PathogenicModerate
    assert "Verification failed due to an error." in pm4_criteria.summary

    assert bp3_criteria.name == "BP3"
    assert bp3_criteria.prediction == AutoACMGPrediction.Failed
    assert bp3_criteria.strength == AutoACMGStrength.BenignSupporting
    assert "Verification failed due to an error." in bp3_criteria.summary
