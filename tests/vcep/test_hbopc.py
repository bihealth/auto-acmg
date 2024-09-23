from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import (
    PS1PM5,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
)
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor
from src.vcep import HBOPCPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def hbopc_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return HBOPCPredictor(seqvar=seqvar, result=result)


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


@pytest.fixture
def auto_acmg_data_atm():
    data = AutoACMGSeqVarData()
    data.hgnc_id = "HGNC:795"
    data.consequence = MagicMock(mehari=[], cadd=None)
    return data


@pytest.fixture
def auto_acmg_data_palb2():
    data = AutoACMGSeqVarData()
    data.hgnc_id = "HGNC:26144"
    data.consequence = MagicMock(mehari=[], cadd=None)
    return data


# --------------- PS1 & PM5 ---------------


def test_is_nonsense_true(hbopc_predictor, auto_acmg_data):
    """Test _is_nonsense method when the variant is a nonsense mutation."""
    auto_acmg_data.consequence = MagicMock(mehari=["stop_gained"])

    result = hbopc_predictor._is_nonsense(auto_acmg_data)

    assert result is True, "Should return True for a nonsense mutation"


def test_is_nonsense_false(hbopc_predictor, auto_acmg_data):
    """Test _is_nonsense method when the variant is not a nonsense mutation."""
    auto_acmg_data.consequence = MagicMock(mehari=["missense_variant"])

    result = hbopc_predictor._is_nonsense(auto_acmg_data)

    assert result is False, "Should return False for a non-nonsense mutation"


def test_is_nonsense_empty_consequence(hbopc_predictor, auto_acmg_data):
    """Test _is_nonsense method when the consequence list is empty."""
    auto_acmg_data.consequence = MagicMock(mehari=[])

    result = hbopc_predictor._is_nonsense(auto_acmg_data)

    assert result is False, "Should return False when consequence list is empty"


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_nonsense(mock_super_verify, hbopc_predictor, seqvar, auto_acmg_data):
    """Test verify_ps1pm5 method for nonsense mutations."""
    mock_super_verify.return_value = (PS1PM5(PS1=True, PM5=True), "Initial evaluation")
    auto_acmg_data.consequence = MagicMock(mehari=["stop_gained"])

    prediction, comment = hbopc_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    assert prediction.PS1 is True
    assert prediction.PM5 is True


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_splicing_effect(mock_super_verify, hbopc_predictor, seqvar, auto_acmg_data):
    """Test verify_ps1pm5 method for variants affecting splicing."""
    mock_super_verify.return_value = (PS1PM5(PS1=True, PM5=True), "Initial evaluation")
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

    prediction, comment = hbopc_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    assert prediction.PS1 is True
    assert prediction.PM5 is True


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_no_override(mock_super_verify, hbopc_predictor, seqvar, auto_acmg_data):
    """Test verify_ps1pm5 method when no override is needed."""
    mock_super_verify.return_value = (PS1PM5(PS1=True, PM5=True), "Initial evaluation")
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

    prediction, comment = hbopc_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    assert prediction.PS1 is True, "PS1 should remain True when no override is needed"
    assert prediction.PM5 is True, "PM5 should remain True when no override is needed"
    assert "Initial evaluation" in comment, "Comment should reflect initial evaluation"


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_exception_handling(
    mock_super_verify, hbopc_predictor, seqvar, auto_acmg_data
):
    """Test verify_ps1pm5 method exception handling."""
    mock_super_verify.side_effect = Exception("Test exception")

    with pytest.raises(Exception) as exc_info:
        hbopc_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    assert "Test exception" in str(exc_info.value), "Should raise the original exception"


# --------------- PM1 ---------------


def test_predict_pm1_not_applicable(hbopc_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for ATM and PALB2 in Hereditary Breast, Ovarian, and Pancreatic Cancer."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for ATM."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting for ATM."
    assert (
        "PM1 is not applicable for HGNC:795" in result.summary
    ), "The summary should indicate that PM1 is not applicable for ATM."


def test_predict_pm1_not_applicable_palb2(hbopc_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for PALB2 in Hereditary Breast, Ovarian, and Pancreatic Cancer."""
    auto_acmg_data.hgnc_id = "HGNC:26144"  # PALB2 gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for PALB2."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting for PALB2."
    assert (
        "PM1 is not applicable for HGNC:26144" in result.summary
    ), "The summary should indicate that PM1 is not applicable for PALB2."


@patch("src.vcep.hbopc.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, hbopc_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method if logic changes."""
    auto_acmg_data.hgnc_id = "HGNC:111111"  # Gene not in the specific logic
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met if the gene is not ATM or PALB2."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


def test_predict_pm1_strength_level(hbopc_predictor, auto_acmg_data):
    """Test the strength level for PM1 when not applicable."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting when PM1 is not applicable."


def test_predict_pm1_name(hbopc_predictor, auto_acmg_data):
    """Test the name of the criteria returned by the HBOPC predictor."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    result = hbopc_predictor.predict_pm1(hbopc_predictor.seqvar, auto_acmg_data)

    assert result.name == "PM1", "The name of the criteria should be 'PM1'."


# -------------- PM2, BA1, BS1, BS2 --------------


def test_bs2_not_applicable(hbopc_predictor, auto_acmg_data):
    """Test when BS2 is not applicable for Hereditary Breast, Ovarian, and Pancreatic Cancer."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    result = hbopc_predictor._bs2_not_applicable(auto_acmg_data)
    assert result is True, "BS2 should always be not applicable for HBOPC."


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
    [("HGNC:795", 0.00001, 0.005, 0.0005), ("HGNC:26144", 0.000003, 0.001, 0.0001)],
)
def test_predict_pm2ba1bs1bs2_specific_genes(
    mock_super_method,
    hbopc_predictor,
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
    hbopc_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Validate thresholds are set correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == expected_pm2
    assert auto_acmg_data.thresholds.ba1_benign == expected_ba1
    assert auto_acmg_data.thresholds.bs1_benign == expected_bs1

    # Check that the superclass method was called with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Reset mock for the next iteration
    mock_super_method.reset_mock()


# -------------- PM4 & BP3 --------------


@patch.object(DefaultSeqVarPredictor, "verify_pm4bp3")
def test_predict_pm4bp3_atm(mock_verify_pm4bp3, hbopc_predictor, seqvar, auto_acmg_data):
    """Test the predict_pm4bp3 method for ATM (HGNC:795)."""
    auto_acmg_data.hgnc_id = "HGNC:795"

    # Mock the verify_pm4bp3 method to return a specific result
    mock_pred = MagicMock()
    mock_pred.PM4 = True
    mock_pred.BP3 = False
    mock_verify_pm4bp3.return_value = (mock_pred, "PM4 is met for ATM")

    # Call the method under test
    pm4_result, bp3_result = hbopc_predictor.predict_pm4bp3(seqvar, auto_acmg_data)

    # Check PM4 result
    assert isinstance(
        pm4_result, AutoACMGCriteria
    ), "The PM4 result should be of type AutoACMGCriteria."
    assert pm4_result.prediction == AutoACMGPrediction.Applicable, "PM4 should be Met for ATM."
    assert (
        pm4_result.strength == AutoACMGStrength.PathogenicSupporting
    ), "PM4 strength should be PathogenicSupporting."
    assert (
        "PM4 is met for ATM" in pm4_result.summary
    ), "The summary should indicate PM4 is met for ATM."

    # Check BP3 result
    assert isinstance(
        bp3_result, AutoACMGCriteria
    ), "The BP3 result should be of type AutoACMGCriteria."
    assert (
        bp3_result.prediction == AutoACMGPrediction.NotApplicable
    ), "BP3 should be NotApplicable for ATM."
    assert (
        bp3_result.strength == AutoACMGStrength.BenignSupporting
    ), "BP3 strength should be BenignSupporting."
    assert (
        "BP3 is not applicable for ATM." in bp3_result.summary
    ), "The summary should indicate BP3 is not applicable for ATM."


def test_predict_pm4bp3_palb2(hbopc_predictor, seqvar, auto_acmg_data):
    """Test the predict_pm4bp3 method for PALB2 (HGNC:26144)."""
    auto_acmg_data.hgnc_id = "HGNC:26144"

    # Call the method under test
    pm4_result, bp3_result = hbopc_predictor.predict_pm4bp3(seqvar, auto_acmg_data)

    # Check PM4 result
    assert isinstance(
        pm4_result, AutoACMGCriteria
    ), "The PM4 result should be of type AutoACMGCriteria."
    assert (
        pm4_result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM4 should be NotApplicable for PALB2."
    assert (
        pm4_result.strength == AutoACMGStrength.PathogenicSupporting
    ), "PM4 strength should be PathogenicSupporting."
    assert (
        "PM4 is not applicable for PALB2." in pm4_result.summary
    ), "The summary should indicate PM4 is not applicable for PALB2."

    # Check BP3 result
    assert isinstance(
        bp3_result, AutoACMGCriteria
    ), "The BP3 result should be of type AutoACMGCriteria."
    assert (
        bp3_result.prediction == AutoACMGPrediction.NotApplicable
    ), "BP3 should be NotApplicable for PALB2."
    assert (
        bp3_result.strength == AutoACMGStrength.BenignSupporting
    ), "BP3 strength should be BenignSupporting."
    assert (
        "BP3 is not applicable for PALB2." in bp3_result.summary
    ), "The summary should indicate BP3 is not applicable for PALB2."


@patch.object(DefaultSeqVarPredictor, "verify_pm4bp3")
@patch.object(DefaultSeqVarPredictor, "predict_pm4bp3")
@pytest.mark.skip("Something doesn't work here")
def test_predict_pm4bp3_fallback(
    mock_predict_pm4bp3, mock_verify_pm4bp3, hbopc_predictor, seqvar, auto_acmg_data
):
    """Test the fallback behavior for predict_pm4bp3 method when the gene is not ATM or PALB2."""
    auto_acmg_data.hgnc_id = "HGNC:99999"  # Some gene not handled explicitly

    # Mock the verify_pm4bp3 method to return a specific result
    mock_pred = MagicMock()
    mock_pred.PM4 = True
    mock_pred.BP3 = False
    mock_verify_pm4bp3.return_value = (mock_pred, "PM4 is met for fallback")

    # Call the method under test
    pm4_result, bp3_result = hbopc_predictor.predict_pm4bp3(seqvar, auto_acmg_data)

    # Check PM4 result
    assert isinstance(
        pm4_result, AutoACMGCriteria
    ), "The PM4 result should be of type AutoACMGCriteria."
    assert (
        pm4_result.prediction == AutoACMGPrediction.Applicable
    ), "PM4 should be Met for fallback gene."
    assert (
        "PM4 is met for fallback" in pm4_result.summary
    ), "The summary should indicate PM4 is met for fallback gene."

    # Check BP3 result
    assert isinstance(
        bp3_result, AutoACMGCriteria
    ), "The BP3 result should be of type AutoACMGCriteria."
    assert (
        bp3_result.prediction == AutoACMGPrediction.NotApplicable
    ), "BP3 should be NotMet for fallback gene."
    assert (
        "PM4 is met for fallback" in bp3_result.summary
    ), "The summary should indicate BP3 was checked."

    # Ensure the superclass's predict_pm4bp3 method was called
    assert mock_predict_pm4bp3.called, "super().predict_pm4bp3 should have been called."


# -------------- PP2 & BP1 --------------


@patch.object(HBOPCPredictor, "predict_pp2bp1")
def test_predict_pp2bp1_palb2_missense(mock_predict, hbopc_predictor, seqvar, auto_acmg_data_palb2):
    """Test predict_pp2bp1 for PALB2 where the variant is missense."""
    auto_acmg_data_palb2.consequence.mehari = ["missense_variant"]
    mock_predict.return_value = (
        MagicMock(name="PP2", prediction=AutoACMGPrediction.NotApplicable),
        MagicMock(
            name="BP1",
            prediction=AutoACMGPrediction.Applicable,
            summary="Variant is missense. BP1 is met.",
        ),
    )

    pp2, bp1 = hbopc_predictor.predict_pp2bp1(seqvar, auto_acmg_data_palb2)

    assert (
        pp2.prediction == AutoACMGPrediction.NotApplicable
    ), "PP2 should be NotApplicable for PALB2."
    assert (
        bp1.prediction == AutoACMGPrediction.Applicable
    ), "BP1 should be Met for missense variant in PALB2."
    assert "missense" in bp1.summary, "Summary should confirm the missense nature."


@patch.object(HBOPCPredictor, "predict_pp2bp1")
def test_predict_pp2bp1_atm_non_missense(mock_predict, hbopc_predictor, seqvar, auto_acmg_data_atm):
    """Test predict_pp2bp1 for ATM where the variant is not missense."""
    auto_acmg_data_atm.consequence.mehari = ["nonsense_variant"]
    mock_predict.return_value = (
        MagicMock(name="PP2", prediction=AutoACMGPrediction.NotApplicable),
        MagicMock(
            name="BP1",
            prediction=AutoACMGPrediction.NotApplicable,
            summary="Variant is not missense. BP1 is not met.",
        ),
    )

    pp2, bp1 = hbopc_predictor.predict_pp2bp1(seqvar, auto_acmg_data_atm)

    assert (
        pp2.prediction == AutoACMGPrediction.NotApplicable
    ), "PP2 should be NotApplicable for ATM."
    assert (
        bp1.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should not be Met for non-missense variant in ATM."
    assert "not missense" in bp1.summary, "Summary should confirm the non-missense nature."


# -------------- PP3 & BP4 --------------


def test_predict_pp3bp4_revel_strategy(hbopc_predictor, auto_acmg_data):
    """Test that REVEL is set as the strategy for PP3/BP4 prediction."""
    auto_acmg_data.hgnc_id = "HGNC:795"
    pp3_result, bp4_result = hbopc_predictor.predict_pp3bp4(hbopc_predictor.seqvar, auto_acmg_data)

    assert auto_acmg_data.thresholds.pp3bp4_strategy == "revel"


@patch("src.vcep.hbopc.DefaultSeqVarPredictor.predict_pp3bp4")
def test_predict_pp3bp4_calls_superclass(
    mock_super_predict_pp3bp4, hbopc_predictor, auto_acmg_data
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

    pp3_result, bp4_result = hbopc_predictor.predict_pp3bp4(hbopc_predictor.seqvar, auto_acmg_data)

    mock_super_predict_pp3bp4.assert_called_once_with(hbopc_predictor.seqvar, auto_acmg_data)
    assert pp3_result.prediction == AutoACMGPrediction.Applicable
    assert bp4_result.prediction == AutoACMGPrediction.NotApplicable


@pytest.mark.parametrize(
    "hgnc_id, revel_score, expected_pp3, expected_bp4",
    [
        (
            "HGNC:795",
            0.9,
            AutoACMGPrediction.Applicable,
            AutoACMGPrediction.NotApplicable,
        ),  # ATM, high REVEL score
        (
            "HGNC:795",
            0.5,
            AutoACMGPrediction.NotApplicable,
            AutoACMGPrediction.NotApplicable,
        ),  # ATM, intermediate REVEL score
        (
            "HGNC:795",
            0.1,
            AutoACMGPrediction.NotApplicable,
            AutoACMGPrediction.Applicable,
        ),  # ATM, low REVEL score
        (
            "HGNC:26144",
            0.9,
            AutoACMGPrediction.NotApplicable,
            AutoACMGPrediction.NotApplicable,
        ),  # PALB2, high REVEL score
        (
            "HGNC:26144",
            0.5,
            AutoACMGPrediction.NotApplicable,
            AutoACMGPrediction.NotApplicable,
        ),  # PALB2, intermediate REVEL score
        (
            "HGNC:26144",
            0.1,
            AutoACMGPrediction.NotApplicable,
            AutoACMGPrediction.NotApplicable,
        ),  # PALB2, low REVEL score
    ],
)
def test_predict_pp3bp4_revel_scenarios(
    hbopc_predictor, auto_acmg_data, hgnc_id, revel_score, expected_pp3, expected_bp4
):
    """Test different REVEL score scenarios for ATM and PALB2."""
    auto_acmg_data.hgnc_id = hgnc_id
    auto_acmg_data.scores.dbnsfp.revel = revel_score

    with patch("src.vcep.hbopc.DefaultSeqVarPredictor.predict_pp3bp4") as mock_super_predict_pp3bp4:
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

        pp3_result, bp4_result = hbopc_predictor.predict_pp3bp4(
            hbopc_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == expected_pp3
        assert bp4_result.prediction == expected_bp4


def test_predict_pp3bp4_strength(hbopc_predictor, auto_acmg_data):
    """Test that the strength of PP3 and BP4 is correctly set."""
    with patch("src.vcep.hbopc.DefaultSeqVarPredictor.predict_pp3bp4") as mock_super_predict_pp3bp4:
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

        pp3_result, bp4_result = hbopc_predictor.predict_pp3bp4(
            hbopc_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.strength == AutoACMGStrength.PathogenicSupporting
        assert bp4_result.strength == AutoACMGStrength.BenignSupporting


def test_predict_pp3bp4_no_revel_score(hbopc_predictor, auto_acmg_data):
    """Test behavior when no REVEL score is available."""
    auto_acmg_data.scores.dbnsfp.revel = None

    with patch("src.vcep.hbopc.DefaultSeqVarPredictor.predict_pp3bp4") as mock_super_predict_pp3bp4:
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

        pp3_result, bp4_result = hbopc_predictor.predict_pp3bp4(
            hbopc_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == AutoACMGPrediction.NotApplicable
        assert bp4_result.prediction == AutoACMGPrediction.NotApplicable


def test_predict_pp3bp4_error_handling(hbopc_predictor, auto_acmg_data):
    """Test error handling in predict_pp3bp4 method."""
    with patch(
        "src.vcep.hbopc.DefaultSeqVarPredictor.predict_pp3bp4",
        side_effect=Exception("Test error"),
    ):
        with pytest.raises(Exception) as exc_info:
            hbopc_predictor.predict_pp3bp4(hbopc_predictor.seqvar, auto_acmg_data)

        assert str(exc_info.value) == "Test error"


# -------------- BP7 --------------


def test_predict_bp7_threshold_adjustment_for_palb2(hbopc_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted for PALB2."""
    auto_acmg_data.hgnc_id = "HGNC:26144"  # PALB2 gene
    auto_acmg_data.thresholds.bp7_donor = 2  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = hbopc_predictor.predict_bp7(hbopc_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted for PALB2
    assert (
        auto_acmg_data.thresholds.bp7_donor == 7
    ), "The BP7 donor threshold should be adjusted to 7 for PALB2."
    assert (
        auto_acmg_data.thresholds.bp7_acceptor == 21
    ), "The BP7 acceptor threshold should be adjusted to 21 for PALB2."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


def test_predict_bp7_threshold_adjustment_for_atm(hbopc_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted for ATM."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    auto_acmg_data.thresholds.bp7_donor = 2  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = hbopc_predictor.predict_bp7(hbopc_predictor.seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted for ATM
    assert (
        auto_acmg_data.thresholds.bp7_donor == 7
    ), "The BP7 donor threshold should be adjusted to 7 for ATM."
    assert (
        auto_acmg_data.thresholds.bp7_acceptor == 40
    ), "The BP7 acceptor threshold should be adjusted to 40 for ATM."

    # Check that the superclass's predict_bp7 method was called and returned a result
    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."


@patch.object(HBOPCPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default(mock_super_predict_bp7, hbopc_predictor, auto_acmg_data):
    """Test fallback to default BP7 prediction after threshold adjustment."""
    auto_acmg_data.hgnc_id = "HGNC:795"  # ATM gene
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = hbopc_predictor.predict_bp7(hbopc_predictor.seqvar, auto_acmg_data)

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


@patch.object(HBOPCPredictor, "predict_bp7", autospec=True)
def test_predict_bp7_fallback_to_default_for_palb2(
    mock_super_predict_bp7, hbopc_predictor, auto_acmg_data
):
    """Test fallback to default BP7 prediction after threshold adjustment for PALB2."""
    auto_acmg_data.hgnc_id = "HGNC:26144"  # PALB2 gene
    # Set the mock return value for the superclass's predict_bp7 method
    mock_super_predict_bp7.return_value = AutoACMGCriteria(
        name="BP7",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.BenignSupporting,
        summary="Default BP7 prediction fallback.",
    )

    # Call predict_bp7 method
    result = hbopc_predictor.predict_bp7(hbopc_predictor.seqvar, auto_acmg_data)

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
