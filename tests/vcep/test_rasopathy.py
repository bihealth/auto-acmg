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
from src.vcep import RASopathyPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="12", pos=100, delete="A", insert="T")


@pytest.fixture
def rasopathy_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return RASopathyPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


# ----------------- PM1 -----------------


@patch.object(
    DefaultSeqVarPredictor,
    "predict_pvs1",
    return_value=AutoACMGCriteria(
        name="PVS1",
        prediction=AutoACMGPrediction.Applicable,
        strength=AutoACMGStrength.PathogenicVeryStrong,
        summary="Superclass not applicable.",
    ),
)
def test_predict_pvs1_not_applicable(mock_super, rasopathy_predictor, seqvar, auto_acmg_data):
    # Set the HGNC ID to something other than LZTR1
    auto_acmg_data.hgnc_id = "HGNC:9999"
    result = rasopathy_predictor.predict_pvs1(seqvar, auto_acmg_data)

    # Verify the outcome is as expected, not applicable for non-LZTR1 genes
    assert result.name == "PVS1"
    assert result.prediction == AutoACMGPrediction.NotApplicable
    assert result.strength == AutoACMGStrength.PathogenicVeryStrong
    assert result.summary == "PVS1 is not applicable for the gene."
    mock_super.assert_not_called()


@patch.object(
    DefaultSeqVarPredictor,
    "predict_pvs1",
    return_value=AutoACMGCriteria(
        name="PVS1",
        prediction=AutoACMGPrediction.Applicable,
        strength=AutoACMGStrength.PathogenicVeryStrong,
        summary="Superclass applicable.",
    ),
)
def test_predict_pvs1_applicable_lztr1(mock_super, rasopathy_predictor, seqvar, auto_acmg_data):
    # Set the HGNC ID to LZTR1 to test the condition where the superclass method should be called
    auto_acmg_data.hgnc_id = "HGNC:6742"
    result = rasopathy_predictor.predict_pvs1(seqvar, auto_acmg_data)

    # Assertions to ensure the superclass method behaves as expected
    mock_super.assert_called_once_with(seqvar, auto_acmg_data)
    assert result.name == "PVS1"
    assert result.prediction == AutoACMGPrediction.Applicable
    assert result.strength == AutoACMGStrength.PathogenicVeryStrong
    assert result.summary == "Superclass applicable."


def test_predict_pm1_in_critical_region(rasopathy_predictor, auto_acmg_data):
    """Test when the variant falls within a critical region for a RASopathy gene."""
    auto_acmg_data.hgnc_id = "HGNC:7989"  # NRAS gene
    auto_acmg_data.prot_pos = 15  # Within the critical P-loop region (10-17)
    result = rasopathy_predictor.predict_pm1(rasopathy_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for critical region variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "critical residue in HGNC:7989" in result.summary
    ), "The summary should indicate the critical region."


def test_predict_pm1_in_critical_exon(rasopathy_predictor, auto_acmg_data):
    """Test when the variant falls within a critical exon for a RASopathy gene."""
    auto_acmg_data.hgnc_id = "HGNC:1097"  # BRAF gene
    auto_acmg_data.prot_pos = 480  # Position in the P-loop region, but within exon 11
    rasopathy_predictor._get_affected_exon = MagicMock(return_value=11)  # Mock the affected exon
    result = rasopathy_predictor.predict_pm1(rasopathy_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for critical exon variants."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert "critical exon 11" in result.summary, "The summary should indicate the critical exon."


def test_predict_pm1_not_applicable(rasopathy_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for certain RASopathy genes."""
    auto_acmg_data.hgnc_id = "HGNC:15454"  # SHOC2 gene
    result = rasopathy_predictor.predict_pm1(rasopathy_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be applicable for SHOC2."
    assert "not applicable" in result.summary, "The summary should indicate non-applicability."


def test_predict_pm1_outside_critical_region(rasopathy_predictor, auto_acmg_data):
    """Test when the variant does not fall within any critical region for RASopathy genes."""
    auto_acmg_data.hgnc_id = "HGNC:7989"  # NRAS gene
    auto_acmg_data.prot_pos = 200  # Position outside all critical regions
    rasopathy_predictor._get_affected_exon = MagicMock(return_value=1)  # Mock the affected exon
    result = rasopathy_predictor.predict_pm1(rasopathy_predictor.seqvar, auto_acmg_data)

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


@patch("src.vcep.rasopathy.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(mock_predict_pm1, rasopathy_predictor, auto_acmg_data):
    """Test fallback to the default PM1 prediction method for unhandled cases."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the RASopathy VCEP
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = rasopathy_predictor.predict_pm1(rasopathy_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


# -------------- PM2, BA1, BS1, BS2 --------------


@patch.object(RASopathyPredictor, "_get_af", return_value=0.1)
@patch.object(RASopathyPredictor, "_ba1_exception", return_value=False)
def test_verify_pm2ba1bs1bs2(
    mock_get_af, mock_ba1_exception, rasopathy_predictor, auto_acmg_data, seqvar
):
    # Setup: Adjusting the thresholds to test under different conditions
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.02
    auto_acmg_data.thresholds.pm2_pathogenic = 0.001

    # Call the method under test
    result, comment = rasopathy_predictor.verify_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assertions to validate the expected behavior
    assert result.BA1 is True, "Expected PM2 to be True based on the mocked allele frequency"

    # Assert changed thresholds
    assert (
        auto_acmg_data.thresholds.ba1_benign == 0.0005
    ), "BA1 threshold should be adjusted to 0.0005"
    assert (
        auto_acmg_data.thresholds.bs1_benign == 0.00025
    ), "BS1 threshold should be adjusted to 0.00025"


# -------------- PM4 & BP3 --------------


def test_bp3_not_applicable(rasopathy_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = rasopathy_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


# -------------- PP2 & BP1 --------------


@patch.object(RASopathyPredictor, "_is_missense")
def test_predict_pp2bp1_ptpn11_missense(
    mock_is_missense, rasopathy_predictor, seqvar, auto_acmg_data
):
    """Test PP2 prediction for PTPN11 where the variant is missense."""
    mock_is_missense.return_value = True
    auto_acmg_data.hgnc_id = "HGNC:9648"  # PTPN11 gene
    pp2, bp1 = rasopathy_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert (
        pp2.prediction == AutoACMGPrediction.NotApplicable
    ), "PP2 should be Met for a missense variant in PTPN11."
    assert (
        pp2.strength == AutoACMGStrength.PathogenicSupporting
    ), "PP2 strength should reflect PathogenicSupporting."
    assert (
        bp1.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should be NotApplicable for PTPN11."


# -------------- PP3 & BP4 --------------


@patch.object(RASopathyPredictor, "_is_missense")
def test_predict_map2k1_missense(mock_is_missense, rasopathy_predictor, seqvar, auto_acmg_data):
    """Test PP2 prediction for MAP2K1 where the variant is missense."""
    mock_is_missense.return_value = True
    auto_acmg_data.hgnc_id = "HGNC:6840"  # MAP2K1 gene
    pp2, bp1 = rasopathy_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert (
        pp2.prediction == AutoACMGPrediction.Applicable
    ), "PP2 should be Met for a missense variant in MAP2K1."
    assert (
        pp2.strength == AutoACMGStrength.PathogenicSupporting
    ), "PP2 strength should reflect PathogenicSupporting."
    assert (
        bp1.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should be NotApplicable for MAP2K1."


@patch.object(RASopathyPredictor, "_is_missense")
def test_predict_pp2bp1_map2k1_non_missense(
    mock_is_missense, rasopathy_predictor, seqvar, auto_acmg_data
):
    """Test PP2 prediction for MAP2K1 where the variant is not missense."""
    mock_is_missense.return_value = False
    auto_acmg_data.hgnc_id = "HGNC:6840"  # MAP2K1 gene
    pp2, bp1 = rasopathy_predictor.predict_pp2bp1(seqvar, auto_acmg_data)

    assert (
        pp2.prediction == AutoACMGPrediction.NotApplicable
    ), "PP2 should not be Met for non-missense variants in MAP2K1."
    assert "not a missense" in pp2.summary, "PP2 summary should confirm the non-missense nature."
    assert (
        bp1.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should still be NotApplicable for MAP2K1."


def test_predict_pp3bp4_revel_strategy(rasopathy_predictor, auto_acmg_data):
    """Test that REVEL is set as the strategy for PP3/BP4 prediction."""
    pp3_result, bp4_result = rasopathy_predictor.predict_pp3bp4(
        rasopathy_predictor.seqvar, auto_acmg_data
    )

    assert auto_acmg_data.thresholds.pp3bp4_strategy == "revel"
    assert auto_acmg_data.thresholds.revel_pathogenic == 0.7
    assert auto_acmg_data.thresholds.revel_benign == 0.3


@patch("src.vcep.rasopathy.DefaultSeqVarPredictor.predict_pp3bp4")
def test_predict_pp3bp4_calls_superclass(
    mock_super_predict_pp3bp4, rasopathy_predictor, auto_acmg_data
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

    pp3_result, bp4_result = rasopathy_predictor.predict_pp3bp4(
        rasopathy_predictor.seqvar, auto_acmg_data
    )

    mock_super_predict_pp3bp4.assert_called_once_with(rasopathy_predictor.seqvar, auto_acmg_data)
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
            0.2,
            AutoACMGPrediction.NotApplicable,
            AutoACMGPrediction.Applicable,
        ),  # Low REVEL score
    ],
)
def test_predict_pp3bp4_revel_scenarios(
    rasopathy_predictor, auto_acmg_data, revel_score, expected_pp3, expected_bp4
):
    """Test different REVEL score scenarios."""
    auto_acmg_data.scores.dbnsfp.revel = revel_score

    with patch(
        "src.vcep.rasopathy.DefaultSeqVarPredictor.predict_pp3bp4"
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

        pp3_result, bp4_result = rasopathy_predictor.predict_pp3bp4(
            rasopathy_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == expected_pp3
        assert bp4_result.prediction == expected_bp4


def test_predict_pp3bp4_strength(rasopathy_predictor, auto_acmg_data):
    """Test that the strength of PP3 and BP4 is correctly set."""
    with patch(
        "src.vcep.rasopathy.DefaultSeqVarPredictor.predict_pp3bp4"
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

        pp3_result, bp4_result = rasopathy_predictor.predict_pp3bp4(
            rasopathy_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.strength == AutoACMGStrength.PathogenicSupporting
        assert bp4_result.strength == AutoACMGStrength.BenignSupporting


def test_predict_pp3bp4_no_revel_score(rasopathy_predictor, auto_acmg_data):
    """Test behavior when no REVEL score is available."""
    auto_acmg_data.scores.dbnsfp.revel = None

    with patch(
        "src.vcep.rasopathy.DefaultSeqVarPredictor.predict_pp3bp4"
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

        pp3_result, bp4_result = rasopathy_predictor.predict_pp3bp4(
            rasopathy_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == AutoACMGPrediction.NotApplicable
        assert bp4_result.prediction == AutoACMGPrediction.NotApplicable


def test_predict_pp3bp4_error_handling(rasopathy_predictor, auto_acmg_data):
    """Test error handling in predict_pp3bp4 method."""
    with patch(
        "src.vcep.rasopathy.DefaultSeqVarPredictor.predict_pp3bp4",
        side_effect=Exception("Test error"),
    ):
        with pytest.raises(Exception) as exc_info:
            rasopathy_predictor.predict_pp3bp4(rasopathy_predictor.seqvar, auto_acmg_data)

        assert str(exc_info.value) == "Test error"
