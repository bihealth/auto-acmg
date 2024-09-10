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
from src.vcep import FamilialHypercholesterolemiaPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="19", pos=100, delete="A", insert="T")


@pytest.fixture
def familial_hypercholesterolemia_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return FamilialHypercholesterolemiaPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


def test_predict_pm1_in_critical_residue(familial_hypercholesterolemia_predictor, auto_acmg_data):
    """Test when variant affects a critical residue in LDLR."""
    auto_acmg_data.hgnc_id = "HGNC:6547"  # LDLR gene
    auto_acmg_data.prot_pos = 75  # Critical residue in LDLR
    result = familial_hypercholesterolemia_predictor.predict_pm1(
        familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for a variant affecting a critical residue."
    assert (
        "affects a critical cysteine residue" in result.summary
    ), "The summary should indicate the critical residue."


def test_predict_pm1_in_exon_4(familial_hypercholesterolemia_predictor, auto_acmg_data):
    """Test when variant affects the 4th exon in LDLR."""
    auto_acmg_data.hgnc_id = "HGNC:6547"  # LDLR gene
    with patch.object(FamilialHypercholesterolemiaPredictor, "_get_affected_exon", return_value=4):
        result = familial_hypercholesterolemia_predictor.predict_pm1(
            familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
        )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met for a variant affecting the 4th exon."
    assert (
        "affects the 4th exon" in result.summary
    ), "The summary should indicate the affected exon."


def test_predict_pm1_outside_critical_region(
    familial_hypercholesterolemia_predictor, auto_acmg_data
):
    """Test when variant does not fall within any critical region or exon 4 for LDLR."""
    auto_acmg_data.hgnc_id = "HGNC:6547"  # LDLR gene
    auto_acmg_data.prot_pos = 500  # Position outside critical residues and not in exon 4
    with patch.object(FamilialHypercholesterolemiaPredictor, "_get_affected_exon", return_value=5):
        result = familial_hypercholesterolemia_predictor.predict_pm1(
            familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
        )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met for a variant outside critical residues and exon 4."
    assert (
        "Variant does not meet the PM1 criteria for LDLR" in result.summary
    ), "The summary should indicate that no criteria were met."


@patch("src.vcep.familial_hypercholesterolemia.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, familial_hypercholesterolemia_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method for genes not in PM1_CLUSTER_LDLR."""
    auto_acmg_data.hgnc_id = "HGNC:99999"  # Gene not in PM1_CLUSTER_LDLR
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = familial_hypercholesterolemia_predictor.predict_pm1(
        familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met if the gene is not in the PM1_CLUSTER_LDLR."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should be from the default fallback."


def test_predict_pm1_edge_case_start_boundary(
    familial_hypercholesterolemia_predictor, auto_acmg_data
):
    """Test when variant falls exactly on the start boundary of a critical residue."""
    auto_acmg_data.hgnc_id = "HGNC:6547"  # LDLR gene
    auto_acmg_data.prot_pos = 27  # Start boundary of the critical residues (27)
    result = familial_hypercholesterolemia_predictor.predict_pm1(
        familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met when on the start boundary of a critical residue."
    assert (
        "affects a critical cysteine residue" in result.summary
    ), "The summary should indicate the critical residue."


def test_predict_pm1_edge_case_end_boundary(
    familial_hypercholesterolemia_predictor, auto_acmg_data
):
    """Test when variant falls exactly on the end boundary of a critical residue."""
    auto_acmg_data.hgnc_id = "HGNC:6547"  # LDLR gene
    auto_acmg_data.prot_pos = 711  # End boundary of the critical residues (711)
    result = familial_hypercholesterolemia_predictor.predict_pm1(
        familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
    )

    assert (
        result.prediction == AutoACMGPrediction.Applicable
    ), "PM1 should be met when on the end boundary of a critical residue."
    assert (
        "affects a critical cysteine residue" in result.summary
    ), "The summary should indicate the critical residue."


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
def test_predict_pm2ba1bs1bs2(
    mock_super_method, familial_hypercholesterolemia_predictor, auto_acmg_data, seqvar
):
    # Default thresholds
    auto_acmg_data.thresholds.pm2_pathogenic = 0.00001
    auto_acmg_data.thresholds.ba1_benign = 0.05
    auto_acmg_data.thresholds.bs1_benign = 0.01

    result = familial_hypercholesterolemia_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Assert that the thresholds were updated correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == 0.0002
    assert auto_acmg_data.thresholds.ba1_benign == 0.005
    assert auto_acmg_data.thresholds.bs1_benign == 0.002

    # Assert that the superclass method was called once with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Assert the response (optional, as we know it's mocked)
    assert all(
        c.name in ["PM2", "BA1", "BS1", "BS2"] for c in result
    ), "Unexpected criteria names returned"


def test_bp3_not_applicable(familial_hypercholesterolemia_predictor, seqvar, auto_acmg_data):
    """Test BP3 is not applicable for ACADVL as overridden."""
    result = familial_hypercholesterolemia_predictor._bp3_not_applicable(seqvar, auto_acmg_data)
    assert result is True, "BP3 should always be not applicable"


def test_predict_pp2bp1(familial_hypercholesterolemia_predictor, seqvar, auto_acmg_data):
    """Test predict_pp2bp1 for Familial Hypercholesterolemia predictor."""

    # Call the method under test
    pp2_result, bp1_result = familial_hypercholesterolemia_predictor.predict_pp2bp1(
        seqvar, auto_acmg_data
    )

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


def test_predict_pp3bp4_revel_strategy(familial_hypercholesterolemia_predictor, auto_acmg_data):
    """Test that REVEL is set as the strategy for PP3/BP4 prediction."""
    familial_hypercholesterolemia_predictor.predict_pp3bp4(
        familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
    )

    assert auto_acmg_data.thresholds.pp3bp4_strategy == "revel"


def test_predict_pp3bp4_revel_thresholds(familial_hypercholesterolemia_predictor, auto_acmg_data):
    """Test that REVEL thresholds are correctly set for PP3/BP4 prediction."""
    familial_hypercholesterolemia_predictor.predict_pp3bp4(
        familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
    )

    assert auto_acmg_data.thresholds.revel_pathogenic == 0.75
    assert auto_acmg_data.thresholds.revel_benign == 0.5


@patch("src.vcep.familial_hypercholesterolemia.DefaultSeqVarPredictor.predict_pp3bp4")
def test_predict_pp3bp4_calls_superclass(
    mock_super_predict_pp3bp4, familial_hypercholesterolemia_predictor, auto_acmg_data
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

    pp3_result, bp4_result = familial_hypercholesterolemia_predictor.predict_pp3bp4(
        familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
    )

    mock_super_predict_pp3bp4.assert_called_once_with(
        familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
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
            0.6,
            AutoACMGPrediction.NotApplicable,
            AutoACMGPrediction.NotApplicable,
        ),  # Intermediate REVEL score
        (
            0.4,
            AutoACMGPrediction.NotApplicable,
            AutoACMGPrediction.Applicable,
        ),  # Low REVEL score
    ],
)
def test_predict_pp3bp4_revel_scenarios(
    familial_hypercholesterolemia_predictor,
    auto_acmg_data,
    revel_score,
    expected_pp3,
    expected_bp4,
):
    """Test different REVEL score scenarios."""
    auto_acmg_data.scores.dbnsfp.revel = revel_score

    with patch(
        "src.vcep.familial_hypercholesterolemia.DefaultSeqVarPredictor.predict_pp3bp4"
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

        pp3_result, bp4_result = familial_hypercholesterolemia_predictor.predict_pp3bp4(
            familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == expected_pp3
        assert bp4_result.prediction == expected_bp4


def test_predict_pp3bp4_strength(familial_hypercholesterolemia_predictor, auto_acmg_data):
    """Test that the strength of PP3 and BP4 is correctly set."""
    with patch(
        "src.vcep.familial_hypercholesterolemia.DefaultSeqVarPredictor.predict_pp3bp4"
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

        pp3_result, bp4_result = familial_hypercholesterolemia_predictor.predict_pp3bp4(
            familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.strength == AutoACMGStrength.PathogenicSupporting
        assert bp4_result.strength == AutoACMGStrength.BenignSupporting


def test_predict_pp3bp4_no_revel_score(familial_hypercholesterolemia_predictor, auto_acmg_data):
    """Test behavior when no REVEL score is available."""
    auto_acmg_data.scores.dbnsfp.revel = None

    with patch(
        "src.vcep.familial_hypercholesterolemia.DefaultSeqVarPredictor.predict_pp3bp4"
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

        pp3_result, bp4_result = familial_hypercholesterolemia_predictor.predict_pp3bp4(
            familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
        )

        assert pp3_result.prediction == AutoACMGPrediction.NotApplicable
        assert bp4_result.prediction == AutoACMGPrediction.NotApplicable


def test_predict_pp3bp4_error_handling(familial_hypercholesterolemia_predictor, auto_acmg_data):
    """Test error handling in predict_pp3bp4 method."""
    with patch(
        "src.vcep.familial_hypercholesterolemia.DefaultSeqVarPredictor.predict_pp3bp4",
        side_effect=Exception("Test error"),
    ):
        with pytest.raises(Exception) as exc_info:
            familial_hypercholesterolemia_predictor.predict_pp3bp4(
                familial_hypercholesterolemia_predictor.seqvar, auto_acmg_data
            )

        assert str(exc_info.value) == "Test error"
