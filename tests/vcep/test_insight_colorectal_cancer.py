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
from src.vcep import InsightColorectalCancerPredictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="5", pos=100, delete="A", insert="T")


@pytest.fixture
def insight_colorectal_cancer_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return InsightColorectalCancerPredictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGSeqVarData()


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_overrides(
    mock_super_verify, insight_colorectal_cancer_predictor, seqvar, auto_acmg_data
):
    """Test that the overridden verify_ps1pm5 method in InsightColorectalCancerPredictor works correctly."""
    # Set up the mock to return PS1 and PM5 as applicable initially
    mock_super_verify.return_value = (PS1PM5(PS1=True, PM5=True), "Initial evaluation")

    # Setup the data
    auto_acmg_data.consequence = MagicMock(mehari=["missense_variant"])
    auto_acmg_data.thresholds = MagicMock(
        spliceAI_acceptor_gain=0.5,
        spliceAI_acceptor_loss=0.5,
        spliceAI_donor_gain=0.5,
        spliceAI_donor_loss=0.5,
    )
    auto_acmg_data.scores = MagicMock(
        cadd=MagicMock(spliceAI_acceptor_gain=0.3, spliceAI_donor_gain=0.3)
    )

    # Run the method under test
    prediction, comment = insight_colorectal_cancer_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that the splicing effect leads to overriding PS1 and PM5 as not applicable
    assert not prediction.PS1, "PS1 should be marked as not applicable due to splicing effect."
    assert not prediction.PM5, "PM5 should be marked as not applicable due to splicing effect."
    assert "Variant affects splicing" in comment, "Comment should note the splicing effect."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
    ), "SpliceAI acceptor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
    ), "SpliceAI acceptor loss threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
    ), "SpliceAI donor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2
    ), "SpliceAI donor loss threshold should be adjusted to 0.2"


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_no_splicing_effect(
    mock_super_verify, insight_colorectal_cancer_predictor, seqvar, auto_acmg_data
):
    """Test that PS1 and PM5 remain applicable when there's no splicing effect."""
    # Set up the mock to return PS1 and PM5 as applicable initially
    mock_super_verify.return_value = (PS1PM5(PS1=True, PM5=True), "Initial evaluation")

    # Setup the data with no splicing effect
    auto_acmg_data.consequence = MagicMock(mehari=["missense_variant"])
    auto_acmg_data.thresholds = MagicMock(
        spliceAI_acceptor_gain=0.5,
        spliceAI_acceptor_loss=0.5,
        spliceAI_donor_gain=0.5,
        spliceAI_donor_loss=0.5,
    )
    auto_acmg_data.scores = MagicMock(
        cadd=MagicMock(
            spliceAI_acceptor_gain=0.1,
            spliceAI_acceptor_loss=0.1,
            spliceAI_donor_gain=0.1,
            spliceAI_donor_loss=0.1,
        )
    )

    # Run the method under test
    prediction, comment = insight_colorectal_cancer_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that PS1 and PM5 remain applicable
    assert prediction.PS1, "PS1 should remain applicable when there's no splicing effect."
    assert prediction.PM5, "PM5 should remain applicable when there's no splicing effect."
    assert "Initial evaluation" in comment, "Comment should reflect the initial evaluation."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
    ), "SpliceAI acceptor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
    ), "SpliceAI acceptor loss threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
    ), "SpliceAI donor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2
    ), "SpliceAI donor loss threshold should be adjusted to 0.2"


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_non_missense(
    mock_super_verify, insight_colorectal_cancer_predictor, seqvar, auto_acmg_data
):
    """Test that PS1 and PM5 remain as per superclass for non-missense variants."""
    # Set up the mock to return PS1 and PM5 as not applicable initially
    mock_super_verify.return_value = (
        PS1PM5(PS1=False, PM5=False),
        "Not applicable for non-missense",
    )

    # Setup the data with a non-missense variant
    auto_acmg_data.consequence = MagicMock(mehari=["synonymous_variant"])

    # Run the method under test
    prediction, comment = insight_colorectal_cancer_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    # Check that PS1 and PM5 remain as per superclass evaluation
    assert not prediction.PS1, "PS1 should remain not applicable for non-missense variants."
    assert not prediction.PM5, "PM5 should remain not applicable for non-missense variants."
    assert (
        "Not applicable for non-missense" in comment
    ), "Comment should reflect superclass evaluation."

    # Ensure that the mock of the superclass method is called to simulate the inherited behavior
    mock_super_verify.assert_called_once_with(seqvar, auto_acmg_data)

    # Check that the thresholds were adjusted
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.2
    ), "SpliceAI acceptor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.2
    ), "SpliceAI acceptor loss threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_gain == 0.2
    ), "SpliceAI donor gain threshold should be adjusted to 0.2"
    assert (
        auto_acmg_data.thresholds.spliceAI_donor_loss == 0.2
    ), "SpliceAI donor loss threshold should be adjusted to 0.2"


@patch.object(DefaultSeqVarPredictor, "verify_ps1pm5")
def test_verify_ps1pm5_exception_handling(
    mock_super_verify, insight_colorectal_cancer_predictor, seqvar, auto_acmg_data
):
    """Test verify_ps1pm5 method exception handling."""
    mock_super_verify.side_effect = Exception("Test exception")

    with pytest.raises(Exception) as exc_info:
        insight_colorectal_cancer_predictor.verify_ps1pm5(seqvar, auto_acmg_data)

    assert "Test exception" in str(exc_info.value), "Should raise the original exception"


@pytest.fixture
def auto_acmg_data_apc():
    data = AutoACMGSeqVarData(hgnc_id="HGNC:583")
    data.consequence = MagicMock(mehari=["missense_variant"], cadd="missense")
    data.cds_start = 200
    data.cds_end = 300
    data.thresholds.pp2bp1_benign = 0.05
    return data


def test_predict_pm1_not_applicable_apc(insight_colorectal_cancer_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for APC."""
    auto_acmg_data.hgnc_id = "HGNC:583"  # APC gene
    result = insight_colorectal_cancer_predictor.predict_pm1(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for APC."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not applicable for HGNC:583" in result.summary
    ), "The summary should indicate PM1 is not applicable."


def test_predict_pm1_not_applicable_mlh1(insight_colorectal_cancer_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for MLH1."""
    auto_acmg_data.hgnc_id = "HGNC:7127"  # MLH1 gene
    result = insight_colorectal_cancer_predictor.predict_pm1(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for MLH1."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not applicable for HGNC:7127" in result.summary
    ), "The summary should indicate PM1 is not applicable."


def test_predict_pm1_not_applicable_msh2(insight_colorectal_cancer_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for MSH2."""
    auto_acmg_data.hgnc_id = "HGNC:7325"  # MSH2 gene
    result = insight_colorectal_cancer_predictor.predict_pm1(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for MSH2."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not applicable for HGNC:7325" in result.summary
    ), "The summary should indicate PM1 is not applicable."


def test_predict_pm1_not_applicable_msh6(insight_colorectal_cancer_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for MSH6."""
    auto_acmg_data.hgnc_id = "HGNC:7329"  # MSH6 gene
    result = insight_colorectal_cancer_predictor.predict_pm1(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for MSH6."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not applicable for HGNC:7329" in result.summary
    ), "The summary should indicate PM1 is not applicable."


def test_predict_pm1_not_applicable_pms2(insight_colorectal_cancer_predictor, auto_acmg_data):
    """Test when PM1 is not applicable for PMS2."""
    auto_acmg_data.hgnc_id = "HGNC:9122"  # PMS2 gene
    result = insight_colorectal_cancer_predictor.predict_pm1(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should be not applicable for PMS2."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "PM1 is not applicable for HGNC:9122" in result.summary
    ), "The summary should indicate PM1 is not applicable."


@patch("src.vcep.insight_colorectal_cancer.DefaultSeqVarPredictor.predict_pm1")
def test_predict_pm1_fallback_to_default(
    mock_predict_pm1, insight_colorectal_cancer_predictor, auto_acmg_data
):
    """Test fallback to the default PM1 prediction method if logic changes."""
    auto_acmg_data.hgnc_id = "HGNC:9999"  # Gene not in the specific logic
    mock_predict_pm1.return_value = AutoACMGCriteria(
        name="PM1",
        prediction=AutoACMGPrediction.NotApplicable,
        strength=AutoACMGStrength.PathogenicModerate,
        summary="Default PM1 prediction fallback.",
    )
    result = insight_colorectal_cancer_predictor.predict_pm1(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotApplicable
    ), "PM1 should not be met in the default fallback."
    assert (
        "Default PM1 prediction fallback." in result.summary
    ), "The summary should indicate the default fallback."


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
        ("HGNC:583", 0.000003, 0.001, 0.0001),
        ("HGNC:7329", 0.00002, 0.0022, 0.00022),
        ("HGNC:9122", 0.00002, 0.0028, 0.0001),
    ],
)
def test_predict_pm2ba1bs1bs2_specific_genes(
    mock_super_method,
    insight_colorectal_cancer_predictor,
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
    insight_colorectal_cancer_predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_data)

    # Validate thresholds are set correctly
    assert auto_acmg_data.thresholds.pm2_pathogenic == expected_pm2
    assert auto_acmg_data.thresholds.ba1_benign == expected_ba1
    assert auto_acmg_data.thresholds.bs1_benign == expected_bs1

    # Check that the superclass method was called with the modified var_data
    mock_super_method.assert_called_once_with(seqvar, auto_acmg_data)

    # Reset mock for the next iteration
    mock_super_method.reset_mock()


def test_predict_pm4bp3_not_applicable(insight_colorectal_cancer_predictor, seqvar, auto_acmg_data):
    """Test that PM4 and BP3 are marked as Not Applicable for the brain malformations VCEP."""
    pm4_result, bp3_result = insight_colorectal_cancer_predictor.predict_pm4bp3(
        seqvar, auto_acmg_data
    )

    # Check PM4 result
    assert isinstance(
        pm4_result, AutoACMGCriteria
    ), "The PM4 result should be of type AutoACMGCriteria."
    assert pm4_result.prediction == AutoACMGPrediction.NotApplicable, "PM4 should be NotApplicable."
    assert (
        pm4_result.strength == AutoACMGStrength.PathogenicModerate
    ), "PM4 strength should be PathogenicModerate."
    assert (
        "PM4 is not applicable" in pm4_result.summary
    ), "The summary should indicate PM4 is not applicable."

    # Check BP3 result
    assert isinstance(
        bp3_result, AutoACMGCriteria
    ), "The BP3 result should be of type AutoACMGCriteria."
    assert bp3_result.prediction == AutoACMGPrediction.NotApplicable, "BP3 should be NotApplicable."
    assert (
        bp3_result.strength == AutoACMGStrength.BenignSupporting
    ), "BP3 strength should be BenignSupporting."
    assert (
        "BP3 is not applicable" in bp3_result.summary
    ), "The summary should indicate BP3 is not applicable."


@patch.object(InsightColorectalCancerPredictor, "_is_missense")
@patch.object(InsightColorectalCancerPredictor, "_get_missense_vars")
def test_predict_pp2bp1_apc_missense_high_benign_ratio(
    mock_get_missense_vars,
    mock_is_missense,
    insight_colorectal_cancer_predictor,
    seqvar,
    auto_acmg_data_apc,
):
    """Test PP2 and BP1 prediction where the APC variant is missense with a high benign ratio."""
    mock_is_missense.return_value = True
    mock_get_missense_vars.return_value = (None, 100, 1000)  # 10% benign ratio

    pp2, bp1 = insight_colorectal_cancer_predictor.predict_pp2bp1(seqvar, auto_acmg_data_apc)

    assert pp2.prediction == AutoACMGPrediction.NotApplicable, "PP2 should be NotApplicable."
    assert (
        bp1.prediction == AutoACMGPrediction.Applicable
    ), "BP1 should be Met due to high benign ratio."
    assert (
        "Benign ratio" in bp1.summary and "is met" in bp1.summary
    ), "BP1 summary should confirm criteria met due to benign ratio."


@patch.object(InsightColorectalCancerPredictor, "_is_missense")
@patch.object(InsightColorectalCancerPredictor, "_get_missense_vars")
def test_predict_pp2bp1_apc_missense_low_benign_ratio(
    mock_get_missense_vars,
    mock_is_missense,
    insight_colorectal_cancer_predictor,
    seqvar,
    auto_acmg_data_apc,
):
    """Test PP2 and BP1 prediction where the APC variant is missense with a low benign ratio."""
    mock_is_missense.return_value = True
    mock_get_missense_vars.return_value = (None, 1, 1000)  # 0.1% benign ratio

    pp2, bp1 = insight_colorectal_cancer_predictor.predict_pp2bp1(seqvar, auto_acmg_data_apc)

    assert pp2.prediction == AutoACMGPrediction.NotApplicable, "PP2 should be NotApplicable."
    assert (
        bp1.prediction == AutoACMGPrediction.NotApplicable
    ), "BP1 should not be Met due to low benign ratio."
    assert (
        "Benign ratio" in bp1.summary and "is not met" in bp1.summary
    ), "BP1 summary should confirm criteria not met due to benign ratio."


def test_predict_bp7_threshold_adjustment(insight_colorectal_cancer_predictor, auto_acmg_data):
    """Test that the BP7 donor and acceptor thresholds are correctly adjusted."""
    auto_acmg_data.thresholds.bp7_donor = 1  # Initial donor threshold value
    auto_acmg_data.thresholds.bp7_acceptor = 2  # Initial acceptor threshold value

    # Call predict_bp7 method
    result = insight_colorectal_cancer_predictor.predict_bp7(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
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
    mock_super_predict_bp7, insight_colorectal_cancer_predictor, auto_acmg_data
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
    result = insight_colorectal_cancer_predictor.predict_bp7(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

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


def test_verify_pp3bp4_thresholds(insight_colorectal_cancer_predictor, auto_acmg_data):
    """Test that the thresholds for PP3/BP4 prediction are correctly set."""
    insight_colorectal_cancer_predictor.verify_pp3bp4(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert auto_acmg_data.thresholds.spliceAI_acceptor_gain == 0.1
    assert auto_acmg_data.thresholds.spliceAI_acceptor_loss == 0.1
    assert auto_acmg_data.thresholds.spliceAI_donor_gain == 0.1
    assert auto_acmg_data.thresholds.spliceAI_donor_loss == 0.1


@patch.object(InsightColorectalCancerPredictor, "_is_pathogenic_score")
@patch.object(InsightColorectalCancerPredictor, "_is_benign_score")
@patch.object(InsightColorectalCancerPredictor, "_affect_spliceAI")
def test_verify_pp3bp4_prediction_logic(
    mock_affect_spliceAI,
    mock_is_benign_score,
    mock_is_pathogenic_score,
    insight_colorectal_cancer_predictor,
    auto_acmg_data,
):
    """Test the prediction logic for PP3 and BP4."""
    mock_is_pathogenic_score.return_value = True
    mock_is_benign_score.return_value = False
    mock_affect_spliceAI.side_effect = [
        True,
        False,
    ]  # First call True, second call False

    prediction, comment = insight_colorectal_cancer_predictor.verify_pp3bp4(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 is True
    assert prediction.BP4 is False
    assert "MetaRNN score" in comment
    assert "BayesDel_noAF score" in comment


@pytest.mark.parametrize(
    "metaRNN_score, bayesDel_score, spliceAI_scores, expected_pp3, expected_bp4",
    [
        (0.9, 0.9, [0.3, 0.3, 0.3, 0.3], True, False),  # High pathogenic scores
        (0.1, 0.1, [0.1, 0.1, 0.1, 0.1], False, True),  # High benign scores
        (0.5, 0.5, [0.15, 0.15, 0.15, 0.15], False, False),  # Intermediate scores
        (0.9, 0.1, [0.3, 0.3, 0.3, 0.3], True, False),  # Mixed scores, high spliceAI
        # (0.1, 0.9, [0.1, 0.1, 0.1, 0.1], False, True),  # Mixed scores, low spliceAI
    ],
)
def test_verify_pp3bp4_various_scenarios(
    insight_colorectal_cancer_predictor,
    auto_acmg_data,
    metaRNN_score,
    bayesDel_score,
    spliceAI_scores,
    expected_pp3,
    expected_bp4,
):
    """Test different scenarios for PP3 and BP4 prediction."""
    auto_acmg_data.scores.dbnsfp.metaRNN = metaRNN_score
    auto_acmg_data.scores.dbnsfp.bayesDel_noAF = bayesDel_score
    auto_acmg_data.scores.cadd.spliceAI_acceptor_gain = spliceAI_scores[0]
    auto_acmg_data.scores.cadd.spliceAI_acceptor_loss = spliceAI_scores[1]
    auto_acmg_data.scores.cadd.spliceAI_donor_gain = spliceAI_scores[2]
    auto_acmg_data.scores.cadd.spliceAI_donor_loss = spliceAI_scores[3]

    prediction, _ = insight_colorectal_cancer_predictor.verify_pp3bp4(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert prediction.PP3 == expected_pp3
    assert prediction.BP4 == expected_bp4


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_missing_scores(insight_colorectal_cancer_predictor, auto_acmg_data):
    """Test behavior when scores are missing."""
    auto_acmg_data.scores.dbnsfp.metaRNN = None
    auto_acmg_data.scores.dbnsfp.bayesDel_noAF = None

    prediction, comment = insight_colorectal_cancer_predictor.verify_pp3bp4(
        insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
    )

    assert prediction is None
    assert "An error occurred during prediction" in comment


@pytest.mark.skip(reason="Fix it")
def test_verify_pp3bp4_error_handling(insight_colorectal_cancer_predictor, auto_acmg_data):
    """Test error handling in verify_pp3bp4 method."""
    with patch.object(
        InsightColorectalCancerPredictor,
        "_is_pathogenic_score",
        side_effect=Exception("Test error"),
    ):
        prediction, comment = insight_colorectal_cancer_predictor.verify_pp3bp4(
            insight_colorectal_cancer_predictor.seqvar, auto_acmg_data
        )

        assert prediction is None
        assert "An error occurred during prediction" in comment
        assert "Test error" in comment
