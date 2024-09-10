from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import PP2BP1, AutoACMGPrediction, AutoACMGStrength
from src.defs.exceptions import AlgorithmError, InvalidAPIResposeError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.auto_pp2_bp1 import AutoPP2BP1


@pytest.fixture
def seqvar():
    return MagicMock(spec=SeqVar)


@pytest.fixture
def auto_pp2bp1():
    return AutoPP2BP1()


# ============= _get_missense_vars ==================


@pytest.mark.skip(reason="Annonars is not mocked properly")
@patch("src.utils.AutoACMGHelper.annonars_client.get_variant_from_range")
def test_get_missense_vars_success(mock_get_variant, auto_pp2bp1, seqvar):
    # Setup the mock response with pathogenic and benign variants
    pathogenic_variant = MagicMock()
    pathogenic_variant.records = [
        MagicMock(
            classifications=MagicMock(germlineClassification=MagicMock(description="Pathogenic"))
        )
    ]
    benign_variant = MagicMock()
    benign_variant.records = [
        MagicMock(classifications=MagicMock(germlineClassification=MagicMock(description="Benign")))
    ]
    response = MagicMock(clinvar=[pathogenic_variant, benign_variant])
    mock_get_variant.return_value = response

    pathogenic, benign, total = auto_pp2bp1._get_missense_vars(seqvar, 100, 200)

    assert pathogenic == 1
    assert benign == 1
    assert total == 2


@pytest.mark.skip(reason="Annonars is not mocked properly")
@patch("src.utils.AutoACMGHelper.annonars_client.get_variant_from_range")
def test_get_missense_vars_empty_response(mock_get_variant, auto_pp2bp1, seqvar):
    # Setup the mock to return an empty response
    mock_get_variant.return_value = MagicMock(clinvar=[])

    with pytest.raises(InvalidAPIResposeError):
        auto_pp2bp1._get_missense_vars(seqvar, 100, 200)


def test_get_missense_vars_invalid_range(auto_pp2bp1, seqvar):
    # Test error handling when the end position is less than the start position
    with pytest.raises(AlgorithmError) as excinfo:
        auto_pp2bp1._get_missense_vars(seqvar, 200, 100)
    assert "End position is less than the start position" in str(excinfo.value)


# =============== _is_missense ==================


@pytest.fixture
def var_data_missense():
    # Simulate data structure for consequences with missense in both 'cadd' and 'mehari'
    cadd = MagicMock(cadd="missense_variant other_info")
    mehari = ["non_effect", "missense_variant", "other_effect"]
    return MagicMock(consequence=MagicMock(cadd=cadd, mehari=mehari))


@pytest.fixture
def var_data_non_missense():
    # Simulate data structure for consequences without any missense
    cadd = MagicMock(cadd="other_variant_info")
    mehari = ["non_effect", "other_effect"]
    return MagicMock(consequence=MagicMock(cadd=cadd, mehari=mehari))


@pytest.fixture
def mock_var_data_missense():
    var_data = MagicMock()
    var_data.consequence = MagicMock()
    var_data.consequence.cadd = "some_effect missense other_effect"
    var_data.consequence.mehari = ["some_effect", "missense", "other_effect"]
    return var_data


@pytest.fixture
def mock_var_data_non_missense():
    var_data = MagicMock()
    var_data.consequence = MagicMock()
    var_data.consequence.cadd = "some_effect other_effect"
    var_data.consequence.mehari = ["some_effect", "other_effect"]
    return var_data


def test_is_missense_true_from_cadd(var_data_missense):
    # Using the var_data_missense fixture to test missense detection from 'cadd'
    assert AutoPP2BP1._is_missense(var_data_missense) is True


def test_is_missense_true_from_mehari(var_data_missense):
    # Using the var_data_missense fixture to test missense detection from 'mehari'
    assert AutoPP2BP1._is_missense(var_data_missense) is True


def test_is_missense_false(var_data_non_missense):
    # Using the var_data_non_missense fixture to confirm no false positives for missense detection
    assert AutoPP2BP1._is_missense(var_data_non_missense) is False


# =============== verify_pp2bp1 ==================


@pytest.fixture
def seqvar_mitochondrial():
    return SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="MT",
        pos=1000,
        delete="A",
        insert="G",
    )


@pytest.fixture
def seqvar_non_mitochondrial():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=1000, delete="A", insert="G")


@pytest.fixture
def var_data_verify():
    cadd = MagicMock()
    mehari = ["missense"]
    return MagicMock(consequence=MagicMock(cadd=cadd, mehari=mehari), cds_start=950, cds_end=1050)


@pytest.fixture
def mock_thresholds():
    return MagicMock(pp2bp1_pathogenic=0.2, pp2bp1_benign=0.1)


def test_verify_pp2bp1_mitochondrial(auto_pp2bp1, seqvar_mitochondrial, var_data_verify):
    """Test verify_pp2bp1 with mitochondrial variant."""
    result, comment = auto_pp2bp1.verify_pp2bp1(seqvar_mitochondrial, var_data_verify)
    assert result.PP2 is False and result.BP1 is False
    assert "mitochondrial DNA" in comment


def test_verify_pp2bp1_non_missense(auto_pp2bp1, seqvar_non_mitochondrial, var_data_verify):
    """Test verify_pp2bp1 with non-missense variant."""
    with patch.object(auto_pp2bp1, "_is_missense", return_value=False):
        result, comment = auto_pp2bp1.verify_pp2bp1(seqvar_non_mitochondrial, var_data_verify)
        assert result.PP2 is False and result.BP1 is False
        assert "not a missense variant" in comment


def test_verify_pp2bp1_with_misZ(
    auto_pp2bp1, seqvar_non_mitochondrial, var_data_verify, mock_thresholds
):
    """Test verify_pp2bp1 with misZ score."""
    var_data_verify.scores = MagicMock(misZ=0.3)
    var_data_verify.thresholds = mock_thresholds
    result, comment = auto_pp2bp1.verify_pp2bp1(seqvar_non_mitochondrial, var_data_verify)
    assert result.PP2 is True
    assert "Z-score" in comment


def test_verify_pp2bp1_without_misZ(
    auto_pp2bp1, seqvar_non_mitochondrial, var_data_verify, mock_thresholds
):
    """Test verify_pp2bp1 without misZ score."""
    var_data_verify.scores = MagicMock(misZ=None)
    var_data_verify.thresholds = mock_thresholds
    with patch.object(auto_pp2bp1, "_get_missense_vars", return_value=(14, 1, 15)):
        result, comment = auto_pp2bp1.verify_pp2bp1(seqvar_non_mitochondrial, var_data_verify)
        assert result.PP2 is True
        assert result.BP1 is False
        assert "Counting missense variants" in comment


def test_verify_pp2bp1_error(auto_pp2bp1, seqvar_non_mitochondrial, var_data_verify):
    """Test verify_pp2bp1 with error from _get_missense_vars method."""
    var_data_verify.scores = MagicMock(misZ=None)
    var_data_verify.thresholds = mock_thresholds
    with patch.object(
        auto_pp2bp1,
        "_get_missense_vars",
        side_effect=InvalidAPIResposeError("API Error"),
    ):
        result, comment = auto_pp2bp1.verify_pp2bp1(seqvar_non_mitochondrial, var_data_verify)
        assert result is None
        assert "Error occurred" in comment


# =============== predict_pp2bp1 ==================


@pytest.fixture
def var_data_pp2bp1():
    # Mock the thresholds for pathogenic and benign missense variants
    thresholds = MagicMock(
        pp2bp1_pathogenic=0.2,  # Hypothetical threshold for pathogenic missense variants
        pp2bp1_benign=0.1,  # Hypothetical threshold for benign missense variants
    )

    # Mock the consequence information in a similar fashion to your spliceAI example
    consequence = MagicMock(
        cadd={"missense": True},  # Simulates CADD indicating a missense variant
        mehari=["missense"],  # Simulates other tools also indicating a missense variant
    )

    # Combine into a mocked AutoACMGData structure
    return MagicMock(consequence=consequence, thresholds=thresholds)


@pytest.fixture
def pp2bp1_result_met():
    return (
        PP2BP1(
            PP2=True,
            BP1=False,
            PP2_strength=AutoACMGStrength.PathogenicSupporting,
            BP1_strength=AutoACMGStrength.BenignSupporting,
        ),
        "PP2 criteria met",
    )


@pytest.fixture
def pp2bp1_result_not_met():
    return (
        PP2BP1(
            PP2=False,
            BP1=True,
            PP2_strength=AutoACMGStrength.PathogenicSupporting,
            BP1_strength=AutoACMGStrength.BenignSupporting,
        ),
        "BP1 criteria not met",
    )


@pytest.fixture
def pp2bp1_result_failed():
    return None, "Error during prediction."


@patch("src.seqvar.auto_pp2_bp1.AutoPP2BP1.verify_pp2bp1")
def test_predict_pp2bp1_met(mock_verify, auto_pp2bp1, seqvar, var_data_pp2bp1, pp2bp1_result_met):
    """Test predict_pp2bp1 where both criteria are correctly met."""
    mock_verify.return_value = pp2bp1_result_met
    result = auto_pp2bp1.predict_pp2bp1(seqvar, var_data_pp2bp1)
    assert result[0].prediction == AutoACMGPrediction.Applicable
    assert result[0].strength == AutoACMGStrength.PathogenicSupporting
    assert "PP2 criteria met" in result[0].summary
    assert result[1].prediction == AutoACMGPrediction.NotApplicable
    assert result[1].strength == AutoACMGStrength.BenignSupporting
    assert "PP2 criteria met" in result[1].summary


@patch("src.seqvar.auto_pp2_bp1.AutoPP2BP1.verify_pp2bp1")
def test_predict_pp2bp1_not_met(
    mock_verify, auto_pp2bp1, seqvar, var_data_pp2bp1, pp2bp1_result_not_met
):
    """Test predict_pp2bp1 where the criteria are not met."""
    mock_verify.return_value = pp2bp1_result_not_met
    result = auto_pp2bp1.predict_pp2bp1(seqvar, var_data_pp2bp1)
    assert result[0].prediction == AutoACMGPrediction.NotApplicable
    assert result[0].strength == AutoACMGStrength.PathogenicSupporting
    assert "BP1 criteria not met" in result[0].summary
    assert result[1].prediction == AutoACMGPrediction.Applicable
    assert result[1].strength == AutoACMGStrength.BenignSupporting
    assert "BP1 criteria not met" in result[1].summary


@patch("src.seqvar.auto_pp2_bp1.AutoPP2BP1.verify_pp2bp1")
def test_predict_pp2bp1_failed(
    mock_verify, auto_pp2bp1, seqvar, var_data_pp2bp1, pp2bp1_result_failed
):
    """Test predict_pp2bp1 when there's a failure to evaluate the criteria."""
    mock_verify.return_value = pp2bp1_result_failed
    result = auto_pp2bp1.predict_pp2bp1(seqvar, var_data_pp2bp1)
    assert result[0].prediction == AutoACMGPrediction.Failed
    assert result[0].strength == AutoACMGStrength.PathogenicSupporting
    assert "Error during prediction." in result[0].summary
    assert result[1].prediction == AutoACMGPrediction.Failed
    assert result[1].strength == AutoACMGStrength.BenignSupporting
    assert "Error during prediction." in result[1].summary
