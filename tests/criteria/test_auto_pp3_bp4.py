from unittest.mock import MagicMock, patch

import pytest

from src.criteria.auto_pp3_bp4 import AutoPP3BP4
from src.defs.auto_acmg import PP3BP4, AutoACMGPrediction, AutoACMGStrength
from src.defs.exceptions import AutoAcmgBaseException, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


@pytest.fixture
def auto_pp3bp4():
    return AutoPP3BP4()


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


# =========== _splice_variant ===========


@pytest.fixture
def var_data_splice():
    consequence = MagicMock(
        cadd={"splice": True},  # Simulates CADD indicating a splice variant
        mehari=["splice", "other_effect"],  # Simulates other tools also indicating a splice variant
    )
    return MagicMock(consequence=consequence)


@pytest.fixture
def var_data_no_splice():
    consequence = MagicMock(cadd={}, mehari=[])  # No splice key  # No splice-related terms
    return MagicMock(consequence=consequence)


def test_splice_variant_positive(auto_pp3bp4: AutoPP3BP4, var_data_splice: MagicMock):
    """Test that _splice_variant correctly identifies a splice variant."""
    assert (
        auto_pp3bp4._is_splice_variant(var_data_splice) is True
    ), "Should return True when splice indicators are present in the data."


def test_splice_variant_negative(auto_pp3bp4: AutoPP3BP4, var_data_no_splice: MagicMock):
    """Test that _splice_variant correctly identifies non-splice variants."""
    assert (
        auto_pp3bp4._is_splice_variant(var_data_no_splice) is False
    ), "Should return False when no splice indicators are present in the data."


def test_splice_variant_with_other_effects(auto_pp3bp4: AutoPP3BP4, var_data_splice: MagicMock):
    """Test _splice_variant with other non-splice related effects."""
    # Adjust the mock to include other non-splice related effects
    var_data_splice.consequence.mehari.append("non_splice_effect")
    assert (
        auto_pp3bp4._is_splice_variant(var_data_splice) is True
    ), "Should still return True as long as one splice indicator is present."


# =========== _is_pathogenic_score ===========


@pytest.fixture
def var_data_pathogenic():
    thresholds = MagicMock(
        metaRNN_pathogenic=0.5,
        bayesDel_noAF_pathogenic=0.5,
    )
    scores_dbnsfp = MagicMock(
        metaRNN=0.7,  # Above pathogenic threshold
        bayesDel_noAF=0.6,  # Above pathogenic threshold
    )
    scores = MagicMock(dbnsfp=scores_dbnsfp)
    return MagicMock(scores=scores, thresholds=thresholds)


@pytest.fixture
def var_data_non_pathogenic():
    thresholds = MagicMock(
        metaRNN_pathogenic=0.5,
        bayesDel_noAF_pathogenic=0.5,
    )
    scores_dbnsfp = MagicMock(
        metaRNN=0.4,  # Below pathogenic threshold
        bayesDel_noAF=0.3,  # Below pathogenic threshold
    )
    scores = MagicMock(dbnsfp=scores_dbnsfp)
    return MagicMock(scores=scores, thresholds=thresholds)


@pytest.fixture
def var_data_missing_pathogenic_scores():
    thresholds = MagicMock(
        metaRNN_pathogenic=0.5,
        bayesDel_noAF_pathogenic=0.5,
    )
    scores_dbnsfp = MagicMock(
        metaRNN=None,  # Missing score
        bayesDel_noAF=None,  # Missing score
    )
    scores = MagicMock(dbnsfp=scores_dbnsfp)
    return MagicMock(scores=scores, thresholds=thresholds)


def test_is_pathogenic_score_true(auto_pp3bp4, var_data_pathogenic):
    """Test when pathogenic scores are above the thresholds."""
    assert (
        auto_pp3bp4._is_pathogenic_score(var_data_pathogenic) is True
    ), "Should return True when any pathogenic score exceeds its threshold."


def test_is_pathogenic_score_false(auto_pp3bp4, var_data_non_pathogenic):
    """Test when pathogenic scores are below the thresholds."""
    assert (
        auto_pp3bp4._is_pathogenic_score(var_data_non_pathogenic) is False
    ), "Should return False when no pathogenic score exceeds its threshold."


def test_is_pathogenic_score_missing_scores(auto_pp3bp4, var_data_missing_pathogenic_scores):
    """Test when pathogenic scores are missing."""
    assert (
        auto_pp3bp4._is_pathogenic_score(var_data_missing_pathogenic_scores) is False
    ), "Should return False when pathogenic scores are missing."


def test_is_pathogenic_score_mixed(auto_pp3bp4, var_data_pathogenic, var_data_non_pathogenic):
    """Test when one pathogenic score is above the threshold and another is below."""
    var_data_pathogenic.scores.dbnsfp.bayesDel_noAF = 0.3  # Below threshold
    assert (
        auto_pp3bp4._is_pathogenic_score(var_data_pathogenic) is True
    ), "Should return True when at least one pathogenic score exceeds its threshold."


# =========== _is_benign_score ===========


@pytest.fixture
def var_data_benign():
    thresholds = MagicMock(
        metaRNN_benign=0.3,
        bayesDel_noAF_benign=0.3,
    )
    scores_dbnsfp = MagicMock(
        metaRNN=0.2,  # Below benign threshold
        bayesDel_noAF=0.1,  # Below benign threshold
    )
    scores = MagicMock(dbnsfp=scores_dbnsfp)
    return MagicMock(scores=scores, thresholds=thresholds)


@pytest.fixture
def var_data_non_benign():
    thresholds = MagicMock(
        metaRNN_benign=0.3,
        bayesDel_noAF_benign=0.3,
    )
    scores_dbnsfp = MagicMock(
        metaRNN=0.4,  # Above benign threshold
        bayesDel_noAF=0.5,  # Above benign threshold
    )
    scores = MagicMock(dbnsfp=scores_dbnsfp)
    return MagicMock(scores=scores, thresholds=thresholds)


@pytest.fixture
def var_data_missing_benign_scores():
    thresholds = MagicMock(
        metaRNN_benign=0.3,
        bayesDel_noAF_benign=0.3,
    )
    scores_dbnsfp = MagicMock(
        metaRNN=None,  # Missing score
        bayesDel_noAF=None,  # Missing score
    )
    scores = MagicMock(dbnsfp=scores_dbnsfp)
    return MagicMock(scores=scores, thresholds=thresholds)


def test_is_benign_score_true(auto_pp3bp4, var_data_benign):
    """Test when benign scores are below the thresholds."""
    assert (
        auto_pp3bp4._is_benign_score(var_data_benign) is True
    ), "Should return True when any benign score is below its threshold."


def test_is_benign_score_false(auto_pp3bp4, var_data_non_benign):
    """Test when benign scores are above the thresholds."""
    assert (
        auto_pp3bp4._is_benign_score(var_data_non_benign) is False
    ), "Should return False when no benign score is below its threshold."


def test_is_benign_score_missing_scores(auto_pp3bp4, var_data_missing_benign_scores):
    """Test when benign scores are missing."""
    assert (
        auto_pp3bp4._is_benign_score(var_data_missing_benign_scores) is False
    ), "Should return False when benign scores are missing."


def test_is_benign_score_mixed(auto_pp3bp4, var_data_benign, var_data_non_benign):
    """Test when one benign score is below the threshold and another is above."""
    var_data_benign.scores.dbnsfp.bayesDel_noAF = 0.4  # Above threshold
    assert (
        auto_pp3bp4._is_benign_score(var_data_benign) is True
    ), "Should return True when at least one benign score is below its threshold."


# =========== _is_pathogenic_splicing ===========


@pytest.fixture
def var_data_pathogenic_splicing():
    thresholds = MagicMock(
        ada=0.6,
        rf=0.7,
    )
    scores_dbscsnv = MagicMock(
        ada=0.7,  # Above ada threshold
        rf=0.8,  # Above rf threshold
    )
    scores = MagicMock(dbscsnv=scores_dbscsnv, cadd=scores_dbscsnv)
    return MagicMock(scores=scores, thresholds=thresholds)


@pytest.fixture
def var_data_non_pathogenic_splicing():
    thresholds = MagicMock(
        ada=0.6,
        rf=0.7,
    )
    scores_dbscsnv = MagicMock(
        ada=0.5,  # Below ada threshold
        rf=0.6,  # Below rf threshold
    )
    scores = MagicMock(dbscsnv=scores_dbscsnv, cadd=scores_dbscsnv)
    return MagicMock(scores=scores, thresholds=thresholds)


@pytest.fixture
def var_data_missing_pathogenic_s_scores():
    thresholds = MagicMock(
        ada=0.6,
        rf=0.7,
    )
    scores_dbscsnv = MagicMock(
        ada=None,  # Missing ada score
        rf=None,  # Missing rf score
    )
    scores = MagicMock(dbscsnv=scores_dbscsnv, cadd=scores_dbscsnv)
    return MagicMock(scores=scores, thresholds=thresholds)


def test_is_pathogenic_splicing_true(auto_pp3bp4, var_data_pathogenic_splicing):
    """Test when splicing scores are above the thresholds."""
    assert (
        auto_pp3bp4._is_pathogenic_splicing(var_data_pathogenic_splicing) is True
    ), "Should return True when any splicing score is above its threshold."


def test_is_pathogenic_splicing_false(auto_pp3bp4, var_data_non_pathogenic_splicing):
    """Test when splicing scores are below the thresholds."""
    assert (
        auto_pp3bp4._is_pathogenic_splicing(var_data_non_pathogenic_splicing) is False
    ), "Should return False when no splicing score is above its threshold."


def test_is_pathogenic_splicing_missing_scores(auto_pp3bp4, var_data_missing_pathogenic_s_scores):
    """Test when splicing scores are missing."""
    with pytest.raises(MissingDataError) as exc_info:
        auto_pp3bp4._is_pathogenic_splicing(var_data_missing_pathogenic_s_scores)
    assert "Missing Ada and RF scores" in str(
        exc_info.value
    ), "Should raise MissingDataError when both ada and rf scores are missing."


def test_is_pathogenic_splicing_mixed(auto_pp3bp4, var_data_non_pathogenic_splicing):
    """Test when one splicing score is above the threshold and another is below."""
    var_data_non_pathogenic_splicing.scores.dbscsnv.ada = 0.8  # Above threshold
    assert (
        auto_pp3bp4._is_pathogenic_splicing(var_data_non_pathogenic_splicing) is True
    ), "Should return True when at least one splicing score is above its threshold."


# =========== _is_benign_splicing ===========


@pytest.fixture
def var_data_benign_splicing():
    thresholds = MagicMock(
        ada=0.6,
        rf=0.7,
    )
    scores_dbscsnv = MagicMock(
        ada=0.5,  # Below ada threshold
        rf=0.6,  # Below rf threshold
    )
    scores = MagicMock(dbscsnv=scores_dbscsnv, cadd=scores_dbscsnv)
    return MagicMock(scores=scores, thresholds=thresholds)


@pytest.fixture
def var_data_non_benign_splicing():
    thresholds = MagicMock(
        ada=0.6,
        rf=0.7,
    )
    scores_dbscsnv = MagicMock(
        ada=0.7,  # Above ada threshold
        rf=0.8,  # Above rf threshold
    )
    scores = MagicMock(dbscsnv=scores_dbscsnv, cadd=scores_dbscsnv)
    return MagicMock(scores=scores, thresholds=thresholds)


@pytest.fixture
def var_data_missing_benign_s_scores():
    thresholds = MagicMock(
        ada=0.6,
        rf=0.7,
    )
    scores_dbscsnv = MagicMock(
        ada=None,  # Missing ada score
        rf=None,  # Missing rf score
    )
    scores = MagicMock(dbscsnv=scores_dbscsnv, cadd=scores_dbscsnv)
    return MagicMock(scores=scores, thresholds=thresholds)


def test_is_benign_splicing_true(auto_pp3bp4, var_data_benign_splicing):
    """Test when splicing scores are below the thresholds."""
    assert (
        auto_pp3bp4._is_benign_splicing(var_data_benign_splicing) is True
    ), "Should return True when any splicing score is below its threshold."


def test_is_benign_splicing_false(auto_pp3bp4, var_data_non_benign_splicing):
    """Test when splicing scores are above the thresholds."""
    assert (
        auto_pp3bp4._is_benign_splicing(var_data_non_benign_splicing) is False
    ), "Should return False when no splicing score is below its threshold."


def test_is_benign_splicing_missing_scores(auto_pp3bp4, var_data_missing_benign_s_scores):
    """Test when splicing scores are missing."""
    with pytest.raises(MissingDataError) as exc_info:
        auto_pp3bp4._is_benign_splicing(var_data_missing_benign_s_scores)
    assert "Missing Ada and RF scores" in str(
        exc_info.value
    ), "Should raise MissingDataError when both ada and rf scores are missing."


def test_is_benign_splicing_mixed(auto_pp3bp4, var_data_non_benign_splicing):
    """Test when one splicing score is below the threshold and another is above."""
    var_data_non_benign_splicing.scores.dbscsnv.ada = 0.5  # Below threshold
    assert (
        auto_pp3bp4._is_benign_splicing(var_data_non_benign_splicing) is True
    ), "Should return True when at least one splicing score is below its threshold."


# =========== verify_pp3bp4 ==============


@pytest.fixture
def seqvar_mt():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="MT", pos=100, delete="A", insert="T")


@pytest.fixture
def var_data_verify():
    thresholds = MagicMock(
        ada=0.6,
        rf=0.7,
    )
    scores_dbscsnv = MagicMock(
        ada=0.7,  # Above ada threshold
        rf=0.8,  # Above rf threshold
    )
    scores_dbnsfp = MagicMock(
        metaRNN=0.9,
        bayesDel_noAF=0.8,
    )
    consequence = MagicMock(
        cadd={"splice": True},
        mehari=["splice"],
    )
    scores = MagicMock(dbscsnv=scores_dbscsnv, cadd=scores_dbscsnv, dbnsfp=scores_dbnsfp)
    return MagicMock(scores=scores, thresholds=thresholds, consequence=consequence)


@pytest.fixture
def var_data_verify_non_splice():
    thresholds = MagicMock(
        metaRNN_pathogenic=0.85,
        bayesDel_noAF_pathogenic=0.75,
    )
    scores_dbscsnv = MagicMock(
        ada=0.4,  # Below ada threshold
        rf=0.5,  # Below rf threshold
    )
    scores_dbnsfp = MagicMock(
        metaRNN=0.9,  # Above pathogenic threshold
        bayesDel_noAF=0.8,  # Above pathogenic threshold
    )
    consequence = MagicMock(
        cadd={"missense": True},
        mehari=["missense"],
    )
    scores = MagicMock(dbscsnv=scores_dbscsnv, cadd=scores_dbscsnv, dbnsfp=scores_dbnsfp)
    return MagicMock(scores=scores, thresholds=thresholds, consequence=consequence)


@patch.object(AutoPP3BP4, "_splice_variant", return_value=True)
@patch.object(AutoPP3BP4, "_is_pathogenic_splicing", return_value=True)
@patch.object(AutoPP3BP4, "_is_benign_splicing", return_value=False)
def test_verify_pp3bp4_splice_variant(
    mock_benign_splicing,
    mock_pathogenic_splicing,
    mock_splice_variant,
    auto_pp3bp4,
    seqvar,
    var_data_verify,
):
    """Test verify_pp3bp4 when the variant is a splice variant."""
    prediction, comment = auto_pp3bp4.verify_pp3bp4(seqvar, var_data_verify)
    assert prediction.PP3 is True
    assert prediction.BP4 is False
    assert "Variant is a splice variant." in comment


@patch.object(AutoPP3BP4, "_splice_variant", return_value=False)
@patch.object(AutoPP3BP4, "_is_pathogenic_score", return_value=True)
@patch.object(AutoPP3BP4, "_is_benign_score", return_value=False)
def test_verify_pp3bp4_non_splice_variant(
    mock_benign_score,
    mock_pathogenic_score,
    mock_splice_variant,
    auto_pp3bp4,
    seqvar,
    var_data_verify_non_splice,
):
    """Test verify_pp3bp4 when the variant is not a splice variant."""
    prediction, comment = auto_pp3bp4.verify_pp3bp4(seqvar, var_data_verify_non_splice)
    assert prediction.PP3 is True
    assert prediction.BP4 is False
    assert "Variant is not a splice variant." in comment


@patch.object(AutoPP3BP4, "_splice_variant", return_value=True)
@patch.object(AutoPP3BP4, "_is_pathogenic_splicing", return_value=False)
@patch.object(AutoPP3BP4, "_is_benign_splicing", return_value=True)
def test_verify_pp3bp4_splice_variant_benign(
    mock_benign_splicing,
    mock_pathogenic_splicing,
    mock_splice_variant,
    auto_pp3bp4,
    seqvar,
    var_data_verify,
):
    """Test verify_pp3bp4 when the variant is a splice variant and benign."""
    prediction, comment = auto_pp3bp4.verify_pp3bp4(seqvar, var_data_verify)
    assert prediction.PP3 is False
    assert prediction.BP4 is True
    assert "Variant is a splice variant." in comment


@patch.object(AutoPP3BP4, "_splice_variant", return_value=False)
@patch.object(AutoPP3BP4, "_is_pathogenic_score", return_value=False)
@patch.object(AutoPP3BP4, "_is_benign_score", return_value=True)
def test_verify_pp3bp4_non_splice_variant_benign(
    mock_benign_score,
    mock_pathogenic_score,
    mock_splice_variant,
    auto_pp3bp4,
    seqvar,
    var_data_verify_non_splice,
):
    """Test verify_pp3bp4 when the variant is not a splice variant and benign."""
    prediction, comment = auto_pp3bp4.verify_pp3bp4(seqvar, var_data_verify_non_splice)
    assert prediction.PP3 is False
    assert prediction.BP4 is True
    assert "Variant is not a splice variant." in comment


@patch.object(
    AutoPP3BP4,
    "_splice_variant",
    side_effect=AutoAcmgBaseException("Error predicting splice variant"),
)
def test_verify_pp3bp4_exception(mock_splice_variant, auto_pp3bp4, seqvar, var_data_verify):
    """Test verify_pp3bp4 when an exception occurs."""
    prediction, comment = auto_pp3bp4.verify_pp3bp4(seqvar, var_data_verify)
    assert prediction is None
    assert "An error occurred during prediction." in comment


def test_verify_pp3bp4_mitochondrial(auto_pp3bp4, seqvar_mt, var_data_verify):
    """Test verify_pp3bp4 when the variant is in mitochondrial DNA."""
    prediction, comment = auto_pp3bp4.verify_pp3bp4(seqvar_mt, var_data_verify)
    assert prediction.PP3 is False
    assert prediction.BP4 is False
    assert "Variant is in mitochondrial DNA" in comment


# =========== predict_pp3bp4 ===========


@pytest.fixture
def var_data():
    thresholds = MagicMock(
        metaRNN_pathogenic=0.85,
        bayesDel_noAF_pathogenic=0.75,
    )
    scores_dbnsfp = MagicMock(
        metaRNN=0.9,  # Above pathogenic threshold
        bayesDel_noAF=0.8,  # Above pathogenic threshold
    )
    consequence = MagicMock(
        cadd={"missense": True},
        mehari=["missense"],
    )
    scores = MagicMock(dbnsfp=scores_dbnsfp)
    return MagicMock(scores=scores, thresholds=thresholds, consequence=consequence)


@pytest.fixture
def pp3bp4_result_met():
    return (
        PP3BP4(
            PP3=True,
            BP4=False,
            PP3_strength=AutoACMGStrength.PathogenicSupporting,
            BP4_strength=AutoACMGStrength.BenignSupporting,
        ),
        "PP3 criteria met",
    )


@pytest.fixture
def pp3bp4_result_not_met():
    return (
        PP3BP4(
            PP3=False,
            BP4=True,
            PP3_strength=AutoACMGStrength.PathogenicSupporting,
            BP4_strength=AutoACMGStrength.BenignSupporting,
        ),
        "BP4 criteria met",
    )


@pytest.fixture
def pp3bp4_result_failed():
    return None, "Failed to evaluate."


@patch("src.criteria.auto_pp3_bp4.AutoPP3BP4.verify_pp3bp4")
def test_predict_pp3bp4_met(mock_verify, auto_pp3bp4, seqvar, var_data, pp3bp4_result_met):
    """Test predict_pp3bp4 where PP3 criterion is met."""
    mock_verify.return_value = pp3bp4_result_met
    pp3, bp4 = auto_pp3bp4.predict_pp3bp4(seqvar, var_data)
    assert pp3.prediction == AutoACMGPrediction.Met
    assert pp3.strength == AutoACMGStrength.PathogenicSupporting
    assert "PP3 criteria met" in pp3.summary
    assert bp4.prediction == AutoACMGPrediction.NotMet
    assert bp4.strength == AutoACMGStrength.BenignSupporting
    assert "PP3 criteria met" in bp4.summary


@patch("src.criteria.auto_pp3_bp4.AutoPP3BP4.verify_pp3bp4")
def test_predict_pp3bp4_not_met(mock_verify, auto_pp3bp4, seqvar, var_data, pp3bp4_result_not_met):
    """Test predict_pp3bp4 where BP4 criterion is met."""
    mock_verify.return_value = pp3bp4_result_not_met
    pp3, bp4 = auto_pp3bp4.predict_pp3bp4(seqvar, var_data)
    assert pp3.prediction == AutoACMGPrediction.NotMet
    assert pp3.strength == AutoACMGStrength.PathogenicSupporting
    assert "BP4 criteria met" in pp3.summary
    assert bp4.prediction == AutoACMGPrediction.Met
    assert bp4.strength == AutoACMGStrength.BenignSupporting
    assert "BP4 criteria met" in bp4.summary


@patch("src.criteria.auto_pp3_bp4.AutoPP3BP4.verify_pp3bp4")
def test_predict_pp3bp4_failed(mock_verify, auto_pp3bp4, seqvar, var_data, pp3bp4_result_failed):
    """Test predict_pp3bp4 when there's a failure to evaluate the criteria."""
    mock_verify.return_value = pp3bp4_result_failed
    pp3, bp4 = auto_pp3bp4.predict_pp3bp4(seqvar, var_data)
    assert pp3.prediction == AutoACMGPrediction.Failed
    assert pp3.strength == AutoACMGStrength.PathogenicSupporting
    assert "Failed to evaluate." in pp3.summary
    assert bp4.prediction == AutoACMGPrediction.Failed
    assert bp4.strength == AutoACMGStrength.BenignSupporting
    assert "Failed to evaluate." in bp4.summary
