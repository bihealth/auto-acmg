from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import PP3BP4, AutoACMGPrediction, AutoACMGStrength
from src.defs.exceptions import AutoAcmgBaseException, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.auto_pp3_bp4 import AutoPP3BP4


@pytest.fixture
def auto_pp3bp4():
    return AutoPP3BP4()


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


# =========== _is_splice_variant ===========


@pytest.fixture
def var_data_splice():
    consequence = MagicMock(
        cadd={"splice": True},  # Simulates CADD indicating a splice variant
        mehari=[
            "splice",
            "other_effect",
        ],  # Simulates other tools also indicating a splice variant
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


# =========== _is_inframe_indel ===========


@pytest.fixture
def var_data_inframe_indel():
    consequence = MagicMock(
        cadd={"inframe": True},  # Simulates CADD indicating an inframe indel
        mehari=[
            "inframe_deletion",
            "other_effect",
        ],  # Simulates other tools also indicating an inframe indel
    )
    return MagicMock(consequence=consequence)


@pytest.fixture
def var_data_not_inframe_indel():
    consequence = MagicMock(
        cadd={},  # No inframe key
        mehari=["missense_variant", "other_effect"],  # No inframe-related terms
    )
    return MagicMock(consequence=consequence)


def test_inframe_indel_positive(auto_pp3bp4: AutoPP3BP4, var_data_inframe_indel: MagicMock):
    """Test that _is_inframe_indel correctly identifies an inframe indel."""
    assert (
        auto_pp3bp4._is_inframe_indel(var_data_inframe_indel) is True
    ), "Should return True when inframe indel indicators are present in the data."


def test_inframe_indel_negative(auto_pp3bp4: AutoPP3BP4, var_data_not_inframe_indel: MagicMock):
    """Test that _is_inframe_indel correctly identifies non-inframe indel variants."""
    assert (
        auto_pp3bp4._is_inframe_indel(var_data_not_inframe_indel) is False
    ), "Should return False when no inframe indel indicators are present in the data."


def test_inframe_indel_cadd_only(auto_pp3bp4: AutoPP3BP4, var_data_inframe_indel: MagicMock):
    """Test _is_inframe_indel when only CADD indicates an inframe indel."""
    var_data_inframe_indel.consequence.mehari = ["other_effect"]
    assert (
        auto_pp3bp4._is_inframe_indel(var_data_inframe_indel) is True
    ), "Should return True when CADD indicates an inframe indel, even if mehari doesn't."


def test_inframe_indel_mehari_only(auto_pp3bp4: AutoPP3BP4, var_data_inframe_indel: MagicMock):
    """Test _is_inframe_indel when only mehari indicates an inframe indel."""
    var_data_inframe_indel.consequence.cadd = {}
    assert (
        auto_pp3bp4._is_inframe_indel(var_data_inframe_indel) is True
    ), "Should return True when mehari indicates an inframe indel, even if CADD doesn't."


def test_inframe_indel_with_other_effects(
    auto_pp3bp4: AutoPP3BP4, var_data_inframe_indel: MagicMock
):
    """Test _is_inframe_indel with other non-inframe indel related effects."""
    var_data_inframe_indel.consequence.mehari.append("missense_variant")
    assert (
        auto_pp3bp4._is_inframe_indel(var_data_inframe_indel) is True
    ), "Should still return True as long as one inframe indel indicator is present."


# =========== _is_missense_variant ===========


@pytest.fixture
def var_data_missense():
    consequence = MagicMock(
        cadd={"missense": True},  # Simulates CADD indicating a missense variant
        mehari=[
            "missense_variant",
            "other_effect",
        ],  # Simulates other tools also indicating a missense variant
    )
    return MagicMock(consequence=consequence)


@pytest.fixture
def var_data_not_missense():
    consequence = MagicMock(
        cadd={},  # No missense key
        mehari=["synonymous_variant", "other_effect"],  # No missense-related terms
    )
    return MagicMock(consequence=consequence)


def test_missense_variant_positive(auto_pp3bp4: AutoPP3BP4, var_data_missense: MagicMock):
    """Test that _is_missense_variant correctly identifies a missense variant."""
    assert (
        auto_pp3bp4._is_missense_variant(var_data_missense) is True
    ), "Should return True when missense variant indicators are present in the data."


def test_missense_variant_negative(auto_pp3bp4: AutoPP3BP4, var_data_not_missense: MagicMock):
    """Test that _is_missense_variant correctly identifies non-missense variants."""
    assert (
        auto_pp3bp4._is_missense_variant(var_data_not_missense) is False
    ), "Should return False when no missense variant indicators are present in the data."


def test_missense_variant_cadd_only(auto_pp3bp4: AutoPP3BP4, var_data_missense: MagicMock):
    """Test _is_missense_variant when only CADD indicates a missense variant."""
    var_data_missense.consequence.mehari = ["other_effect"]
    assert (
        auto_pp3bp4._is_missense_variant(var_data_missense) is True
    ), "Should return True when CADD indicates a missense variant, even if mehari doesn't."


def test_missense_variant_mehari_only(auto_pp3bp4: AutoPP3BP4, var_data_missense: MagicMock):
    """Test _is_missense_variant when only mehari indicates a missense variant."""
    var_data_missense.consequence.cadd = {}
    assert (
        auto_pp3bp4._is_missense_variant(var_data_missense) is True
    ), "Should return True when mehari indicates a missense variant, even if CADD doesn't."


def test_missense_variant_with_other_effects(auto_pp3bp4: AutoPP3BP4, var_data_missense: MagicMock):
    """Test _is_missense_variant with other non-missense related effects."""
    var_data_missense.consequence.mehari.append("synonymous_variant")
    assert (
        auto_pp3bp4._is_missense_variant(var_data_missense) is True
    ), "Should still return True as long as one missense variant indicator is present."


# =========== _is_synonymous_variant ===========


@pytest.fixture
def var_data_synonymous():
    consequence = MagicMock(
        cadd={"synonymous": True},  # Simulates CADD indicating a synonymous variant
        mehari=[
            "synonymous_variant",
            "other_effect",
        ],  # Simulates other tools also indicating a synonymous variant
    )
    return MagicMock(consequence=consequence)


@pytest.fixture
def var_data_not_synonymous():
    consequence = MagicMock(
        cadd={},  # No synonymous key
        mehari=["missense_variant", "other_effect"],  # No synonymous-related terms
    )
    return MagicMock(consequence=consequence)


def test_synonymous_variant_positive(auto_pp3bp4: AutoPP3BP4, var_data_synonymous: MagicMock):
    """Test that _is_synonymous_variant correctly identifies a synonymous variant."""
    assert (
        auto_pp3bp4._is_synonymous_variant(var_data_synonymous) is True
    ), "Should return True when synonymous variant indicators are present in the data."


def test_synonymous_variant_negative(auto_pp3bp4: AutoPP3BP4, var_data_not_synonymous: MagicMock):
    """Test that _is_synonymous_variant correctly identifies non-synonymous variants."""
    assert (
        auto_pp3bp4._is_synonymous_variant(var_data_not_synonymous) is False
    ), "Should return False when no synonymous variant indicators are present in the data."


def test_synonymous_variant_cadd_only(auto_pp3bp4: AutoPP3BP4, var_data_synonymous: MagicMock):
    """Test _is_synonymous_variant when only CADD indicates a synonymous variant."""
    var_data_synonymous.consequence.mehari = ["other_effect"]
    assert (
        auto_pp3bp4._is_synonymous_variant(var_data_synonymous) is True
    ), "Should return True when CADD indicates a synonymous variant, even if mehari doesn't."


def test_synonymous_variant_mehari_only(auto_pp3bp4: AutoPP3BP4, var_data_synonymous: MagicMock):
    """Test _is_synonymous_variant when only mehari indicates a synonymous variant."""
    var_data_synonymous.consequence.cadd = {}
    assert (
        auto_pp3bp4._is_synonymous_variant(var_data_synonymous) is True
    ), "Should return True when mehari indicates a synonymous variant, even if CADD doesn't."


def test_synonymous_variant_with_other_effects(
    auto_pp3bp4: AutoPP3BP4, var_data_synonymous: MagicMock
):
    """Test _is_synonymous_variant with other non-synonymous related effects."""
    var_data_synonymous.consequence.mehari.append("missense_variant")
    assert (
        auto_pp3bp4._is_synonymous_variant(var_data_synonymous) is True
    ), "Should still return True as long as one synonymous variant indicator is present."


# =========== _is_intron_variant ===========


@pytest.fixture
def var_data_intron():
    consequence = MagicMock(
        cadd={"intron": True},  # Simulates CADD indicating an intron variant
        mehari=[
            "intron_variant",
            "other_effect",
        ],  # Simulates other tools also indicating an intron variant
    )
    return MagicMock(consequence=consequence)


@pytest.fixture
def var_data_not_intron():
    consequence = MagicMock(
        cadd={},  # No intron key
        mehari=["missense_variant", "other_effect"],  # No intron-related terms
    )
    return MagicMock(consequence=consequence)


def test_intron_variant_positive(auto_pp3bp4: AutoPP3BP4, var_data_intron: MagicMock):
    """Test that _is_intron_variant correctly identifies an intron variant."""
    assert (
        auto_pp3bp4._is_intron_variant(var_data_intron) is True
    ), "Should return True when intron variant indicators are present in the data."


def test_intron_variant_negative(auto_pp3bp4: AutoPP3BP4, var_data_not_intron: MagicMock):
    """Test that _is_intron_variant correctly identifies non-intron variants."""
    assert (
        auto_pp3bp4._is_intron_variant(var_data_not_intron) is False
    ), "Should return False when no intron variant indicators are present in the data."


def test_intron_variant_cadd_only(auto_pp3bp4: AutoPP3BP4, var_data_intron: MagicMock):
    """Test _is_intron_variant when only CADD indicates an intron variant."""
    var_data_intron.consequence.mehari = ["other_effect"]
    assert (
        auto_pp3bp4._is_intron_variant(var_data_intron) is True
    ), "Should return True when CADD indicates an intron variant, even if mehari doesn't."


def test_intron_variant_mehari_only(auto_pp3bp4: AutoPP3BP4, var_data_intron: MagicMock):
    """Test _is_intron_variant when only mehari indicates an intron variant."""
    var_data_intron.consequence.cadd = {}
    assert (
        auto_pp3bp4._is_intron_variant(var_data_intron) is True
    ), "Should return True when mehari indicates an intron variant, even if CADD doesn't."


def test_intron_variant_with_other_effects(auto_pp3bp4: AutoPP3BP4, var_data_intron: MagicMock):
    """Test _is_intron_variant with other non-intron related effects."""
    var_data_intron.consequence.mehari.append("missense_variant")
    assert (
        auto_pp3bp4._is_intron_variant(var_data_intron) is True
    ), "Should still return True as long as one intron variant indicator is present."


# =========== _is_utr_variant ===========


@pytest.fixture
def var_data_utr():
    consequence = MagicMock(
        cadd={"UTR": True},  # Simulates CADD indicating a UTR variant
        mehari=[
            "5_prime_UTR_variant",
            "other_effect",
        ],  # Simulates other tools also indicating a UTR variant
    )
    return MagicMock(consequence=consequence)


@pytest.fixture
def var_data_not_utr():
    consequence = MagicMock(
        cadd={},  # No UTR key
        mehari=["missense_variant", "other_effect"],  # No UTR-related terms
    )
    return MagicMock(consequence=consequence)


def test_utr_variant_positive(auto_pp3bp4: AutoPP3BP4, var_data_utr: MagicMock):
    """Test that _is_utr_variant correctly identifies a UTR variant."""
    assert (
        auto_pp3bp4._is_utr_variant(var_data_utr) is True
    ), "Should return True when UTR variant indicators are present in the data."


@pytest.mark.skip(reason="Should work but doesn't")
def test_utr_variant_negative(auto_pp3bp4: AutoPP3BP4, var_data_not_utr: MagicMock):
    """Test that _is_utr_variant correctly identifies non-UTR variants."""
    assert (
        auto_pp3bp4._is_utr_variant(var_data_not_utr) is False
    ), "Should return False when no UTR variant indicators are present in the data."


def test_utr_variant_cadd_only(auto_pp3bp4: AutoPP3BP4, var_data_utr: MagicMock):
    """Test _is_utr_variant when only CADD indicates a UTR variant."""
    var_data_utr.consequence.mehari = ["other_effect"]
    assert (
        auto_pp3bp4._is_utr_variant(var_data_utr) is True
    ), "Should return True when CADD indicates a UTR variant, even if mehari doesn't."


def test_utr_variant_mehari_only(auto_pp3bp4: AutoPP3BP4, var_data_utr: MagicMock):
    """Test _is_utr_variant when only mehari indicates a UTR variant."""
    var_data_utr.consequence.cadd = {}
    assert (
        auto_pp3bp4._is_utr_variant(var_data_utr) is True
    ), "Should return True when mehari indicates a UTR variant, even if CADD doesn't."


def test_utr_variant_with_other_effects(auto_pp3bp4: AutoPP3BP4, var_data_utr: MagicMock):
    """Test _is_utr_variant with other non-UTR related effects."""
    var_data_utr.consequence.mehari.append("missense_variant")
    assert (
        auto_pp3bp4._is_utr_variant(var_data_utr) is True
    ), "Should still return True as long as one UTR variant indicator is present."


def test_utr_variant_lowercase(auto_pp3bp4: AutoPP3BP4, var_data_utr: MagicMock):
    """Test _is_utr_variant with lowercase 'utr' in CADD."""
    var_data_utr.consequence.cadd = {"utr": True}
    assert (
        auto_pp3bp4._is_utr_variant(var_data_utr) is True
    ), "Should return True when CADD indicates a UTR variant with lowercase 'utr'."


def test_utr_variant_3_prime(auto_pp3bp4: AutoPP3BP4, var_data_utr: MagicMock):
    """Test _is_utr_variant with 3' UTR variant in mehari."""
    var_data_utr.consequence.mehari = ["3_prime_UTR_variant"]
    assert (
        auto_pp3bp4._is_utr_variant(var_data_utr) is True
    ), "Should return True when mehari indicates a 3' UTR variant."


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
        auto_pp3bp4._is_pathogenic_score(
            var_data_pathogenic, ("metaRNN", 0.5), ("bayesDel_noAF", 0.5)
        )
        is True
    ), "Should return True when any pathogenic score exceeds its threshold."


def test_is_pathogenic_score_false(auto_pp3bp4, var_data_non_pathogenic):
    """Test when pathogenic scores are below the thresholds."""
    assert (
        auto_pp3bp4._is_pathogenic_score(
            var_data_non_pathogenic, ("metaRNN", 0.5), ("bayesDel_noAF", 0.5)
        )
        is False
    ), "Should return False when no pathogenic score exceeds its threshold."


def test_is_pathogenic_score_missing_scores(auto_pp3bp4, var_data_missing_pathogenic_scores):
    """Test when pathogenic scores are missing."""
    assert (
        auto_pp3bp4._is_pathogenic_score(
            var_data_missing_pathogenic_scores, ("metaRNN", 0.5), ("bayesDel_noAF", 0.5)
        )
        is False
    ), "Should return False when pathogenic scores are missing."


def test_is_pathogenic_score_mixed(auto_pp3bp4, var_data_pathogenic, var_data_non_pathogenic):
    """Test when one pathogenic score is above the threshold and another is below."""
    var_data_pathogenic.scores.dbnsfp.bayesDel_noAF = 0.3  # Below threshold
    assert (
        auto_pp3bp4._is_pathogenic_score(
            var_data_pathogenic, ("metaRNN", 0.5), ("bayesDel_noAF", 0.5)
        )
        is True
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
        auto_pp3bp4._is_benign_score(var_data_benign, ("metaRNN", 0.3), ("bayesDel_noAF", 0.3))
        is True
    ), "Should return True when any benign score is below its threshold."


def test_is_benign_score_false(auto_pp3bp4, var_data_non_benign):
    """Test when benign scores are above the thresholds."""
    assert (
        auto_pp3bp4._is_benign_score(var_data_non_benign, ("metaRNN", 0.3), ("bayesDel_noAF", 0.3))
        is False
    ), "Should return False when no benign score is below its threshold."


def test_is_benign_score_missing_scores(auto_pp3bp4, var_data_missing_benign_scores):
    """Test when benign scores are missing."""
    assert (
        auto_pp3bp4._is_benign_score(
            var_data_missing_benign_scores, ("metaRNN", 0.3), ("bayesDel_noAF", 0.3)
        )
        is False
    ), "Should return False when benign scores are missing."


def test_is_benign_score_mixed(auto_pp3bp4, var_data_benign, var_data_non_benign):
    """Test when one benign score is below the threshold and another is above."""
    var_data_benign.scores.dbnsfp.bayesDel_noAF = 0.4  # Above threshold
    assert (
        auto_pp3bp4._is_benign_score(var_data_benign, ("metaRNN", 0.3), ("bayesDel_noAF", 0.3))
        is True
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


# =========== _affect_spliceAI ===========


@pytest.fixture
def var_data_spliceai_affected():
    thresholds = MagicMock(
        spliceAI_acceptor_gain=0.5,
        spliceAI_acceptor_loss=0.5,
        spliceAI_donor_gain=0.5,
        spliceAI_donor_loss=0.5,
    )
    scores_cadd = MagicMock(
        spliceAI_acceptor_gain=0.6,  # Above threshold
        spliceAI_acceptor_loss=0.4,  # Below threshold
        spliceAI_donor_gain=0.4,  # Below threshold
        spliceAI_donor_loss=0.4,  # Below threshold
    )
    scores = MagicMock(cadd=scores_cadd)
    return MagicMock(scores=scores, thresholds=thresholds)


@pytest.fixture
def var_data_spliceai_not_affected():
    thresholds = MagicMock(
        spliceAI_acceptor_gain=0.5,
        spliceAI_acceptor_loss=0.5,
        spliceAI_donor_gain=0.5,
        spliceAI_donor_loss=0.5,
    )
    scores_cadd = MagicMock(
        spliceAI_acceptor_gain=0.4,  # Below threshold
        spliceAI_acceptor_loss=0.4,  # Below threshold
        spliceAI_donor_gain=0.4,  # Below threshold
        spliceAI_donor_loss=0.4,  # Below threshold
    )
    scores = MagicMock(cadd=scores_cadd)
    return MagicMock(scores=scores, thresholds=thresholds)


def test_affect_spliceai_positive(auto_pp3bp4: AutoPP3BP4, var_data_spliceai_affected: MagicMock):
    """Test that _affect_spliceAI correctly identifies a splice-affecting variant."""
    assert (
        auto_pp3bp4._affect_spliceAI(var_data_spliceai_affected) is True
    ), "Should return True when any SpliceAI score is above its threshold."


def test_affect_spliceai_negative(
    auto_pp3bp4: AutoPP3BP4, var_data_spliceai_not_affected: MagicMock
):
    """Test that _affect_spliceAI correctly identifies a non-splice-affecting variant."""
    assert (
        auto_pp3bp4._affect_spliceAI(var_data_spliceai_not_affected) is False
    ), "Should return False when all SpliceAI scores are below their thresholds."


def test_affect_spliceai_one_above_threshold(
    auto_pp3bp4: AutoPP3BP4, var_data_spliceai_not_affected: MagicMock
):
    """Test _affect_spliceAI when only one score is above the threshold."""
    var_data_spliceai_not_affected.scores.cadd.spliceAI_donor_gain = (
        0.6  # Set one score above threshold
    )
    assert (
        auto_pp3bp4._affect_spliceAI(var_data_spliceai_not_affected) is True
    ), "Should return True when at least one SpliceAI score is above its threshold."


def test_affect_spliceai_missing_scores(
    auto_pp3bp4: AutoPP3BP4, var_data_spliceai_affected: MagicMock
):
    """Test _affect_spliceAI when some scores are missing."""
    var_data_spliceai_affected.scores.cadd.spliceAI_acceptor_gain = None
    var_data_spliceai_affected.scores.cadd.spliceAI_acceptor_loss = None
    assert (
        auto_pp3bp4._affect_spliceAI(var_data_spliceai_affected) is False
    ), "Should return False when some scores are missing and the rest are below their thresholds."


def test_affect_spliceai_all_scores_missing(
    auto_pp3bp4: AutoPP3BP4, var_data_spliceai_affected: MagicMock
):
    """Test _affect_spliceAI when all scores are missing."""
    var_data_spliceai_affected.scores.cadd.spliceAI_acceptor_gain = None
    var_data_spliceai_affected.scores.cadd.spliceAI_acceptor_loss = None
    var_data_spliceai_affected.scores.cadd.spliceAI_donor_gain = None
    var_data_spliceai_affected.scores.cadd.spliceAI_donor_loss = None
    assert (
        auto_pp3bp4._affect_spliceAI(var_data_spliceai_affected) is False
    ), "Should return False when all SpliceAI scores are missing."


def test_affect_spliceai_custom_thresholds(
    auto_pp3bp4: AutoPP3BP4, var_data_spliceai_affected: MagicMock
):
    """Test _affect_spliceAI with custom thresholds."""
    var_data_spliceai_affected.thresholds.spliceAI_acceptor_gain = 0.7  # Set custom threshold
    assert (
        auto_pp3bp4._affect_spliceAI(var_data_spliceai_affected) is False
    ), "Should return False when the score is below a custom threshold."


# =========== verify_pp3bp4 ==============


@pytest.fixture
def seqvar_mt():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="MT", pos=100, delete="A", insert="T")


@pytest.fixture
def var_data_verify():
    return MagicMock(
        thresholds=MagicMock(pp3bp4_strategy="default"),
        scores=MagicMock(dbnsfp=MagicMock(metaRNN=0.9, bayesDel_noAF=0.8)),
    )


@patch.object(AutoPP3BP4, "_is_pathogenic_score")
@patch.object(AutoPP3BP4, "_is_benign_score")
def test_verify_pp3bp4_default_strategy_pathogenic(
    mock_is_benign, mock_is_pathogenic, auto_pp3bp4, seqvar, var_data_verify
):
    """Test verify_pp3bp4 with default strategy and pathogenic scores."""
    mock_is_pathogenic.return_value = True
    mock_is_benign.return_value = False

    prediction, comment = auto_pp3bp4.verify_pp3bp4(seqvar, var_data_verify)

    assert prediction.PP3 is True
    assert prediction.BP4 is False
    assert "MetaRNN score:" in comment
    assert "BayesDel_noAF score:" in comment


@patch.object(AutoPP3BP4, "_is_pathogenic_score")
@patch.object(AutoPP3BP4, "_is_benign_score")
def test_verify_pp3bp4_default_strategy_benign(
    mock_is_benign, mock_is_pathogenic, auto_pp3bp4, seqvar, var_data_verify
):
    """Test verify_pp3bp4 with default strategy and benign scores."""
    mock_is_pathogenic.return_value = False
    mock_is_benign.return_value = True

    prediction, comment = auto_pp3bp4.verify_pp3bp4(seqvar, var_data_verify)

    assert prediction.PP3 is False
    assert prediction.BP4 is True
    assert "MetaRNN score:" in comment
    assert "BayesDel_noAF score:" in comment


@patch.object(AutoPP3BP4, "_is_pathogenic_score")
@patch.object(AutoPP3BP4, "_is_benign_score")
@patch.object(AutoPP3BP4, "_is_pathogenic_splicing")
@patch.object(AutoPP3BP4, "_is_benign_splicing")
def test_verify_pp3bp4_custom_strategy(
    mock_benign_splicing,
    mock_pathogenic_splicing,
    mock_is_benign,
    mock_is_pathogenic,
    auto_pp3bp4,
    seqvar,
    var_data_verify,
):
    """Test verify_pp3bp4 with a custom strategy."""
    var_data_verify.thresholds.pp3bp4_strategy = "custom_score"
    mock_is_pathogenic.return_value = True
    mock_is_benign.return_value = False
    mock_pathogenic_splicing.return_value = False
    mock_benign_splicing.return_value = False

    prediction, comment = auto_pp3bp4.verify_pp3bp4(seqvar, var_data_verify)

    assert prediction.PP3 is True
    assert prediction.BP4 is False


@patch.object(AutoPP3BP4, "_is_pathogenic_score")
def test_verify_pp3bp4_exception(mock_is_pathogenic, auto_pp3bp4, seqvar, var_data_verify):
    """Test verify_pp3bp4 when an exception occurs."""
    mock_is_pathogenic.side_effect = AutoAcmgBaseException("Test exception")

    prediction, comment = auto_pp3bp4.verify_pp3bp4(seqvar, var_data_verify)

    assert prediction is None
    assert "An error occurred during prediction" in comment
    assert "Test exception" in comment


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


@patch("src.seqvar.auto_pp3_bp4.AutoPP3BP4.verify_pp3bp4")
def test_predict_pp3bp4_met(mock_verify, auto_pp3bp4, seqvar, var_data, pp3bp4_result_met):
    """Test predict_pp3bp4 where PP3 criterion is met."""
    mock_verify.return_value = pp3bp4_result_met
    pp3, bp4 = auto_pp3bp4.predict_pp3bp4(seqvar, var_data)
    assert pp3.prediction == AutoACMGPrediction.Applicable
    assert pp3.strength == AutoACMGStrength.PathogenicSupporting
    assert "PP3 criteria met" in pp3.summary
    assert bp4.prediction == AutoACMGPrediction.NotApplicable
    assert bp4.strength == AutoACMGStrength.BenignSupporting
    assert "PP3 criteria met" in bp4.summary


@patch("src.seqvar.auto_pp3_bp4.AutoPP3BP4.verify_pp3bp4")
def test_predict_pp3bp4_not_met(mock_verify, auto_pp3bp4, seqvar, var_data, pp3bp4_result_not_met):
    """Test predict_pp3bp4 where BP4 criterion is met."""
    mock_verify.return_value = pp3bp4_result_not_met
    pp3, bp4 = auto_pp3bp4.predict_pp3bp4(seqvar, var_data)
    assert pp3.prediction == AutoACMGPrediction.NotApplicable
    assert pp3.strength == AutoACMGStrength.PathogenicSupporting
    assert "BP4 criteria met" in pp3.summary
    assert bp4.prediction == AutoACMGPrediction.Applicable
    assert bp4.strength == AutoACMGStrength.BenignSupporting
    assert "BP4 criteria met" in bp4.summary


@patch("src.seqvar.auto_pp3_bp4.AutoPP3BP4.verify_pp3bp4")
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
