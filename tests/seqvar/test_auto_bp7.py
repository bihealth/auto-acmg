from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import (
    BP7,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    GenomicStrand,
)
from src.defs.exceptions import AutoAcmgBaseException, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.auto_bp7 import AutoBP7


@pytest.fixture
def auto_bp7():
    return AutoBP7()


@pytest.fixture
def seqvar_mt():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="MT", pos=100, delete="A", insert="T")


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


# =========== _spliceai_impact ===========


@pytest.fixture
def var_data_spliceai_impact():
    thresholds = MagicMock(
        spliceAI_acceptor_gain=0.2,
        spliceAI_acceptor_loss=0.2,
        spliceAI_donor_gain=0.2,
        spliceAI_donor_loss=0.2,
    )
    scores = MagicMock(
        spliceAI_acceptor_gain=0.3,
        spliceAI_acceptor_loss=0.1,
        spliceAI_donor_gain=0.05,
        spliceAI_donor_loss=0.25,
    )
    cadd = MagicMock(cadd=scores)
    return MagicMock(scores=cadd, thresholds=thresholds)


def test_spliceai_impact_true(auto_bp7, var_data_spliceai_impact):
    """Test when spliceAI scores are above the thresholds"""
    assert (
        auto_bp7._spliceai_impact(var_data_spliceai_impact) is True
    ), "Should return True if any SpliceAI score exceeds its threshold"


def test_spliceai_impact_false(auto_bp7, var_data_spliceai_impact):
    """Test when spliceAI scores are below the thresholds"""
    # Adjusting the scores to be below the threshold for testing
    var_data_spliceai_impact.scores.cadd.spliceAI_acceptor_gain = 0.1
    var_data_spliceai_impact.scores.cadd.spliceAI_donor_loss = 0.1
    assert (
        auto_bp7._spliceai_impact(var_data_spliceai_impact) is False
    ), "Should return False if no SpliceAI score exceeds its threshold"


def test_spliceai_impact_missing_scores(auto_bp7, var_data_spliceai_impact):
    """Test when some spliceAI scores are None"""
    var_data_spliceai_impact.scores.cadd.spliceAI_acceptor_gain = None
    var_data_spliceai_impact.scores.cadd.spliceAI_donor_loss = None
    # Even with missing scores, the function should handle None values gracefully
    assert (
        auto_bp7._spliceai_impact(var_data_spliceai_impact) is False
    ), "Should handle None values and return False if no remaining score exceeds its threshold"


def test_spliceai_impact_missing_scores_1(auto_bp7, var_data_spliceai_impact):
    """Test when some spliceAI scores are None, but some are above the threshold"""
    var_data_spliceai_impact.scores.cadd.spliceAI_acceptor_gain = None
    var_data_spliceai_impact.scores.cadd.spliceAI_donor_loss = 0.4
    # Even with missing scores, the function should handle None values gracefully
    assert (
        auto_bp7._spliceai_impact(var_data_spliceai_impact) is True
    ), "Should handle None values and return True if any remaining score exceeds its threshold"


# =========== _is_conserved ===========


@pytest.fixture
def var_data_conserve():
    thresholds = MagicMock(phyloP100=1.0)
    scores_cadd = MagicMock(phyloP100=1.1)
    scores_dbnsfp = MagicMock(phyloP100=0.9)
    scores = MagicMock(cadd=scores_cadd, dbnsfp=scores_dbnsfp)
    return MagicMock(scores=scores, thresholds=thresholds)


def test_is_conserved_true(auto_bp7, var_data_conserve):
    """Test when phyloP100 score meets or exceeds the threshold."""
    assert (
        auto_bp7._is_conserved(var_data_conserve) is True
    ), "Should return True if phyloP100 meets or exceeds threshold"


def test_is_conserved_false(auto_bp7, var_data_conserve):
    """Test when phyloP100 score is below the threshold."""
    var_data_conserve.scores.cadd.phyloP100 = 0.8
    var_data_conserve.scores.dbnsfp.phyloP100 = None
    assert (
        auto_bp7._is_conserved(var_data_conserve) is False
    ), "Should return False if phyloP100 is below threshold"


def test_is_conserved_missing_scores(auto_bp7, var_data_conserve):
    """Test handling when phyloP100 score is missing."""
    var_data_conserve.scores.cadd.phyloP100 = None
    var_data_conserve.scores.dbnsfp.phyloP100 = None
    with pytest.raises(MissingDataError) as exc_info:
        auto_bp7._is_conserved(var_data_conserve)
    assert "Missing phyloP100 score" in str(
        exc_info.value
    ), "Should raise MissingDataError when phyloP100 score is missing"


def test_is_conserved_dbnsfp_priority(auto_bp7, var_data_conserve):
    """Test the prioritization of dbnsfp phyloP100 when cadd phyloP100 is missing."""
    var_data_conserve.scores.cadd.phyloP100 = None
    var_data_conserve.scores.dbnsfp.phyloP100 = 1.2
    assert (
        auto_bp7._is_conserved(var_data_conserve) is True
    ), "Should return True when dbnsfp phyloP100 is valid and meets the threshold"


def test_is_conserved_cadd_priority(auto_bp7, var_data_conserve):
    """Test the prioritization of cadd phyloP100 when dbnsfp phyloP100 is missing."""
    var_data_conserve.scores.cadd.phyloP100 = 0.5
    var_data_conserve.scores.dbnsfp.phyloP100 = None
    assert (
        auto_bp7._is_conserved(var_data_conserve) is False
    ), "Should return False when cadd phyloP100 is valid but below the threshold"


# =========== _is_synonymous ===========


@pytest.fixture
def auto_acmg_data_synonymous_mehari():
    data = AutoACMGSeqVarData()
    consequence = MagicMock()
    consequence.mehari = ["synonymous_variant"]
    consequence.cadd = "nonsynonymous"
    data.consequence = consequence
    return data


@pytest.fixture
def auto_acmg_data_synonymous_cadd():
    data = AutoACMGSeqVarData()
    consequence = MagicMock()
    consequence.mehari = ["missense_variant"]
    consequence.cadd = "synonymous"
    data.consequence = consequence
    return data


@pytest.fixture
def auto_acmg_data_nonsynonymous():
    data = AutoACMGSeqVarData()
    consequence = MagicMock()
    consequence.mehari = ["missense_variant"]
    consequence.cadd = "nonsynonymous"
    data.consequence = consequence
    return data


def test_is_synonymous_mehari(auto_acmg_data_synonymous_mehari):
    """Test when the variant is predicted as synonymous based on Mehari consequence."""
    result = AutoBP7._is_synonymous(auto_acmg_data_synonymous_mehari)
    assert (
        result is True
    ), "The variant should be identified as synonymous based on Mehari consequence."


def test_is_synonymous_cadd(auto_acmg_data_synonymous_cadd):
    """Test when the variant is predicted as synonymous based on CADD consequence."""
    result = AutoBP7._is_synonymous(auto_acmg_data_synonymous_cadd)
    assert (
        result is True
    ), "The variant should be identified as synonymous based on CADD consequence."


def test_is_nonsynonymous(auto_acmg_data_nonsynonymous):
    """Test when the variant is not synonymous based on both Mehari and CADD consequences."""
    result = AutoBP7._is_synonymous(auto_acmg_data_nonsynonymous)
    assert (
        result is False
    ), "The variant should not be identified as synonymous when neither Mehari nor CADD indicate it."


# =========== _is_intronic ===========


@pytest.fixture
def auto_acmg_data_intronic_mehari():
    data = AutoACMGSeqVarData()
    consequence = MagicMock()
    consequence.mehari = ["intron_variant"]
    consequence.cadd = "nonsynonymous"
    data.consequence = consequence
    return data


@pytest.fixture
def auto_acmg_data_intronic_cadd():
    data = AutoACMGSeqVarData()
    consequence = MagicMock()
    consequence.mehari = ["missense_variant"]
    consequence.cadd = "intron"
    data.consequence = consequence
    return data


@pytest.fixture
def auto_acmg_data_intronic_splice_cadd():
    data = AutoACMGSeqVarData()
    consequence = MagicMock()
    consequence.mehari = ["missense_variant"]
    consequence.cadd = "splice"
    data.consequence = consequence
    return data


@pytest.fixture
def auto_acmg_data_non_intronic():
    data = AutoACMGSeqVarData()
    consequence = MagicMock()
    consequence.mehari = ["missense_variant"]
    consequence.cadd = "nonsynonymous"
    data.consequence = consequence
    return data


def test_is_intronic_mehari(auto_acmg_data_intronic_mehari):
    """Test when the variant is predicted as intronic based on Mehari consequence."""
    result = AutoBP7._is_intronic(auto_acmg_data_intronic_mehari)
    assert (
        result is True
    ), "The variant should be identified as intronic based on Mehari consequence."


def test_is_intronic_cadd(auto_acmg_data_intronic_cadd):
    """Test when the variant is predicted as intronic based on CADD consequence."""
    result = AutoBP7._is_intronic(auto_acmg_data_intronic_cadd)
    assert result is True, "The variant should be identified as intronic based on CADD consequence."


def test_is_intronic_splice_cadd(auto_acmg_data_intronic_splice_cadd):
    """Test when the variant is predicted as intronic based on splice-related CADD consequence."""
    result = AutoBP7._is_intronic(auto_acmg_data_intronic_splice_cadd)
    assert (
        result is True
    ), "The variant should be identified as intronic based on CADD splice consequence."


def test_is_non_intronic(auto_acmg_data_non_intronic):
    """Test when the variant is not intronic based on both Mehari and CADD consequences."""
    result = AutoBP7._is_intronic(auto_acmg_data_non_intronic)
    assert (
        result is False
    ), "The variant should not be identified as intronic when neither Mehari nor CADD indicate it."


# =========== _affect_canonical_ss ===========


@pytest.fixture
def auto_acmg_data():
    """Fixture for an AutoACMGData object."""
    thresholds = MagicMock(
        bp7_acceptor=2, bp7_donor=2
    )  # Mock thresholds for acceptor and donor sites
    exons = [
        MagicMock(altStartI=98, altEndI=102),  # Example exon data
        MagicMock(altStartI=150, altEndI=160),
    ]
    var_data = MagicMock()
    var_data.strand = GenomicStrand.Plus
    var_data.exons = exons
    var_data.thresholds = thresholds
    return var_data


def test_affect_canonical_ss_acceptor_site(seqvar, auto_acmg_data):
    """Test when the variant affects the canonical splice site (acceptor) on the plus strand."""
    seqvar.pos = 97  # Position just within the acceptor site window
    bp7_predictor = AutoBP7()

    result = bp7_predictor._affect_canonical_ss(seqvar, auto_acmg_data)

    assert (
        result is True
    ), "The method should return True when the variant affects the acceptor site on the plus strand."


def test_affect_canonical_ss_donor_site(seqvar, auto_acmg_data):
    """Test when the variant affects the canonical splice site (donor) on the plus strand."""
    seqvar.pos = 103  # Position just within the donor site window
    bp7_predictor = AutoBP7()

    result = bp7_predictor._affect_canonical_ss(seqvar, auto_acmg_data)

    assert (
        result is True
    ), "The method should return True when the variant affects the donor site on the plus strand."


def test_affect_canonical_ss_minus_strand(seqvar, auto_acmg_data):
    """Test when the variant affects the canonical splice site on the minus strand."""
    auto_acmg_data.strand = GenomicStrand.Minus
    seqvar.pos = 97  # Position within the donor site window on the minus strand
    bp7_predictor = AutoBP7()

    result = bp7_predictor._affect_canonical_ss(seqvar, auto_acmg_data)

    assert (
        result is True
    ), "The method should return True when the variant affects the canonical splice site on the minus strand."


def test_not_affect_canonical_ss(seqvar, auto_acmg_data):
    """Test when the variant does not affect any canonical splice site."""
    seqvar.pos = 110  # Position outside the canonical splice site window
    bp7_predictor = AutoBP7()

    result = bp7_predictor._affect_canonical_ss(seqvar, auto_acmg_data)

    assert (
        result is False
    ), "The method should return False when the variant does not affect any canonical splice site."


def test_missing_strand_info(seqvar, auto_acmg_data):
    """Test when the strand information is missing."""
    auto_acmg_data.strand = None  # Missing strand information
    bp7_predictor = AutoBP7()

    with pytest.raises(MissingDataError):
        bp7_predictor._affect_canonical_ss(seqvar, auto_acmg_data)


# =========== _is_bp7_exception ===========


def test_is_bp7_exception(seqvar):
    """Test that the method correctly identifies no exceptions by default."""
    bp7_predictor = AutoBP7()
    auto_acmg_data = AutoACMGSeqVarData()
    result = bp7_predictor._is_bp7_exception(seqvar, auto_acmg_data)
    assert (
        result is False
    ), "The method should return False, indicating no exceptions for BP7 by default."


# =========== verify_bp7 ===========


@pytest.fixture
def var_data_verify():
    thresholds = MagicMock(phyloP100=1.0)
    scores_cadd = MagicMock(
        phyloP100=0.8,
        spliceAI_acceptor_gain=0.1,
        spliceAI_acceptor_loss=0.1,
        spliceAI_donor_gain=0.1,
        spliceAI_donor_loss=0.1,
    )
    cadd = MagicMock(cadd=scores_cadd)
    consequence = MagicMock(mehari=["synonymous_variant"], cadd="synonymous")
    return MagicMock(scores=cadd, thresholds=thresholds, consequence=consequence)


def test_verify_bp7_mitochondrial(auto_bp7, seqvar_mt, var_data_verify):
    """Test the handling of mitochondrial variants which should not meet BP7."""
    prediction, comment = auto_bp7.verify_bp7(seqvar_mt, var_data_verify)
    assert prediction.BP7 is False
    assert "mitochondrial genome" in comment


@patch.object(AutoBP7, "_is_conserved", return_value=False)
@patch.object(AutoBP7, "_spliceai_impact", return_value=False)
def test_verify_bp7_not_conserved_not_splicing(
    mock_conserved, mock_splicing, auto_bp7, seqvar, var_data_verify
):
    """Test non-conserved and non-splicing impacting variants which should meet BP7."""
    prediction, comment = auto_bp7.verify_bp7(seqvar, var_data_verify)
    assert prediction.BP7 is True
    assert "BP7 is met" in comment


@patch.object(AutoBP7, "_is_conserved", return_value=True)
@patch.object(AutoBP7, "_spliceai_impact", return_value=True)
def test_verify_bp7_conserved_splicing(
    mock_conserved, mock_splicing, auto_bp7, seqvar, var_data_verify
):
    """Test conserved or splicing impacting variants which should not meet BP7."""
    prediction, comment = auto_bp7.verify_bp7(seqvar, var_data_verify)
    assert prediction.BP7 is False
    assert "BP7 is not met" in comment


@patch.object(
    AutoBP7,
    "_is_conserved",
    side_effect=AutoAcmgBaseException("Error calculating conservation"),
)
def test_verify_bp7_exception(mock_conserved, auto_bp7, seqvar, var_data_verify):
    """Test error handling when an exception occurs during the prediction process."""
    prediction, comment = auto_bp7.verify_bp7(seqvar, var_data_verify)
    assert prediction is None
    assert "Failed to predict BP7 criterion" in comment


# =========== predict_bp7 ===========


@pytest.fixture
def bp7_result_met():
    return BP7(BP7=True, BP7_strength=AutoACMGStrength.BenignSupporting), "Criterion met."


@pytest.fixture
def bp7_result_not_met():
    return BP7(BP7=False, BP7_strength=AutoACMGStrength.BenignSupporting), "Criterion not met."


@pytest.fixture
def bp7_result_failed():
    return None, "Failed to evaluate."


@patch("src.seqvar.auto_bp7.AutoBP7.verify_bp7")
def test_predict_bp7_met(mock_verify, auto_bp7, seqvar, bp7_result_met):
    """Test predict_bp7 where BP7 criterion is met."""
    mock_verify.return_value = bp7_result_met
    result = auto_bp7.predict_bp7(seqvar, MagicMock())
    assert result.prediction == AutoACMGPrediction.Applicable
    assert result.strength == AutoACMGStrength.BenignSupporting
    assert "Criterion met." in result.summary


@patch("src.seqvar.auto_bp7.AutoBP7.verify_bp7")
def test_predict_bp7_not_met(mock_verify, auto_bp7, seqvar, bp7_result_not_met):
    """Test predict_bp7 where BP7 criterion is not met."""
    mock_verify.return_value = bp7_result_not_met
    result = auto_bp7.predict_bp7(seqvar, MagicMock())
    assert result.prediction == AutoACMGPrediction.NotApplicable
    assert result.strength == AutoACMGStrength.BenignSupporting
    assert "Criterion not met." in result.summary


@patch("src.seqvar.auto_bp7.AutoBP7.verify_bp7")
def test_predict_bp7_failed(mock_verify, auto_bp7, seqvar, bp7_result_failed):
    """Test predict_bp7 when there's a failure to evaluate the criterion."""
    mock_verify.return_value = bp7_result_failed
    result = auto_bp7.predict_bp7(seqvar, MagicMock())
    assert result.prediction == AutoACMGPrediction.Failed
    assert result.strength == AutoACMGStrength.BenignSupporting
    assert "Failed to evaluate." in result.summary
