from unittest.mock import MagicMock, patch

import pytest

from src.criteria.auto_pm2_ba1_bs1_bs2 import AutoPM2BA1BS1BS2
from src.defs.auto_acmg import AlleleCondition, AutoACMGPrediction, AutoACMGStrength
from src.defs.exceptions import MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


@pytest.fixture
def auto_pm2ba1bs1bs2():
    return AutoPM2BA1BS1BS2()


# =========== _get_control_af ===========


@pytest.fixture
def gnomad_exomes_data():
    return MagicMock(
        alleleCounts=[
            MagicMock(cohort="controls", afGrpmax=0.01),  # Control population
            MagicMock(cohort="non_controls", afGrpmax=0.05),  # Non-control population
        ]
    )


@pytest.fixture
def gnomad_exomes_no_controls():
    return MagicMock(
        alleleCounts=[MagicMock(cohort="non_controls", afGrpmax=0.05)]  # Only non-control data
    )


@pytest.fixture
def gnomad_exomes_empty():
    return MagicMock(alleleCounts=[])  # No data


def test_get_control_af_success(auto_pm2ba1bs1bs2, gnomad_exomes_data):
    """Test successful retrieval of control allele frequency."""
    result = auto_pm2ba1bs1bs2._get_control_af(gnomad_exomes_data)
    assert result is not None
    assert result.afGrpmax == 0.01
    assert result.cohort == "controls"


def test_get_control_af_no_controls(auto_pm2ba1bs1bs2, gnomad_exomes_no_controls):
    """Test the case where no control data is available."""
    result = auto_pm2ba1bs1bs2._get_control_af(gnomad_exomes_no_controls)
    assert result is None


def test_get_control_af_empty_data(auto_pm2ba1bs1bs2, gnomad_exomes_empty):
    """Test the case where allele counts are empty."""
    result = auto_pm2ba1bs1bs2._get_control_af(gnomad_exomes_empty)
    assert result is None


# =========== _get_any_af ===========


@pytest.fixture
def gnomad_exomes_mixed_data():
    return MagicMock(
        alleleCounts=[
            MagicMock(cohort="non_controls", afGrpmax=0.02),
            MagicMock(cohort="controls", afGrpmax=0.03),
            MagicMock(cohort="non_controls", afGrpmax=0.04),
        ]
    )


@pytest.fixture
def gnomad_exomes_no_controls_higher_af():
    return MagicMock(
        alleleCounts=[
            MagicMock(cohort="non_controls", afGrpmax=0.05),
            MagicMock(cohort="non_controls", afGrpmax=0.07),
        ]
    )


def test_get_any_af_control_preference(auto_pm2ba1bs1bs2, gnomad_exomes_mixed_data):
    """Test retrieval of control data when mixed data is available."""
    result = auto_pm2ba1bs1bs2._get_any_af(gnomad_exomes_mixed_data)
    assert result is not None
    assert result.afGrpmax == 0.03
    assert result.cohort == "controls"


def test_get_any_af_highest_non_control(auto_pm2ba1bs1bs2, gnomad_exomes_no_controls_higher_af):
    """Test retrieval of the highest allele frequency from non-control data."""
    result = auto_pm2ba1bs1bs2._get_any_af(gnomad_exomes_no_controls_higher_af)
    assert result is not None
    assert result.afGrpmax == 0.07
    assert result.cohort == "non_controls"


def test_get_any_af_no_data(auto_pm2ba1bs1bs2, gnomad_exomes_empty):
    """Test behavior when no allele count data is available."""
    result = auto_pm2ba1bs1bs2._get_any_af(gnomad_exomes_empty)
    assert result is None


# =========== _get_af ==================


@pytest.fixture
def seqvar_mitochondrial():
    return MagicMock(chrom="MT", position=12345)


@pytest.fixture
def seqvar_non_mitochondrial():
    return MagicMock(chrom="1", position=67890)


@pytest.fixture
def gnomad_mtdna_data():
    return MagicMock(afHet=0.015)


@pytest.fixture
def gnomad_exomes_data_af():
    return MagicMock(
        alleleCounts=[
            MagicMock(cohort="controls", afGrpmax=0.02),
            MagicMock(cohort="non_controls", afGrpmax=0.03),
        ]
    )


@pytest.fixture
def gnomad_exomes_empty_af():
    return MagicMock(alleleCounts=[])


def test_get_af_mitochondrial_variant(
    auto_pm2ba1bs1bs2, seqvar_mitochondrial, gnomad_mtdna_data, gnomad_exomes_data_af
):
    """Test retrieval of allele frequency for mitochondrial variant."""
    result = auto_pm2ba1bs1bs2._get_af(
        seqvar_mitochondrial, gnomad_mtdna_data, gnomad_exomes_data_af
    )
    assert result == 0.015


def test_get_af_non_mitochondrial_with_controls(
    auto_pm2ba1bs1bs2, seqvar_non_mitochondrial, gnomad_exomes_data_af
):
    """Test retrieval of allele frequency for non-mitochondrial variant with controls data."""
    result = auto_pm2ba1bs1bs2._get_af(seqvar_non_mitochondrial, None, gnomad_exomes_data_af)
    assert result == 0.02


def test_get_af_non_mitochondrial_without_controls(
    auto_pm2ba1bs1bs2, seqvar_non_mitochondrial, gnomad_exomes_data_af
):
    """Test retrieval of allele frequency for non-mitochondrial variant without controls data."""
    result = auto_pm2ba1bs1bs2._get_af(seqvar_non_mitochondrial, None, gnomad_exomes_data_af)
    assert result == 0.02


def test_get_af_missing_mtdna_data(auto_pm2ba1bs1bs2, seqvar_mitochondrial, gnomad_exomes_data_af):
    """Test handling of missing mitochondrial gnomad data."""
    with pytest.raises(MissingDataError):
        auto_pm2ba1bs1bs2._get_af(seqvar_mitochondrial, None, gnomad_exomes_data_af)


def test_get_af_missing_exomes_data(auto_pm2ba1bs1bs2, seqvar_non_mitochondrial, gnomad_mtdna_data):
    """Test handling of missing gnomad exomes data."""
    with pytest.raises(MissingDataError):
        auto_pm2ba1bs1bs2._get_af(seqvar_non_mitochondrial, gnomad_mtdna_data, None)


# =========== _get_allele_cond ===========

# TODO: AnnonarsClient is not patched properly. Fix this tests!

# @pytest.fixture
# def seqvar():
#     return MagicMock(chrom="1", pos=100, delete="A", insert="T")


# @pytest.fixture
# def gene_transcript_data():
#     return MagicMock(geneId="GENE123")


# @pytest.fixture
# def gene_info_with_clingen():
#     gene_info = MagicMock()
#     gene_info.genes.root.get.return_value = {
#         "clingen": MagicMock(haploinsufficiencyScore="Dominant")
#     }
#     return gene_info


# @pytest.fixture
# def gene_info_without_clingen():
#     gene_info = MagicMock()
#     gene_info.genes.root.get.return_value = None
#     return gene_info


# @pytest.fixture
# def gene_info_with_decipher_high_pHi():
#     gene_info = MagicMock()
#     gene_info.genes.root.get.return_value = {"decipherHi": MagicMock(pHi=0.95)}
#     return gene_info


# @pytest.fixture
# def gene_info_with_domino_high_score():
#     gene_info = MagicMock()
#     gene_info.genes.root.get.return_value = {"domino": MagicMock(score=0.6)}
#     return gene_info


# @patch("src.utils.SeqVarTranscriptsHelper")
# @patch("src.api.annonars.AnnonarsClient.get_gene_info")
# def test_get_allele_cond_with_clingen(
#     mock_get_gene_info,
#     mock_transcripts_helper,
#     auto_pm2ba1bs1bs2,
#     seqvar,
#     gene_transcript_data,
#     gene_info_with_clingen,
# ):
#     """Test allele condition retrieval with Clingen data."""
#     mock_transcripts_helper.return_value.get_ts_info.return_value = (
#         None,
#         gene_transcript_data,
#         None,
#         None,
#         None,
#     )
#     mock_get_gene_info.return_value = gene_info_with_clingen
#     condition = auto_pm2ba1bs1bs2._get_allele_cond(seqvar)
#     assert condition == AlleleCondition.Dominant


# @patch("src.utils.SeqVarTranscriptsHelper")
# @patch("src.api.annonars.AnnonarsClient.get_gene_info")
# def test_get_allele_cond_with_decipher_high_pHi(
#     mock_annonars_client,
#     mock_transcripts_helper,
#     auto_pm2ba1bs1bs2,
#     seqvar,
#     gene_transcript_data,
#     gene_info_with_decipher_high_pHi,
# ):
#     """Test allele condition retrieval with high Decipher pHi."""
#     mock_transcripts_helper.return_value.get_ts_info.return_value = (
#         None,
#         gene_transcript_data,
#         None,
#         None,
#         None,
#     )
#     mock_annonars_client.get_gene_info.return_value = gene_info_with_decipher_high_pHi
#     condition = auto_pm2ba1bs1bs2._get_allele_cond(seqvar)
#     assert condition == AlleleCondition.Dominant


# @patch("src.utils.SeqVarTranscriptsHelper")
# @patch("src.api.annonars.AnnonarsClient.get_gene_info")
# def test_get_allele_cond_with_domino_high_score(
#     mock_annonars_client,
#     mock_transcripts_helper,
#     auto_pm2ba1bs1bs2,
#     seqvar,
#     gene_transcript_data,
#     gene_info_with_domino_high_score,
# ):
#     """Test allele condition retrieval with high Domino score."""
#     mock_transcripts_helper.return_value.get_ts_info.return_value = (
#         None,
#         gene_transcript_data,
#         None,
#         None,
#         None,
#     )
#     mock_annonars_client.get_gene_info.return_value = gene_info_with_domino_high_score
#     condition = auto_pm2ba1bs1bs2._get_allele_cond(seqvar)
#     assert condition == AlleleCondition.Dominant


# ========== _check_zyg ===========


@pytest.fixture
def seqvar_mitochondrial_zyg():
    return MagicMock(chrom="M", pos=100, delete="A", insert="T")


@pytest.fixture
def seqvar_x_zyg():
    return MagicMock(chrom="X", pos=100, delete="A", insert="T")


@pytest.fixture
def seqvar_autosomal_zyg():
    return MagicMock(chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def gnomad_exomes_zyg():
    # Define specific attributes directly with integer values.
    controls_af = MagicMock(ac=100, nhomalt=2, cohort="controls")
    any_af = MagicMock(ac=150, nhomalt=5)
    by_sex = MagicMock(
        xx=MagicMock(ac=100, nhomalt=1),  # Define ac and nhomalt explicitly
        xy=MagicMock(ac=80, nhomalt=0),  # Define ac and nhomalt explicitly
        overall=MagicMock(
            ac=180, nhomalt=3
        ),  # Define ac and nhomalt explicitly for autosomal example
    )
    return MagicMock(alleleCounts=[controls_af, any_af], bySex=by_sex)


@patch(
    "src.criteria.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._get_allele_cond",
    return_value=AlleleCondition.Dominant,
)
@patch("src.criteria.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._get_control_af")
@patch("src.criteria.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._get_any_af")
def test_check_zyg_mitochondrial(
    mock_get_any_af,
    mock_get_control_af,
    mock_get_allele_cond,
    auto_pm2ba1bs1bs2,
    seqvar_mitochondrial_zyg,
):
    """Test that mitochondrial variants are correctly ignored for BS2 criteria."""
    result = auto_pm2ba1bs1bs2._check_zyg(seqvar_mitochondrial_zyg, None)
    assert not result
    assert (
        "Mitochondrial variants are not considered for BS2 criteria."
        in auto_pm2ba1bs1bs2.comment_pm2ba1bs1bs2
    )


# @patch(
#     "src.criteria.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._get_allele_cond",
#     return_value=AlleleCondition.Dominant,
# )
# @patch("src.criteria.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._get_control_af")
# @patch("src.criteria.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._get_any_af")
# def test_check_zyg_x_dominant(
#     mock_get_any_af,
#     mock_get_control_af,
#     mock_get_allele_cond,
#     auto_pm2ba1bs1bs2,
#     seqvar_x_zyg,
#     gnomad_exomes_zyg,
# ):
#     """Test X-linked dominant allele condition."""
#     mock_get_control_af.return_value = gnomad_exomes_zyg.bySex.xx
#     mock_get_any_af.return_value = gnomad_exomes_zyg.bySex.xy
#     result = auto_pm2ba1bs1bs2._check_zyg(seqvar_x_zyg, gnomad_exomes_zyg)
#     assert result
#     assert (
#         "The variant is in a dominant X-linked disorder." in auto_pm2ba1bs1bs2.comment_pm2ba1bs1bs2
#     )


# @patch(
#     "src.criteria.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._get_allele_cond",
#     return_value=AlleleCondition.Recessive,
# )
# @patch("src.criteria.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._get_control_af")
# @patch("src.criteria.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._get_any_af")
# def test_check_zyg_autosomal_recessive(
#     mock_get_any_af,
#     mock_get_control_af,
#     mock_get_allele_cond,
#     auto_pm2ba1bs1bs2,
#     seqvar_autosomal_zyg,
#     gnomad_exomes_zyg,
# ):
#     """Test autosomal recessive allele condition."""
#     mock_get_control_af.return_value = gnomad_exomes_zyg.bySex.overall
#     mock_get_any_af.return_value = gnomad_exomes_zyg.bySex.overall
#     result = auto_pm2ba1bs1bs2._check_zyg(seqvar_autosomal_zyg, gnomad_exomes_zyg)
#     assert result
#     assert (
#         "The variant is in a recessive (homozygous) disorder."
#         in auto_pm2ba1bs1bs2.comment_pm2ba1bs1bs2
#     )


# =========== _ba1_exception ==================


@pytest.fixture
def seqvar_in_exception():
    return SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="chr3",
        pos=128598490,
        delete="C",
        insert="CTAAG",
        user_repr="NM_014049.4:c.-44_-41dupTAAG",
    )


# Fixture to create a sequence variant that is not in the BA1 exception list
@pytest.fixture
def seqvar_not_in_exception():
    return SeqVar(
        genome_release=GenomeRelease.GRCh38,
        chrom="chr22",
        pos=50505050,
        delete="G",
        insert="A",
        user_repr="NM_999999.1:c.999G>A",
    )


def test_ba1_exception_true(seqvar_in_exception):
    """Test that _ba1_exception returns True for a variant in the exception list."""
    assert AutoPM2BA1BS1BS2._ba1_exception(seqvar_in_exception) is True


def test_ba1_exception_false(seqvar_not_in_exception):
    """Test that _ba1_exception returns False for a variant not in the exception list."""
    assert AutoPM2BA1BS1BS2._ba1_exception(seqvar_not_in_exception) is False


# =========== verify_pm2ba1bs1bs2 ==============

# TODO: Patching doesn't work...

# @pytest.fixture
# def seqvar():
#     return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")

# @pytest.fixture
# def var_data(mock_thresholds):
#     return MagicMock(gnomad_mtdna=None, gnomad_exomes=None, thresholds=mock_thresholds)

# @pytest.fixture
# def mock_thresholds():
#     return MagicMock(ba1_benign=0.05, bs1_benign=0.01, pm2_pathogenic=0.005)

# # Patch the _get_af method correctly by providing the full path where it's being used
# @patch('src.criteria.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._get_af', return_value=None)
# def test_no_allele_frequency_data_found(mock_get_af, auto_pm2ba1bs1bs2, seqvar, var_data):
#     prediction, comment = auto_pm2ba1bs1bs2.verify_pm2ba1bs1bs2(seqvar, var_data)
#     assert "No allele frequency data found" in comment
#     assert prediction is None

# # Patching multiple methods to test interactions
# @patch('src.criteria.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._ba1_exception', return_value=True)
# @patch('src.criteria.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._get_af', return_value=0.06)
# def test_variant_in_ba1_exception_list(mock_get_af, mock_ba1_exception, auto_pm2ba1bs1bs2, seqvar, var_data):
#     prediction, comment = auto_pm2ba1bs1bs2.verify_pm2ba1bs1bs2(seqvar, var_data)
#     assert "The variant is in the exception list for BA1 criteria" in comment
#     assert not prediction.BA1 and not prediction.BS1

# # Ensure patches are applied where the methods are called and use pytest fixtures efficiently
# @pytest.mark.parametrize("af,expected", [
#     (0.06, True),
#     (0.02, False),
#     (0.01, False),
#     (0.004, False)
# ])
# @patch('src.criteria.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._get_af')
# def test_ba1_criteria(mock_get_af, af, expected, auto_pm2ba1bs1bs2, seqvar, var_data):
#     mock_get_af.return_value = af
#     prediction, comment = auto_pm2ba1bs1bs2.verify_pm2ba1bs1bs2(seqvar, var_data)
#     assert (prediction.BA1 is expected)


# =========== predict_pm2ba1bs1bs2 ==============


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def var_data():
    return MagicMock(thresholds={"ba1_benign": 0.05, "bs1_benign": 0.01, "pm2_pathogenic": 0.005})


def create_pred_object(pm2=False, ba1=False, bs1=False, bs2=False):
    pred = MagicMock()
    pred.PM2 = pm2
    pred.BA1 = ba1
    pred.BS1 = bs1
    pred.BS2 = bs2
    pred.PM2_strength = AutoACMGStrength.PathogenicModerate
    pred.BA1_strength = AutoACMGStrength.BenignStandAlone
    pred.BS1_strength = AutoACMGStrength.BenignStrong
    pred.BS2_strength = AutoACMGStrength.BenignStrong
    return pred


@patch.object(AutoPM2BA1BS1BS2, "verify_pm2ba1bs1bs2")
def test_all_criteria_met(mock_verify, auto_pm2ba1bs1bs2, seqvar, var_data):
    mock_verify.return_value = (create_pred_object(True, True, True, True), "All criteria met.")
    results = auto_pm2ba1bs1bs2.predict_pm2ba1bs1bs2(seqvar, var_data)
    assert all([r.prediction == AutoACMGPrediction.Met for r in results])
    assert all(["All criteria met." in r.summary for r in results])


@patch.object(AutoPM2BA1BS1BS2, "verify_pm2ba1bs1bs2")
def test_mixed_conditions(mock_verify, auto_pm2ba1bs1bs2, seqvar, var_data):
    mock_verify.return_value = (create_pred_object(True, False, True, False), "Mixed conditions.")
    results = auto_pm2ba1bs1bs2.predict_pm2ba1bs1bs2(seqvar, var_data)
    assert results[0].prediction == AutoACMGPrediction.Met
    assert results[1].prediction == AutoACMGPrediction.NotMet
    assert results[2].prediction == AutoACMGPrediction.Met
    assert results[3].prediction == AutoACMGPrediction.NotMet


@patch.object(AutoPM2BA1BS1BS2, "verify_pm2ba1bs1bs2")
def test_all_criteria_not_met(mock_verify, auto_pm2ba1bs1bs2, seqvar, var_data):
    mock_verify.return_value = (create_pred_object(False, False, False, False), "No criteria met.")
    results = auto_pm2ba1bs1bs2.predict_pm2ba1bs1bs2(seqvar, var_data)
    assert all([r.prediction == AutoACMGPrediction.NotMet for r in results])


@patch.object(AutoPM2BA1BS1BS2, "verify_pm2ba1bs1bs2")
def test_failure_condition(mock_verify, auto_pm2ba1bs1bs2, seqvar, var_data):
    mock_verify.return_value = (None, "Error occurred.")
    results = auto_pm2ba1bs1bs2.predict_pm2ba1bs1bs2(seqvar, var_data)
    assert all([r.prediction == AutoACMGPrediction.Failed for r in results])
    assert all(["Error occurred." in r.summary for r in results])
