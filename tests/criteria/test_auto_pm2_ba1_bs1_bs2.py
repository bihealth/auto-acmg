from unittest.mock import MagicMock, patch

import pytest

from src.criteria.auto_pm2_ba1_bs1_bs2 import AutoPM2BA1BS1BS2
from src.defs.auto_acmg import AlleleCondition
from src.defs.exceptions import MissingDataError


@pytest.fixture
def auto_pm2ba1bs1bs2_instance():
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


def test_get_control_af_success(auto_pm2ba1bs1bs2_instance, gnomad_exomes_data):
    """Test successful retrieval of control allele frequency."""
    result = auto_pm2ba1bs1bs2_instance._get_control_af(gnomad_exomes_data)
    assert result is not None
    assert result.afGrpmax == 0.01
    assert result.cohort == "controls"


def test_get_control_af_no_controls(auto_pm2ba1bs1bs2_instance, gnomad_exomes_no_controls):
    """Test the case where no control data is available."""
    result = auto_pm2ba1bs1bs2_instance._get_control_af(gnomad_exomes_no_controls)
    assert result is None


def test_get_control_af_empty_data(auto_pm2ba1bs1bs2_instance, gnomad_exomes_empty):
    """Test the case where allele counts are empty."""
    result = auto_pm2ba1bs1bs2_instance._get_control_af(gnomad_exomes_empty)
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


def test_get_any_af_control_preference(auto_pm2ba1bs1bs2_instance, gnomad_exomes_mixed_data):
    """Test retrieval of control data when mixed data is available."""
    result = auto_pm2ba1bs1bs2_instance._get_any_af(gnomad_exomes_mixed_data)
    assert result is not None
    assert result.afGrpmax == 0.03
    assert result.cohort == "controls"


def test_get_any_af_highest_non_control(
    auto_pm2ba1bs1bs2_instance, gnomad_exomes_no_controls_higher_af
):
    """Test retrieval of the highest allele frequency from non-control data."""
    result = auto_pm2ba1bs1bs2_instance._get_any_af(gnomad_exomes_no_controls_higher_af)
    assert result is not None
    assert result.afGrpmax == 0.07
    assert result.cohort == "non_controls"


def test_get_any_af_no_data(auto_pm2ba1bs1bs2_instance, gnomad_exomes_empty):
    """Test behavior when no allele count data is available."""
    result = auto_pm2ba1bs1bs2_instance._get_any_af(gnomad_exomes_empty)
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
    auto_pm2ba1bs1bs2_instance, seqvar_mitochondrial, gnomad_mtdna_data, gnomad_exomes_data_af
):
    """Test retrieval of allele frequency for mitochondrial variant."""
    result = auto_pm2ba1bs1bs2_instance._get_af(
        seqvar_mitochondrial, gnomad_mtdna_data, gnomad_exomes_data_af
    )
    assert result == 0.015


def test_get_af_non_mitochondrial_with_controls(
    auto_pm2ba1bs1bs2_instance, seqvar_non_mitochondrial, gnomad_exomes_data_af
):
    """Test retrieval of allele frequency for non-mitochondrial variant with controls data."""
    result = auto_pm2ba1bs1bs2_instance._get_af(
        seqvar_non_mitochondrial, None, gnomad_exomes_data_af
    )
    assert result == 0.02


def test_get_af_non_mitochondrial_without_controls(
    auto_pm2ba1bs1bs2_instance, seqvar_non_mitochondrial, gnomad_exomes_data_af
):
    """Test retrieval of allele frequency for non-mitochondrial variant without controls data."""
    result = auto_pm2ba1bs1bs2_instance._get_af(
        seqvar_non_mitochondrial, None, gnomad_exomes_data_af
    )
    assert result == 0.02


def test_get_af_missing_mtdna_data(
    auto_pm2ba1bs1bs2_instance, seqvar_mitochondrial, gnomad_exomes_data_af
):
    """Test handling of missing mitochondrial gnomad data."""
    with pytest.raises(MissingDataError):
        auto_pm2ba1bs1bs2_instance._get_af(seqvar_mitochondrial, None, gnomad_exomes_data_af)


def test_get_af_missing_exomes_data(
    auto_pm2ba1bs1bs2_instance, seqvar_non_mitochondrial, gnomad_mtdna_data
):
    """Test handling of missing gnomad exomes data."""
    with pytest.raises(MissingDataError):
        auto_pm2ba1bs1bs2_instance._get_af(seqvar_non_mitochondrial, gnomad_mtdna_data, None)


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
#     auto_pm2ba1bs1bs2_instance,
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
#     condition = auto_pm2ba1bs1bs2_instance._get_allele_cond(seqvar)
#     assert condition == AlleleCondition.Dominant


# @patch("src.utils.SeqVarTranscriptsHelper")
# @patch("src.api.annonars.AnnonarsClient.get_gene_info")
# def test_get_allele_cond_with_decipher_high_pHi(
#     mock_annonars_client,
#     mock_transcripts_helper,
#     auto_pm2ba1bs1bs2_instance,
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
#     condition = auto_pm2ba1bs1bs2_instance._get_allele_cond(seqvar)
#     assert condition == AlleleCondition.Dominant


# @patch("src.utils.SeqVarTranscriptsHelper")
# @patch("src.api.annonars.AnnonarsClient.get_gene_info")
# def test_get_allele_cond_with_domino_high_score(
#     mock_annonars_client,
#     mock_transcripts_helper,
#     auto_pm2ba1bs1bs2_instance,
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
#     condition = auto_pm2ba1bs1bs2_instance._get_allele_cond(seqvar)
#     assert condition == AlleleCondition.Dominant
