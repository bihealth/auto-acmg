from unittest.mock import MagicMock, patch

import pytest

from src.api.reev.annonars import AnnonarsClient
from src.defs.annonars_gene import (
    AnnonarsGeneResponse,
    Clingen,
    DecipherHi,
    Domino,
    GeneInfo,
    Genes,
)
from src.defs.auto_acmg import (
    AlleleCondition,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGSeqVarTresholds,
    AutoACMGStrength,
)
from src.defs.exceptions import AlgorithmError, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.auto_pm2_ba1_bs1_bs2 import AutoPM2BA1BS1BS2
from src.utils import SeqVarTranscriptsHelper


@pytest.fixture
def auto_pm2ba1bs1bs2():
    return AutoPM2BA1BS1BS2()


# =========== _get_control_af ===========


@pytest.fixture
def gnomad_exomes_data():
    return MagicMock(
        alleleCounts=[
            MagicMock(cohort="controls", afGrpmax=0.01, anGrpmax=3000),  # Control population
            MagicMock(
                cohort="non_controls", afGrpmax=0.05, anGrpmax=3000
            ),  # Non-control population
        ]
    )


@pytest.fixture
def gnomad_exomes_no_controls():
    return MagicMock(
        alleleCounts=[
            MagicMock(cohort="non_controls", afGrpmax=0.05, anGrpmax=3000)
        ]  # Only non-control data
    )


@pytest.fixture
def gnomad_exomes_empty():
    return MagicMock(alleleCounts=[])  # No data


@pytest.fixture
def var_data():
    thresholds = AutoACMGSeqVarTresholds(ba1_benign=0.05, bs1_benign=0.01, pm2_pathogenic=0.005)
    return AutoACMGSeqVarData(thresholds=thresholds)


def test_get_control_af_success(auto_pm2ba1bs1bs2, gnomad_exomes_data, var_data):
    """Test successful retrieval of control allele frequency."""
    var_data.gnomad_exomes = gnomad_exomes_data
    result = auto_pm2ba1bs1bs2._get_control_af(var_data)
    assert result is not None
    assert result.afGrpmax == 0.01
    assert result.cohort == "controls"


def test_get_control_af_no_controls(auto_pm2ba1bs1bs2, gnomad_exomes_no_controls, var_data):
    """Test the case where no control data is available."""
    var_data.gnomad_exomes = gnomad_exomes_no_controls
    result = auto_pm2ba1bs1bs2._get_control_af(var_data)
    assert result is None


def test_get_control_af_empty_data(auto_pm2ba1bs1bs2, gnomad_exomes_empty, var_data):
    """Test the case where allele counts are empty."""
    var_data.gnomad_exomes = gnomad_exomes_empty
    result = auto_pm2ba1bs1bs2._get_control_af(var_data)
    assert result is None


# =========== _get_any_af ===========


@pytest.fixture
def gnomad_exomes_mixed_data():
    return MagicMock(
        alleleCounts=[
            MagicMock(cohort="non_controls", afGrpmax=0.02, anGrpmax=3000),
            MagicMock(cohort="controls", afGrpmax=0.03, anGrpmax=3000),
            MagicMock(cohort="non_controls", afGrpmax=0.04, anGrpmax=3000),
        ]
    )


@pytest.fixture
def gnomad_exomes_no_controls_higher_af():
    return MagicMock(
        alleleCounts=[
            MagicMock(cohort="non_controls", afGrpmax=0.05, anGrpmax=3000),
            MagicMock(cohort="non_controls", afGrpmax=0.07, anGrpmax=3000),
        ]
    )


def test_get_any_af_control_preference(auto_pm2ba1bs1bs2, gnomad_exomes_mixed_data, var_data):
    """Test retrieval of control data when mixed data is available."""
    var_data.gnomad_exomes = gnomad_exomes_mixed_data
    result = auto_pm2ba1bs1bs2._get_any_af(var_data)
    assert result is not None
    assert result.afGrpmax == 0.03
    assert result.cohort == "controls"


def test_get_any_af_highest_non_control(
    auto_pm2ba1bs1bs2, gnomad_exomes_no_controls_higher_af, var_data
):
    """Test retrieval of the highest allele frequency from non-control data."""
    var_data.gnomad_exomes = gnomad_exomes_no_controls_higher_af
    result = auto_pm2ba1bs1bs2._get_any_af(var_data)
    assert result is not None
    assert result.afGrpmax == 0.07
    assert result.cohort == "non_controls"


def test_get_any_af_no_data(auto_pm2ba1bs1bs2, gnomad_exomes_empty, var_data):
    """Test behavior when no allele count data is available."""
    var_data.gnomad_exomes = gnomad_exomes_empty
    result = auto_pm2ba1bs1bs2._get_any_af(var_data)
    assert result is None


# =========== _get_af ==================


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
            MagicMock(cohort="controls", afGrpmax=0.02, anGrpmax=3000),
            MagicMock(cohort="non_controls", afGrpmax=0.03, anGrpmax=3000),
        ]
    )


@pytest.fixture
def gnomad_exomes_empty_af():
    return MagicMock(alleleCounts=[])


def test_get_af_non_mitochondrial_with_controls(
    auto_pm2ba1bs1bs2, seqvar_non_mitochondrial, gnomad_exomes_data_af, var_data
):
    """Test retrieval of allele frequency for non-mitochondrial variant with controls data."""
    var_data.gnomad_exomes = gnomad_exomes_data_af
    result = auto_pm2ba1bs1bs2._get_af(seqvar_non_mitochondrial, var_data)
    assert result == 0.02


def test_get_af_non_mitochondrial_without_controls(
    auto_pm2ba1bs1bs2, seqvar_non_mitochondrial, gnomad_exomes_data_af, var_data
):
    """Test retrieval of allele frequency for non-mitochondrial variant without controls data."""
    var_data.gnomad_exomes = gnomad_exomes_data_af
    result = auto_pm2ba1bs1bs2._get_af(seqvar_non_mitochondrial, var_data)
    assert result == 0.02


def test_get_af_missing_exomes_data(
    auto_pm2ba1bs1bs2, seqvar_non_mitochondrial, gnomad_mtdna_data, var_data
):
    """Test handling of missing gnomad exomes data."""
    var_data.gnomad_mtdna = gnomad_mtdna_data
    var_data.gnomad_exomes = None
    with pytest.raises(MissingDataError):
        auto_pm2ba1bs1bs2._get_af(seqvar_non_mitochondrial, var_data)


# =========== _get_m_af ==================


@pytest.fixture
def seqvar_mitochondrial():
    return MagicMock(chrom="MT", position=12345)


def test_get_m_af_mitochondrial(
    auto_pm2ba1bs1bs2, seqvar_mitochondrial, gnomad_mtdna_data, var_data
):
    """Test retrieval of allele frequency for mitochondrial variant."""
    var_data.gnomad_mtdna = gnomad_mtdna_data
    result = auto_pm2ba1bs1bs2._get_m_af(var_data)
    assert result == 0.015


def test_get_m_af_missing_mtdna_data(
    auto_pm2ba1bs1bs2, seqvar_mitochondrial, gnomad_exomes_data_af, var_data
):
    """Test handling of missing gnomad mtdna data."""
    var_data.gnomad_mtdna = None
    var_data.gnomad_exomes = gnomad_exomes_data_af
    with pytest.raises(MissingDataError):
        auto_pm2ba1bs1bs2._get_m_af(var_data)


# =========== _get_allele_cond ===========


@pytest.fixture
def gene_info():
    # Base mock structure for gene information responses
    return {
        "Gene123": GeneInfo(
            clingen=Clingen(
                haploinsufficiencyScore=None  # No score set initially
            ),
            domino=Domino(score=None),  # No score set initially
            decipherHi=DecipherHi(pHi=None),  # No score set initially
        )
    }


@patch.object(SeqVarTranscriptsHelper, "initialize")
@patch.object(SeqVarTranscriptsHelper, "get_ts_info")
@patch.object(AnnonarsClient, "get_gene_info")
def test_get_allele_cond_clingen_dosage(
    mock_get_gene_info, mock_get_ts_info, mock_initialize, auto_pm2ba1bs1bs2, seqvar, gene_info
):
    # Modify gene_info for the test case to simulate dominant condition
    gene_info[
        "Gene123"
    ].clingen.haploinsufficiencyScore = "CLINGEN_DOSAGE_SCORE_SUFFICIENT_EVIDENCE_AVAILABLE"
    mock_get_gene_info.return_value = AnnonarsGeneResponse(genes=Genes(root=gene_info))

    mock_get_ts_info.return_value = (None, MagicMock(geneId="Gene123"), None, None, None)
    mock_initialize.return_value = None

    result = auto_pm2ba1bs1bs2._get_allele_cond(seqvar)
    assert result == AlleleCondition.Dominant


@patch.object(SeqVarTranscriptsHelper, "initialize")
@patch.object(SeqVarTranscriptsHelper, "get_ts_info")
@patch.object(AnnonarsClient, "get_gene_info")
def test_get_allele_cond_domino_data(
    mock_get_gene_info, mock_get_ts_info, mock_initialize, auto_pm2ba1bs1bs2, seqvar, gene_info
):
    # Modify gene_info for the test case to simulate Domino data
    gene_info["Gene123"].domino.score = 0.6
    mock_get_gene_info.return_value = AnnonarsGeneResponse(genes=Genes(root=gene_info))

    mock_get_ts_info.return_value = (None, MagicMock(geneId="Gene123"), None, None, None)
    mock_initialize.return_value = None

    result = auto_pm2ba1bs1bs2._get_allele_cond(seqvar)
    assert result == AlleleCondition.Dominant


@pytest.mark.parametrize(
    "clingen_score, domino_score, decipher_score, expected_condition",
    [
        (
            "CLINGEN_DOSAGE_SCORE_SUFFICIENT_EVIDENCE_AVAILABLE",
            None,
            None,
            AlleleCondition.Dominant,
        ),
        ("CLINGEN_DOSAGE_SCORE_RECESSIVE", None, None, AlleleCondition.Recessive),
        ("CLINGEN_DOSAGE_SCORE_NO_EVIDENCE_AVAILABLE", 0.6, None, AlleleCondition.Dominant),
        (None, 0.1, None, AlleleCondition.Recessive),
        (None, None, 0.95, AlleleCondition.Dominant),
        (None, None, 0.85, AlleleCondition.Unknown),
    ],
)
@patch.object(SeqVarTranscriptsHelper, "initialize")
@patch.object(SeqVarTranscriptsHelper, "get_ts_info")
@patch.object(AnnonarsClient, "get_gene_info")
def test_get_allele_cond_various_conditions(
    mock_get_gene_info,
    mock_get_ts_info,
    mock_initialize,
    auto_pm2ba1bs1bs2,
    seqvar,
    clingen_score,
    domino_score,
    decipher_score,
    expected_condition,
    gene_info,
):
    # Modify gene_info based on parameters
    if clingen_score:
        gene_info["Gene123"].clingen.haploinsufficiencyScore = clingen_score
    if domino_score is not None:
        gene_info["Gene123"].domino.score = domino_score
    if decipher_score is not None:
        gene_info["Gene123"].decipherHi.pHi = decipher_score

    mock_get_gene_info.return_value = AnnonarsGeneResponse(genes=Genes(root=gene_info))
    mock_get_ts_info.return_value = (None, MagicMock(geneId="Gene123"), None, None, None)
    mock_initialize.return_value = None

    result = auto_pm2ba1bs1bs2._get_allele_cond(seqvar)
    assert result == expected_condition


@patch("src.utils.SeqVarTranscriptsHelper")
def test_get_allele_cond_no_gene_found(mock_seqvar_transcripts_helper, auto_pm2ba1bs1bs2, seqvar):
    mock_helper_instance = mock_seqvar_transcripts_helper.return_value
    mock_helper_instance.initialize.return_value = None
    mock_helper_instance.get_ts_info.return_value = (None, None, None, None, None)

    with pytest.raises(AlgorithmError, match="No gene found for the transcript."):
        auto_pm2ba1bs1bs2._get_allele_cond(seqvar)


# ========== _check_zyg ===========


@pytest.fixture
def seqvar_mitochondrial_zyg():
    return MagicMock(chrom="M", pos=100, delete="A", insert="T")


@patch(
    "src.seqvar.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._get_allele_cond",
    return_value=AlleleCondition.Dominant,
)
@patch("src.seqvar.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._get_control_af")
@patch("src.seqvar.auto_pm2_ba1_bs1_bs2.AutoPM2BA1BS1BS2._get_any_af")
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
    assert AutoPM2BA1BS1BS2()._ba1_exception(seqvar_in_exception) is True


def test_ba1_exception_false(seqvar_not_in_exception):
    """Test that _ba1_exception returns False for a variant not in the exception list."""
    assert AutoPM2BA1BS1BS2()._ba1_exception(seqvar_not_in_exception) is False


# =========== _bs2_not_applicable ==============


def test_bs2_not_applicable(auto_pm2ba1bs1bs2, var_data):
    """Test that BS2 is applicable per default."""
    assert auto_pm2ba1bs1bs2._bs2_not_applicable(var_data) is False


# =========== verify_pm2ba1bs1bs2 ==============


@pytest.fixture
def seqvar_mt():
    return SeqVar(genome_release=GenomeRelease.GRCh38, chrom="MT", pos=1000, delete="A", insert="G")


@patch.object(AutoPM2BA1BS1BS2, "_get_af", return_value=0.0001)
@patch.object(AutoPM2BA1BS1BS2, "_ba1_exception", return_value=False)
@patch.object(AutoPM2BA1BS1BS2, "_bs2_not_applicable", return_value=False)
@patch.object(AutoPM2BA1BS1BS2, "_check_zyg", return_value=True)
def test_verify_pm2ba1bs1bs2(
    mock_get_af,
    mock_ba1_exception,
    mock_bs2_na,
    mock_check_zyg,
    auto_pm2ba1bs1bs2,
    seqvar,
    var_data,
):
    result, comment = auto_pm2ba1bs1bs2.verify_pm2ba1bs1bs2(seqvar, var_data)
    assert result.PM2 is True
    assert "PM2 is met" in comment


@patch.object(AutoPM2BA1BS1BS2, "_get_m_af", return_value=0.00002)
def test_verify_pm2ba1bs1bs2_mitochondrial(mock_get_m_af, auto_pm2ba1bs1bs2, seqvar_mt, var_data):
    result, comment = auto_pm2ba1bs1bs2.verify_pm2ba1bs1bs2(seqvar_mt, var_data)
    assert result.PM2 is True
    assert "Allele frequency <= 0.002%: PM2 is met" in comment


# =========== predict_pm2ba1bs1bs2 ==============


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


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
    mock_verify.return_value = (
        create_pred_object(True, True, True, True),
        "All criteria met.",
    )
    results = auto_pm2ba1bs1bs2.predict_pm2ba1bs1bs2(seqvar, var_data)
    assert all([r.prediction == AutoACMGPrediction.Applicable for r in results])
    assert all(["All criteria met." in r.summary for r in results])


@patch.object(AutoPM2BA1BS1BS2, "verify_pm2ba1bs1bs2")
def test_mixed_conditions(mock_verify, auto_pm2ba1bs1bs2, seqvar, var_data):
    mock_verify.return_value = (
        create_pred_object(True, False, True, False),
        "Mixed conditions.",
    )
    results = auto_pm2ba1bs1bs2.predict_pm2ba1bs1bs2(seqvar, var_data)
    assert results[0].prediction == AutoACMGPrediction.Applicable
    assert results[1].prediction == AutoACMGPrediction.NotApplicable
    assert results[2].prediction == AutoACMGPrediction.Applicable
    assert results[3].prediction == AutoACMGPrediction.NotApplicable


@patch.object(AutoPM2BA1BS1BS2, "verify_pm2ba1bs1bs2")
def test_all_criteria_not_met(mock_verify, auto_pm2ba1bs1bs2, seqvar, var_data):
    mock_verify.return_value = (
        create_pred_object(False, False, False, False),
        "No criteria met.",
    )
    results = auto_pm2ba1bs1bs2.predict_pm2ba1bs1bs2(seqvar, var_data)
    assert all([r.prediction == AutoACMGPrediction.NotApplicable for r in results])


@patch.object(AutoPM2BA1BS1BS2, "verify_pm2ba1bs1bs2")
def test_failure_condition(mock_verify, auto_pm2ba1bs1bs2, seqvar, var_data):
    mock_verify.return_value = (None, "Error occurred.")
    results = auto_pm2ba1bs1bs2.predict_pm2ba1bs1bs2(seqvar, var_data)
    assert all([r.prediction == AutoACMGPrediction.Failed for r in results])
    assert all(["Error occurred." in r.summary for r in results])
