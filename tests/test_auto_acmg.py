from unittest.mock import MagicMock, patch

import pytest

from src.auto_acmg import AutoACMG
from src.defs.auto_acmg import AutoACMGResult
from src.defs.exceptions import AutoAcmgBaseException
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


@pytest.fixture
def auto_acmg():
    return AutoACMG(variant_name="test_variant")


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh38, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def mock_variant_result():
    variant_result = MagicMock()
    variant_result.cadd = MagicMock(
        ConsDetail="cadd_consequence",
        Consequence="some_consequence",
        verPhyloP=0.8,
        SpliceAI_acc_gain=0.2,
        SpliceAI_acc_loss=0.3,
        SpliceAI_don_gain=0.4,
        SpliceAI_don_loss=0.5,
        dbscSNV_ada=0.6,
        dbscSNV_rf=0.7,
    )
    variant_result.dbnsfp = MagicMock(
        HGVSp_VEP="p.Val600Glu",
        AlphaMissense_score="0.1",
        MetaRNN_score="0.2",
        BayesDel_noAF_score="0.3",
        REVEL_score="0.4",
        phyloP100way_vertebrate="0.5",
    )
    variant_result.dbscsnv = MagicMock(ada_score=0.8, rf_score=0.9)
    variant_result.gnomad_exomes = MagicMock()
    variant_result.gnomad_mtdna = MagicMock()
    return variant_result


@pytest.fixture
def mock_seqvar_transcript():
    transcript = MagicMock()
    transcript.consequences = ["consequence"]
    transcript.gene_symbol = "GENE"
    transcript.gene_id = "HGNC:1234"
    transcript.feature_id = "ENST00000367770"
    transcript.feature_tag = "MANE"
    transcript.tx_pos = MagicMock(ord=10)
    transcript.protein_pos = MagicMock(ord=11, total=20)
    return transcript


# Mocking dependencies and inherited methods
@patch("src.auto_acmg.SeqVarTranscriptsHelper.get_ts_info")
@patch("src.auto_acmg.SeqVarResolver.resolve_seqvar")
@patch("src.auto_acmg.AutoACMG.parse_data")
@patch("src.auto_acmg.SeqVarTranscriptsHelper.initialize")
@patch("src.auto_acmg.AutoACMG._get_variant_info")
@patch("src.auto_acmg.AutoACMG.predict_pvs1")
@patch("src.auto_acmg.AutoACMG.predict_ps1pm5")
@patch("src.auto_acmg.AutoACMG.predict_pm1")
@patch("src.auto_acmg.AutoACMG.predict_pm2ba1bs1bs2")
@patch("src.auto_acmg.AutoACMG.predict_pm4bp3")
@patch("src.auto_acmg.AutoACMG.predict_pp2bp1")
@patch("src.auto_acmg.AutoACMG.predict_pp3bp4")
@patch("src.auto_acmg.AutoACMG.predict_bp7")
def test_predict(
    mock_predict_bp7,
    mock_predict_pp3bp4,
    mock_predict_pp2bp1,
    mock_predict_pm4bp3,
    mock_predict_pm2ba1bs1bs2,
    mock_predict_pm1,
    mock_predict_ps1pm5,
    mock_predict_pvs1,
    mock_get_variant_info,
    mock_initialize,
    mock_get_ts_info,
    mock_resolve_seqvar,
    mock_parse_data,
    auto_acmg,
    seqvar,
    mock_variant_result,
    mock_seqvar_transcript,
):
    """Test the full prediction flow of AutoACMG."""
    mock_resolve_seqvar.return_value = seqvar
    mock_get_variant_info.return_value = mock_variant_result
    mock_get_ts_info.return_value = (mock_seqvar_transcript, mock_seqvar_transcript, None, [], None)
    mock_parse_data.return_value = MagicMock()

    # Set return values for all the prediction methods
    mock_predict_pvs1.return_value = MagicMock()
    mock_predict_ps1pm5.return_value = (MagicMock(), MagicMock())
    mock_predict_pm1.return_value = MagicMock()
    mock_predict_pm2ba1bs1bs2.return_value = (MagicMock(), MagicMock(), MagicMock(), MagicMock())
    mock_predict_pm4bp3.return_value = (MagicMock(), MagicMock())
    mock_predict_pp2bp1.return_value = (MagicMock(), MagicMock())
    mock_predict_pp3bp4.return_value = (MagicMock(), MagicMock())
    mock_predict_bp7.return_value = MagicMock()

    result = auto_acmg.predict()

    assert isinstance(result, AutoACMGResult), "Result should be of type AutoACMGResult."

    # Verify predictions for different criteria
    assert result.criteria.pvs1 is not None, "PVS1 should have a prediction."
    assert result.criteria.ps1 is not None, "PS1 should have a prediction."
    assert result.criteria.pm1 is not None, "PM1 should have a prediction."
    assert result.criteria.pm2 is not None, "PM2 should have a prediction."
    assert result.criteria.bp3 is not None, "BP3 should have a prediction."
    assert result.criteria.bp7 is not None, "BP7 should have a prediction."


@patch(
    "src.auto_acmg.SeqVarResolver.resolve_seqvar", side_effect=AutoAcmgBaseException("ParseError")
)
def test_resolve_variant_exception(mock_resolve_seqvar, auto_acmg):
    """Test handling of exceptions during variant resolution."""
    result = auto_acmg.resolve_variant()
    assert result is None, "Should return None if an exception occurs during resolution."


@patch("src.auto_acmg.SeqVarResolver.resolve_seqvar", return_value=None)
def test_predict_variant_resolution_failure(mock_resolve_seqvar, auto_acmg):
    """Test predict method when variant resolution fails."""
    result = auto_acmg.predict()
    assert result is None, "Should return None if variant resolution fails."


@patch("src.auto_acmg.SeqVarTranscriptsHelper.initialize")
@patch("src.auto_acmg.AutoACMG._get_variant_info", return_value=None)
def test_parse_data_variant_info_failure(mock_get_variant_info, mock_initialize, auto_acmg, seqvar):
    """Test parse_data method when getting variant information fails."""
    with pytest.raises(Exception):
        auto_acmg.parse_data(seqvar)


@patch("src.auto_acmg.AutoACMG._convert_score_val", return_value=0.5)
def test_convert_score_val(mock_convert_score_val, auto_acmg):
    """Test convert score value functionality."""
    score_value = "0.1;0.3;0.5"
    result = auto_acmg._convert_score_val(score_value)
    assert result == 0.5, "Should return the maximum score value."


# Test individual methods like _get_variant_info, _convert_score_val, resolve_variant, etc.
@patch("src.auto_acmg.AnnonarsClient.get_variant_info")
def test_get_variant_info(mock_get_variant_info, auto_acmg, seqvar, mock_variant_result):
    """Test _get_variant_info method."""
    mock_get_variant_info.return_value.result = mock_variant_result
    result = auto_acmg._get_variant_info(seqvar)
    assert result == mock_variant_result, "Should return the variant information from Annonars."


@patch("src.auto_acmg.AnnonarsClient.get_variant_info", side_effect=AutoAcmgBaseException("Error"))
def test_get_variant_info_failure(mock_get_variant_info, auto_acmg, seqvar):
    """Test _get_variant_info method failure case."""
    result = auto_acmg._get_variant_info(seqvar)
    assert result is None, "Should return None if getting variant info fails."
