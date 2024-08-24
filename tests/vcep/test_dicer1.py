from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.vcep import DICER1Predictor


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="14", pos=100, delete="A", insert="T")


@pytest.fixture
def dicer1_predictor(seqvar):
    result = MagicMock()  # Mocking the AutoACMGResult object
    return DICER1Predictor(seqvar=seqvar, result=result, config=MagicMock())


@pytest.fixture
def auto_acmg_data():
    return AutoACMGData()


def test_predict_pm1_moderate_criteria_residue(dicer1_predictor, auto_acmg_data):
    """Test when the variant affects a metal ion-binding residue in DICER1."""
    auto_acmg_data.prot_pos = 1705  # Metal ion-binding residue
    auto_acmg_data.hgnc_id = "HGNC:17098"  # DICER1 gene
    result = dicer1_predictor.predict_pm1(dicer1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for metal ion-binding residues."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "affects a metal ion-binding residue" in result.summary
    ), "The summary should indicate the metal ion-binding residue."


def test_predict_pm1_supporting_criteria_domain(dicer1_predictor, auto_acmg_data):
    """Test when the variant affects a residue within the RNase IIIb domain but outside metal ion-binding residues."""
    auto_acmg_data.prot_pos = (
        1720  # Within the RNase IIIb domain but not a metal ion-binding residue
    )
    auto_acmg_data.hgnc_id = "HGNC:17098"  # DICER1 gene
    result = dicer1_predictor.predict_pm1(dicer1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for residues within the RNase IIIb domain."
    assert (
        result.strength == AutoACMGStrength.PathogenicSupporting
    ), "The strength should be PathogenicSupporting."
    assert (
        "affects a residue in the RNase IIIb domain" in result.summary
    ), "The summary should indicate the RNase IIIb domain."


def test_predict_pm1_outside_critical_region(dicer1_predictor, auto_acmg_data):
    """Test when the variant does not affect any critical residues or domains in DICER1."""
    auto_acmg_data.prot_pos = 1600  # Outside any critical region or domain
    auto_acmg_data.hgnc_id = "HGNC:17098"  # DICER1 gene
    result = dicer1_predictor.predict_pm1(dicer1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met for non-critical regions."
    assert (
        result.strength == AutoACMGStrength.PathogenicModerate
    ), "The strength should be PathogenicModerate."
    assert (
        "does not affect a critical domain" in result.summary
    ), "The summary should indicate no critical domain was affected."


def test_predict_pm1_not_applicable_gene(dicer1_predictor, auto_acmg_data):
    """Test when PM1 is not applicable because the gene is not DICER1."""
    auto_acmg_data.hgnc_id = "HGNC:99999"  # Gene not in PM1_CLUSTER mapping
    result = dicer1_predictor.predict_pm1(dicer1_predictor.seqvar, auto_acmg_data)

    assert isinstance(result, AutoACMGCriteria), "The result should be of type AutoACMGCriteria."
    assert (
        result.prediction == AutoACMGPrediction.NotMet
    ), "PM1 should not be met if gene is not in PM1_CLUSTER."
    assert (
        "does not affect a critical domain for DICER1" in result.summary
    ), "The summary should indicate gene is not applicable."


def test_predict_pm1_edge_case_start_boundary_moderate(dicer1_predictor, auto_acmg_data):
    """Test when the variant falls exactly on the start boundary of a moderate criteria residue."""
    auto_acmg_data.prot_pos = 1344  # Start boundary of a metal ion-binding residue
    auto_acmg_data.hgnc_id = "HGNC:17098"  # DICER1 gene
    result = dicer1_predictor.predict_pm1(dicer1_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for the start boundary of a metal ion-binding residue."
    assert (
        "affects a metal ion-binding residue" in result.summary
    ), "The summary should indicate the metal ion-binding residue."


def test_predict_pm1_edge_case_end_boundary_supporting(dicer1_predictor, auto_acmg_data):
    """Test when the variant falls exactly on the end boundary of a supporting criteria domain."""
    auto_acmg_data.prot_pos = 1846  # End boundary of the RNase IIIb domain
    auto_acmg_data.hgnc_id = "HGNC:17098"  # DICER1 gene
    result = dicer1_predictor.predict_pm1(dicer1_predictor.seqvar, auto_acmg_data)

    assert (
        result.prediction == AutoACMGPrediction.Met
    ), "PM1 should be met for the end boundary of the RNase IIIb domain."
    assert (
        "affects a residue in the RNase IIIb domain" in result.summary
    ), "The summary should indicate the RNase IIIb domain."
