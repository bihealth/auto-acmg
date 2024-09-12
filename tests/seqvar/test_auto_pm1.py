from unittest.mock import MagicMock, patch

import pytest
import tabix

from src.api.reev.annonars import AnnonarsClient
from src.defs.auto_acmg import (
    PM1,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    GenomicStrand,
)
from src.defs.exceptions import AlgorithmError, InvalidAPIResposeError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.auto_pm1 import AutoPM1


@pytest.fixture
def seqvar():
    return SeqVar(genome_release=GenomeRelease.GRCh37, chrom="1", pos=100, delete="A", insert="T")


@pytest.fixture
def auto_pm1():
    return AutoPM1()


# ============== _count_vars =================


@pytest.fixture
def mock_response():
    response = MagicMock()
    variant1 = MagicMock(
        records=[
            MagicMock(
                classifications=MagicMock(
                    germlineClassification=MagicMock(description="Pathogenic")
                ),
                variationType="VARIATION_TYPE_SNV",
            )
        ]
    )
    variant2 = MagicMock(
        records=[
            MagicMock(
                classifications=MagicMock(
                    germlineClassification=MagicMock(description="Likely benign")
                ),
                variationType="VARIATION_TYPE_SNV",
            )
        ]
    )
    response.clinvar = [variant1, variant2]
    return response


@pytest.fixture
def auto_acmg_data_plus():
    exons = [
        MagicMock(altStartI=10000, altEndI=20000),
        MagicMock(altStartI=30000, altEndI=40000),
        MagicMock(altStartI=50000, altEndI=60000),
    ]
    data = AutoACMGSeqVarData()
    data.exons = exons  # type: ignore
    data.strand = GenomicStrand.Plus
    return data


@pytest.fixture
def auto_acmg_data_minus():
    exons = [
        MagicMock(altStartI=10000, altEndI=20000),
        MagicMock(altStartI=30000, altEndI=40000),
        MagicMock(altStartI=50000, altEndI=60000),
    ]
    data = AutoACMGSeqVarData()
    data.exons = exons  # type: ignore
    data.strand = GenomicStrand.Minus
    return data


def test_get_affected_exon_plus_strand_in_exon(seqvar, auto_pm1, auto_acmg_data_plus):
    """Test for a variant within an exon on the plus strand."""
    seqvar.pos = 35000  # Position within the second exon
    result = auto_pm1._get_affected_exon(auto_acmg_data_plus, seqvar)
    assert result == 2, "Exon number should be 2 for a position within the second exon."


def test_get_affected_exon_plus_strand_before_exon(seqvar, auto_pm1, auto_acmg_data_plus):
    """Test for a variant before any exon on the plus strand."""
    seqvar.pos = 9000  # Position before the first exon
    result = auto_pm1._get_affected_exon(auto_acmg_data_plus, seqvar)
    assert result == 0, "Exon number should be 0 for a position before any exon."


def test_get_affected_exon_plus_strand_after_exon(seqvar, auto_pm1, auto_acmg_data_plus):
    """Test for a variant after the last exon on the plus strand."""
    seqvar.pos = 70000  # Position after the last exon
    result = auto_pm1._get_affected_exon(auto_acmg_data_plus, seqvar)
    assert result == 3, "Exon number should be 3 for a position after the last exon."


def test_get_affected_exon_minus_strand_in_exon(seqvar, auto_pm1, auto_acmg_data_minus):
    """Test for a variant within an exon on the minus strand."""
    seqvar.pos = 35000  # Position within the second exon (reverse order for minus strand)
    result = auto_pm1._get_affected_exon(auto_acmg_data_minus, seqvar)
    assert (
        result == 2
    ), "Exon number should be 2 for a position within the second exon (minus strand)."


def test_get_affected_exon_minus_strand_before_exon(seqvar, auto_pm1, auto_acmg_data_minus):
    """Test for a variant before any exon on the minus strand."""
    seqvar.pos = 65000  # Position before the first exon (reverse order for minus strand)
    result = auto_pm1._get_affected_exon(auto_acmg_data_minus, seqvar)
    assert result == 0, "Exon number should be 0 for a position before any exon (minus strand)."


def test_get_affected_exon_minus_strand_after_exon(seqvar, auto_pm1, auto_acmg_data_minus):
    """Test for a variant after the last exon on the minus strand."""
    seqvar.pos = 5000  # Position after the last exon (reverse order for minus strand)
    result = auto_pm1._get_affected_exon(auto_acmg_data_minus, seqvar)
    assert result == 3, "Exon number should be 3 for a position after the last exon (minus strand)."


def test_get_affected_exon_not_set_strand(seqvar, auto_pm1, auto_acmg_data_plus):
    """Test for an invalid strand scenario."""
    auto_acmg_data_plus.strand = GenomicStrand.NotSet
    with pytest.raises(AlgorithmError):
        auto_pm1._get_affected_exon(auto_acmg_data_plus, seqvar)


@patch.object(AnnonarsClient, "get_variant_from_range")
def test_count_vars_success(mock_get_variant_from_range, auto_pm1, seqvar, mock_response):
    """Test counting pathogenic and benign variants successfully."""
    mock_get_variant_from_range.return_value = mock_response
    pathogenic, benign = auto_pm1._count_vars(seqvar, 75, 125)
    assert pathogenic == 1
    assert benign == 1


def test_count_vars_invalid_range(auto_pm1, seqvar):
    """Test error handling when the end position is less than the start position."""
    with pytest.raises(AlgorithmError) as excinfo:
        auto_pm1._count_vars(seqvar, 150, 100)
    assert "End position is less than the start position" in str(excinfo.value)


@patch.object(AnnonarsClient, "get_variant_from_range")
def test_count_vars_api_failure(mock_get_variant_from_range, auto_pm1, seqvar):
    """Test handling API response failures."""
    mock_get_variant_from_range.return_value = None
    with pytest.raises(InvalidAPIResposeError) as excinfo:
        auto_pm1._count_vars(seqvar, 75, 125)
    assert "Failed to get variant from range. No ClinVar data." in str(excinfo.value)


@patch.object(AnnonarsClient, "get_variant_from_range")
def test_count_vars_no_clinvar_data(mock_get_variant_from_range, auto_pm1, seqvar):
    """Test handling missing ClinVar data in the response."""
    response = MagicMock(clinvar=[])
    mock_get_variant_from_range.return_value = response
    with pytest.raises(InvalidAPIResposeError) as excinfo:
        auto_pm1._count_vars(seqvar, 75, 125)
    assert "Failed to get variant from range. No ClinVar data." in str(excinfo.value)


# ============== _get_uniprot_domain ==============


@patch("src.seqvar.auto_pm1.tabix")
def test_get_uniprot_domain_found(mock_tabix, auto_pm1, seqvar):
    """Test retrieving UniProt domain successfully."""
    # Setting up the mock
    mock_query = MagicMock()
    mock_query.__iter__.return_value = iter([["1", "50", "150"]])
    mock_open = mock_tabix.open.return_value
    mock_open.query.return_value = mock_query

    domain = auto_pm1._get_uniprot_domain(seqvar)
    assert domain == (50, 150), "Should return the correct start and end positions"


@patch("src.seqvar.auto_pm1.tabix")
def test_get_uniprot_domain_none_found(mock_tabix, auto_pm1, seqvar):
    """Test no UniProt domain found for the variant."""
    # Setting up the mock to return an empty iterator
    mock_query = MagicMock()
    mock_query.__iter__.return_value = iter([])
    mock_open = mock_tabix.open.return_value
    mock_open.query.return_value = mock_query

    domain = auto_pm1._get_uniprot_domain(seqvar)
    assert domain is None, "Should return None when no domain is found"


@patch("src.seqvar.auto_pm1.tabix.open")
def test_get_uniprot_domain_exception(mock_tabix_open, auto_pm1, seqvar):
    """Test handling exceptions when querying UniProt domains."""
    mock_tabix_open.side_effect = tabix.TabixError("Test error")

    with pytest.raises(AlgorithmError) as excinfo:
        auto_pm1._get_uniprot_domain(seqvar)
    assert "Failed to check if the variant is in a UniProt domain." in str(excinfo.value)


# ============== verify_pm1 ==============


@pytest.fixture
def var_data():
    thresholds = MagicMock(pm1_pathogenic=2)
    return MagicMock(thresholds=thresholds)


def test_verify_pm1_mitochondrial(auto_pm1, var_data):
    """Test mitochondrial chromosome handling"""
    seqvar_mt = SeqVar(
        genome_release=GenomeRelease.GRCh37, chrom="MT", pos=100, delete="A", insert="T"
    )
    result, comment = auto_pm1.verify_pm1(seqvar_mt, var_data)
    assert not result.PM1, "PM1 should not be met for mitochondrial variants"
    assert (
        "mitochondrial genome" in comment
    ), "Should note that the variant is in the mitochondrial genome"


@patch.object(AutoPM1, "_count_vars", return_value=(3, 1))  # More than threshold
@patch.object(AutoPM1, "_get_uniprot_domain", return_value=None)  # No domain found
def test_verify_pm1_pathogenic_count_meets(
    mock_count_vars, mock_get_uniprot_domain, auto_pm1, seqvar, var_data
):
    """Test pathogenic count meets the threshold"""
    result, comment = auto_pm1.verify_pm1(seqvar, var_data)
    assert result.PM1, "PM1 should be met when the pathogenic variant count is sufficient"
    assert "PM1 is met" in comment, "Comment should indicate PM1 is met"


@patch.object(
    AutoPM1, "_count_vars", side_effect=[(0, 1), (3, 1)]
)  # First call no pathogenics, second call meets criteria
@patch.object(AutoPM1, "_get_uniprot_domain", return_value=(100, 105))  # Domain found
def test_verify_pm1_uniprot_domain_handling(
    mock_count_vars, mock_get_uniprot_domain, auto_pm1, seqvar, var_data
):
    """Test UniProt domain handling and pathogenic count check"""
    result, comment = auto_pm1.verify_pm1(seqvar, var_data)
    assert result.PM1, "PM1 should be met when the pathogenic count in the domain is sufficient"
    assert "in Uniprot Domain" in comment, "Comment should mention checking within a UniProt domain"


# ============== predict_pm1 ==============


@pytest.fixture
def var_data_predict():
    thresholds = MagicMock(pm1_pathogenic=3)
    data = MagicMock(thresholds=thresholds)
    return data


@patch.object(AutoPM1, "verify_pm1")
def test_predict_pm1_met(mock_verify_pm1, auto_pm1, seqvar, var_data_predict):
    """Test predicting PM1 criteria met"""
    mock_verify_pm1.return_value = (
        MagicMock(PM1=True, PM1_strength=AutoACMGStrength.PathogenicModerate),
        "Criterion met.",
    )
    result = auto_pm1.predict_pm1(seqvar, var_data_predict)
    assert result.prediction == AutoACMGPrediction.Applicable
    assert result.strength == AutoACMGStrength.PathogenicModerate
    assert "Criterion met." in result.summary


@patch.object(AutoPM1, "verify_pm1")
def test_predict_pm1_not_met(mock_verify_pm1, auto_pm1, seqvar, var_data_predict):
    """Test predicting PM1 criteria not met"""
    mock_verify_pm1.return_value = (
        PM1(PM1=False, PM1_strength=AutoACMGStrength.PathogenicModerate),
        "Criterion not met.",
    )
    result = auto_pm1.predict_pm1(seqvar, var_data_predict)
    assert result.prediction == AutoACMGPrediction.NotApplicable
    assert result.strength == AutoACMGStrength.PathogenicModerate
    assert "Criterion not met." in result.summary


@patch.object(AutoPM1, "verify_pm1")
def test_predict_pm1_failed(mock_verify_pm1, auto_pm1, seqvar, var_data_predict):
    """Test predicting PM1 criteria failed"""
    mock_verify_pm1.return_value = (None, "Error during prediction.")
    result = auto_pm1.predict_pm1(seqvar, var_data_predict)
    assert result.prediction == AutoACMGPrediction.Failed
    assert result.strength == AutoACMGStrength.PathogenicModerate
    assert "Error during prediction." in result.summary
