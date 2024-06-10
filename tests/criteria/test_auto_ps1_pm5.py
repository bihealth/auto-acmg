from unittest.mock import MagicMock, patch

import pytest

from src.api.annonars import AnnonarsClient
from src.criteria.auto_ps1_pm5 import AutoPS1PM5
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.auto_acmg import AminoAcid
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from tests.utils import get_json_object


@pytest.fixture
def seqvar():
    return SeqVar(GenomeRelease.GRCh38, "1", 1000, "A", "T", "1:1000A>T")


@pytest.fixture
def variant_info():
    return AnnonarsVariantResponse.model_validate(
        get_json_object("annonars/annonars_variant_success.json")
    )


@pytest.fixture
def auto_ps1_pm5(seqvar, variant_info):
    return AutoPS1PM5(seqvar, GenomeRelease.GRCh38, variant_info.result)


@patch.object(AnnonarsClient, "get_variant_info")
def test_auto_ps1_pm5_get_variant_info(mock_get_variant_info, variant_info, seqvar):
    """Test getting variant information."""
    mock_get_variant_info.return_value = variant_info

    auto_ps1_pm5 = AutoPS1PM5(
        seqvar=seqvar, genome_release=GenomeRelease.GRCh38, variant_info=variant_info.result
    )
    response = auto_ps1_pm5._get_variant_info(auto_ps1_pm5.seqvar)
    assert response is not None
    assert response == variant_info


@patch.object(AnnonarsClient, "get_variant_info")
def test_auto_ps1_pm5_get_variant_info_none(mock_get_variant_info, variant_info, seqvar):
    """Test getting variant information returns None."""
    mock_get_variant_info.return_value = None

    auto_ps1_pm5 = AutoPS1PM5(
        seqvar=seqvar, genome_release=GenomeRelease.GRCh38, variant_info=variant_info.result
    )
    response = auto_ps1_pm5._get_variant_info(auto_ps1_pm5.seqvar)
    assert response is None


@patch.object(AnnonarsClient, "get_variant_info")
def test_auto_ps1_pm5_get_variant_info_auto_acmg_exception(
    mock_get_variant_info, variant_info, seqvar
):
    """Test getting variant information raises AutoAcmgBaseException."""
    mock_get_variant_info.side_effect = AutoAcmgBaseException("An error occurred")

    auto_ps1_pm5 = AutoPS1PM5(
        seqvar=seqvar, genome_release=GenomeRelease.GRCh38, variant_info=variant_info.result
    )
    response = auto_ps1_pm5._get_variant_info(auto_ps1_pm5.seqvar)
    assert response is None


@pytest.mark.parametrize(
    "pHGVSp, expected_result",
    [
        ("p.Thr1399Pro", AminoAcid.Pro),  # Valid change
        ("p.Gly12Ser", AminoAcid.Ser),  # Valid simple change
        ("p.Ala10Ala", AminoAcid.Ala),  # No change in amino acid
        ("p.?1234XYZ", None),  # Invalid format
        ("p.Val100*", None),  # To stop codon
        ("", None),  # Empty string
    ],
)
def test_parse_HGVSp(pHGVSp, expected_result, auto_ps1_pm5):
    """Test parsing of protein changes from HGVS format."""
    result = auto_ps1_pm5._parse_HGVSp(pHGVSp)
    assert result == expected_result


@pytest.mark.parametrize(
    "clinical_significance, expected_result",
    [
        ("CLINICAL_SIGNIFICANCE_PATHOGENIC", True),
        ("CLINICAL_SIGNIFICANCE_BENIGN", False),
        (None, False),  # No clinical significance data
    ],
)
def test_is_pathogenic(clinical_significance, expected_result, auto_ps1_pm5):
    """Test determination of pathogenicity based on variant info."""
    variant_info = MagicMock()
    variant_info.clinvar = MagicMock()
    variant_info.clinvar.referenceAssertions = [MagicMock()]
    variant_info.clinvar.referenceAssertions[0].clinicalSignificance = clinical_significance

    result = auto_ps1_pm5._is_pathogenic(variant_info)
    assert (
        result == expected_result
    ), f"Expected {expected_result} for clinical significance: {clinical_significance}"
