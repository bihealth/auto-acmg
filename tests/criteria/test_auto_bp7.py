from unittest.mock import patch

import pytest

from src.api.annonars import AnnonarsClient
from src.criteria.auto_bp7 import AutoBP7
from src.defs.annonars_range import AnnonarsRangeResponse
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from tests.utils import get_json_object


@pytest.fixture
def seqvar():
    return SeqVar(GenomeRelease.GRCh38, "1", 1000, "A", "T", "1:1000A>T")


@pytest.fixture
def variant_info():
    return AnnonarsVariantResponse.model_validate(get_json_object("annonars/RP11_variant.json"))


@pytest.fixture
def auto_bp7(seqvar, variant_info):
    return AutoBP7(seqvar, variant_info.result)


@pytest.mark.parametrize(
    "file_path, expected",
    [
        ("annonars/GAA_range.json", True),
        ("annonars/PCSK9_range.json", False),
        ("annonars/CDH1_range.json", False),
    ],
)
def test_check_proximity_to_pathogenic_variants(file_path, expected, seqvar, variant_info):
    """Test check_proximity_to_pathogenic_variants."""
    range_response = AnnonarsRangeResponse.model_validate(get_json_object(file_path))
    with patch.object(AnnonarsClient, "get_variant_from_range", return_value=range_response):
        auto_bp7 = AutoBP7(seqvar=seqvar, variant_info=variant_info.result)
        response = auto_bp7._check_proximity_to_pathogenic_variants(seqvar)
        assert response == expected
