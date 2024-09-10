import pytest

from src.defs.annonars_range import AnnonarsRangeResponse
from src.defs.annonars_variant import AnnonarsVariantResponse
from tests.utils import get_json_object

# ----------- annoanrs_range.py -----------


@pytest.mark.parametrize(
    "json_file",
    [
        "annonars/PCSK9_range.json",
        "annonars/range_failure.json",
        "annonars/GAA_range.json",
        "annonars/CDH1_range.json",
        "annonars/SHH_range.json",
    ],
)
def test_annonars_range_response_model(json_file):
    """Test AnnonarsRangeResponse model for various example responses."""
    annonars_response = get_json_object(json_file)
    assert AnnonarsRangeResponse.model_validate(annonars_response)


# ----------- annonars_variant.py -----------


@pytest.mark.parametrize(
    "json_file",
    [
        "annonars/RP11_variant.json",
        "annonars/variant_failure.json",
        "annonars/LARP7_variant.json",
        "annonars/PAH_variant.json",
    ],
)
def test_annonars_variant_response_model(json_file):
    """Test AnnonarsVariantResponse model for various example responses."""
    annonars_response = get_json_object(json_file)
    assert AnnonarsVariantResponse.model_validate(annonars_response)
