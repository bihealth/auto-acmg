import pytest

from src.defs.annonars_range import AnnonarsRangeResponse
from src.defs.annonars_variant import AnnonarsVariantResponse
from tests.utils import get_json_object


@pytest.mark.parametrize(
    "json_file",
    [
        "annonars/annonars_range_success.json",
        "annonars/annonars_range_failure.json",
        "annonars/GAA_range.json",
        "annonars/CDH1_range.json",
        "annonars/SHH_annonars_range_success.json",
    ],
)
def test_annonars_range_response_model(json_file):
    """Test AnnonarsRangeResponse model for various example responses."""
    annonars_response = get_json_object(json_file)
    assert AnnonarsRangeResponse.model_validate(annonars_response)


@pytest.mark.parametrize(
    "json_file",
    [
        "annonars/annonars_variant_success.json",
        "annonars/annonars_variant_failure.json",
        "annonars/larp7_variant.json",
        "annonars/PAH_annonars_variant.json",
    ],
)
def test_annonars_variant_response_model(json_file):
    """Test AnnonarsVariantResponse model for various example responses."""
    annonars_response = get_json_object(json_file)
    assert AnnonarsVariantResponse.model_validate(annonars_response)
