import pytest

from src.defs.annonars_range import AnnonarsRangeResponse
from tests.utils import get_json_object


@pytest.mark.parametrize(
    "json_file",
    [
        "annonars/annonars_range_success.json",
        "annonars/NM_000152.4:c.1A>G_annonars_range.json",
    ],
)
def test_annonars_range_response_model(json_file):
    """Test AnnonarsRangeResponse model for various example responses."""
    annonars_response = get_json_object(json_file)
    assert AnnonarsRangeResponse.model_validate(annonars_response)
