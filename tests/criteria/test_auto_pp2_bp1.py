from unittest.mock import MagicMock, patch

import pytest

from src.criteria.auto_pp2_bp1 import AutoPP2BP1
from src.defs.exceptions import AlgorithmError, InvalidAPIResposeError
from src.defs.seqvar import SeqVar


@pytest.fixture
def seqvar():
    return MagicMock(spec=SeqVar)


@pytest.fixture
def auto_pp2bp1():
    return AutoPP2BP1()


@pytest.mark.skip(reason="Annonars is not mocked properly")
@patch("src.utils.AutoACMGHelper.annonars_client.get_variant_from_range")
def test_get_missense_vars_success(mock_get_variant, auto_pp2bp1, seqvar):
    # Setup the mock response with pathogenic and benign variants
    pathogenic_variant = MagicMock()
    pathogenic_variant.records = [
        MagicMock(
            classifications=MagicMock(germlineClassification=MagicMock(description="Pathogenic"))
        )
    ]
    benign_variant = MagicMock()
    benign_variant.records = [
        MagicMock(classifications=MagicMock(germlineClassification=MagicMock(description="Benign")))
    ]
    response = MagicMock(clinvar=[pathogenic_variant, benign_variant])
    mock_get_variant.return_value = response

    pathogenic, benign, total = auto_pp2bp1._get_missense_vars(seqvar, 100, 200)

    assert pathogenic == 1
    assert benign == 1
    assert total == 2


@pytest.mark.skip(reason="Annonars is not mocked properly")
@patch("src.utils.AutoACMGHelper.annonars_client.get_variant_from_range")
def test_get_missense_vars_empty_response(mock_get_variant, auto_pp2bp1, seqvar):
    # Setup the mock to return an empty response
    mock_get_variant.return_value = MagicMock(clinvar=[])

    with pytest.raises(InvalidAPIResposeError):
        auto_pp2bp1._get_missense_vars(seqvar, 100, 200)


def test_get_missense_vars_invalid_range(auto_pp2bp1, seqvar):
    # Test error handling when the end position is less than the start position
    with pytest.raises(AlgorithmError) as excinfo:
        auto_pp2bp1._get_missense_vars(seqvar, 200, 100)
    assert "End position is less than the start position" in str(excinfo.value)
