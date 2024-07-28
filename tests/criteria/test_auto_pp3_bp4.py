from typing import Optional, Union
from unittest.mock import patch

import pytest

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.criteria.auto_pp3_bp4 import AutoPP3BP4
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.auto_acmg import PP3BP4
from src.defs.exceptions import AlgorithmError, MissingDataError
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
def auto_pp3bp4(seqvar, variant_info):
    return AutoPP3BP4(seqvar, variant_info.result)


@pytest.mark.parametrize(
    "score_value, expected",
    [
        (None, None),
        (5.0, 5.0),
        ("3.5", 3.5),
        ("3.5;4.2;2.1", 4.2),
        ("0.5;.;1.2", 1.2),
        ("-1.0;-0.5;-0.8", -0.5),
    ],
)
def test_convert_score_val(score_value, expected, auto_pp3bp4):
    """Test _convert_score_val method."""
    result = auto_pp3bp4._convert_score_val(score_value)
    assert result == expected, f"Expected {expected}, but got {result}"


def test_convert_score_val_invalid_string(auto_pp3bp4):
    """Test _convert_score_val method with invalid string input."""
    with pytest.raises(AlgorithmError):
        auto_pp3bp4._convert_score_val("invalid")


def test_auto_pp3bp4_initialization(seqvar, variant_info):
    """Test the initialization of AutoPP3BP4 class."""
    auto_pp3bp4 = AutoPP3BP4(seqvar=seqvar, variant_info=variant_info.result)
    assert auto_pp3bp4.seqvar == seqvar
    assert auto_pp3bp4.variant_info == variant_info.result
    assert isinstance(auto_pp3bp4.config, Config)
    assert isinstance(auto_pp3bp4.annonars_client, AnnonarsClient)
    assert auto_pp3bp4.prediction is None
    assert auto_pp3bp4.comment == ""


@pytest.mark.parametrize(
    "score_value, expected",
    [
        (None, None),
        ("3.5;.;1.2", 3.5),
    ],
)
def test_convert_score_val_various(score_value, expected, auto_pp3bp4):
    """Test _convert_score_val method with various valid and invalid inputs."""
    result = auto_pp3bp4._convert_score_val(score_value)
    assert result == expected


@pytest.mark.parametrize(
    "file_path, expected",
    [
        ("annonars/PAH_variant.json", True),
    ],
)
def test_is_pathogenic_score(file_path, expected, seqvar, variant_info):
    """Test is_pathogenic_score."""
    variant_response = AnnonarsVariantResponse.model_validate(get_json_object(file_path))
    auto_pp3bp4 = AutoPP3BP4(seqvar=seqvar, variant_info=variant_info.result)
    is_pathogenic = auto_pp3bp4._is_pathogenic_score(variant_response.result)
    assert is_pathogenic == expected


@pytest.mark.parametrize(
    "file_path",
    [
        ("annonars/RP11_variant.json"),
        ("annonars/LARP7_variant.json"),
    ],
)
def test_is_pathogenic_score_missing_data(file_path, seqvar, variant_info):
    """Test is_pathogenic_score with missing data."""
    variant_response = AnnonarsVariantResponse.model_validate(get_json_object(file_path))
    auto_pp3bp4 = AutoPP3BP4(seqvar=seqvar, variant_info=variant_info.result)
    with pytest.raises(MissingDataError):
        auto_pp3bp4._is_pathogenic_score(variant_response.result)


@pytest.mark.parametrize(
    "file_path, expected",
    [
        ("annonars/PAH_variant.json", False),
    ],
)
def test_is_benign_score(file_path, expected, seqvar, variant_info):
    """Test is_benign_score."""
    variant_response = AnnonarsVariantResponse.model_validate(get_json_object(file_path))
    auto_pp3bp4 = AutoPP3BP4(seqvar=seqvar, variant_info=variant_info.result)
    is_benign = auto_pp3bp4._is_benign_score(variant_response.result)
    assert is_benign == expected


@pytest.mark.parametrize(
    "file_path",
    [
        ("annonars/RP11_variant.json"),
        ("annonars/LARP7_variant.json"),
    ],
)
def test_is_benign_score_missing_data(file_path, seqvar, variant_info):
    """Test is_benign_score with missing data."""
    variant_response = AnnonarsVariantResponse.model_validate(get_json_object(file_path))
    auto_pp3bp4 = AutoPP3BP4(seqvar=seqvar, variant_info=variant_info.result)
    with pytest.raises(MissingDataError):
        auto_pp3bp4._is_benign_score(variant_response.result)


@pytest.mark.parametrize(
    "file_path, expected",
    [
        ("annonars/PAH_variant.json", False),
    ],
)
def test_is_pathogenic_spliceai(file_path, expected, seqvar, variant_info):
    """Test is_pathogenic_spliceai."""
    variant_response = AnnonarsVariantResponse.model_validate(get_json_object(file_path))
    auto_pp3bp4 = AutoPP3BP4(seqvar=seqvar, variant_info=variant_info.result)
    is_pathogenic = auto_pp3bp4._is_pathogenic_spliceai(variant_response.result)
    assert is_pathogenic == expected


@pytest.mark.parametrize(
    "file_path",
    [
        ("annonars/RP11_variant.json"),
        ("annonars/LARP7_variant.json"),
    ],
)
def test_is_pathogenic_spliceai_missing_data(file_path, seqvar, variant_info):
    """Test is_pathogenic_spliceai with missing data."""
    variant_response = AnnonarsVariantResponse.model_validate(get_json_object(file_path))
    auto_pp3bp4 = AutoPP3BP4(seqvar=seqvar, variant_info=variant_info.result)
    with pytest.raises(MissingDataError):
        auto_pp3bp4._is_pathogenic_spliceai(variant_response.result)


@pytest.mark.parametrize(
    "file_path, expected",
    [
        ("annonars/PAH_variant.json", True),
    ],
)
def test_is_benign_spliceai(file_path, expected, seqvar, variant_info):
    """Test is_benign_spliceai."""
    variant_response = AnnonarsVariantResponse.model_validate(get_json_object(file_path))
    auto_pp3bp4 = AutoPP3BP4(seqvar=seqvar, variant_info=variant_info.result)
    is_benign = auto_pp3bp4._is_benign_spliceai(variant_response.result)
    assert is_benign == expected


@pytest.mark.parametrize(
    "file_path",
    [
        ("annonars/RP11_variant.json"),
        ("annonars/LARP7_variant.json"),
    ],
)
def test_is_benign_spliceai_missing_data(file_path, seqvar, variant_info):
    """Test is_benign_spliceai with missing data."""
    variant_response = AnnonarsVariantResponse.model_validate(get_json_object(file_path))
    auto_pp3bp4 = AutoPP3BP4(seqvar=seqvar, variant_info=variant_info.result)
    with pytest.raises(MissingDataError):
        auto_pp3bp4._is_benign_spliceai(variant_response.result)
