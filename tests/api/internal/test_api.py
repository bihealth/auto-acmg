from unittest.mock import patch

import pytest
from fastapi.testclient import TestClient

from src.auto_acmg import AutoACMG
from src.core.config import settings
from src.defs.auto_acmg import AutoACMGSeqVarResult, AutoACMGStrucVarResult
from tests.utils import get_json_object

# ------------------- resolve_variant -------------------


async def test_resolve_variant_valid_data(client: TestClient):
    """Test resolving a variant with valid data."""
    # act:
    response = client.get(
        f"{settings.API_V1_STR}/resolve",
        params={"variant_name": "chr1:123456:G:A", "genome_release": "GRCh38"},
    )
    # assert:
    assert response.status_code == 200
    assert "variant_type" in response.json()


async def test_resolve_variant_invalid_genome_release(client: TestClient):
    """Test resolving a variant with an invalid genome release."""
    # act:
    response = client.get(
        f"{settings.API_V1_STR}/resolve",
        params={"variant_name": "BRCA1:c.5075G>A", "genome_release": "GRCh37"},
    )
    # assert:
    assert response.status_code == 400
    assert response.json()["detail"] == "Failed to resolve the variant"


async def test_resolve_variant_not_found(client: TestClient):
    """Test resolving a variant that cannot be found."""
    # act:
    response = client.get(
        f"{settings.API_V1_STR}/resolve",
        params={"variant_name": "BRCA1:c.9999G>A", "genome_release": "GRCh38"},
    )
    # assert:
    assert response.status_code == 400
    assert response.json()["detail"] == "Failed to resolve the variant"


# ------------------- predict_seqvar -------------------


@pytest.mark.asyncio
async def test_predict_seqvar_success(client: TestClient):
    """Test predicting a sequence variant successfully."""
    # Arrange
    variant_name = "chr1:228282272:G:A"
    mock_response = AutoACMGSeqVarResult.model_validate(
        get_json_object("var_data/example_seqvar_pred.json")
    )

    # Act
    with patch.object(AutoACMG, "predict", return_value=mock_response):
        response = client.get(
            f"{settings.API_V1_STR}/predict/seqvar",
            params={"variant_name": variant_name},
        )

    # Assert
    assert response.status_code == 200
    result = response.json()
    assert "seqvar" in result["prediction"]
    assert result["prediction"]["seqvar"]["user_repr"] == variant_name
    assert "criteria" in result["prediction"]


@pytest.mark.asyncio
async def test_predict_seqvar_invalid_genome_release(client: TestClient):
    """Test predicting a sequence variant with an invalid genome release."""
    # act:
    response = client.get(
        f"{settings.API_V1_STR}/predict/seqvar",
        params={"variant_name": "chr1:123456:G:A", "genome_release": "InvalidRelease"},
    )
    # assert:
    assert response.status_code == 400
    assert "Invalid genome release" in response.json()["detail"]


@pytest.mark.asyncio
async def test_predict_seqvar_invalid_variant(client: TestClient):
    """Test predicting an invalid sequence variant."""
    # act:
    response = client.get(
        f"{settings.API_V1_STR}/predict/seqvar",
        params={"variant_name": "invalid_variant", "genome_release": "GRCh38"},
    )
    # assert:
    assert response.status_code == 400
    assert "No valid sequence variant prediction was made" in response.json()["detail"]


@pytest.mark.asyncio
async def test_predict_seqvar_missing_variant_name(client: TestClient):
    """Test predicting a sequence variant without providing a variant name."""
    # act:
    response = client.get(
        f"{settings.API_V1_STR}/predict/seqvar",
        params={"genome_release": "GRCh38"},
    )
    # assert:
    assert response.status_code == 422  # Unprocessable Entity
    assert "variant_name" in response.json()["detail"][0]["loc"]


@pytest.mark.asyncio
async def test_predict_seqvar_structural_variant(client: TestClient):
    """Test predicting a structural variant using the sequence variant endpoint."""
    # act:
    response = client.get(
        f"{settings.API_V1_STR}/predict/seqvar",
        params={"variant_name": "1-100000-200000-DEL", "genome_release": "GRCh38"},
    )
    # assert:
    assert response.status_code == 400
    assert "No valid sequence variant prediction was made" in response.json()["detail"]


# ------------------- predict_strucvar -------------------


@pytest.mark.asyncio
async def test_predict_strucvar_success(client: TestClient):
    """Test predicting a structural variant successfully."""
    # Arrange
    variant_name = "DEL:chr17:41176312:41277500"
    mock_response = AutoACMGStrucVarResult.model_validate(
        get_json_object("var_data/example_strucvar_pred.json")
    )

    # Act
    with patch.object(AutoACMG, "predict", return_value=mock_response):
        response = client.get(
            f"{settings.API_V1_STR}/predict/strucvar",
            params={"variant_name": variant_name},
        )

    # Assert
    assert response.status_code == 200
    result = response.json()
    assert "strucvar" in result["prediction"]
    assert result["prediction"]["strucvar"]["user_repr"] == variant_name
    assert "criteria" in result["prediction"]


@pytest.mark.asyncio
async def test_predict_strucvar_invalid_genome_release(client: TestClient):
    """Test predicting a structural variant with an invalid genome release."""
    # Act
    response = client.get(
        f"{settings.API_V1_STR}/predict/strucvar",
        params={"variant_name": "chr1:g.1000000_2000000del", "genome_release": "InvalidRelease"},
    )
    # Assert
    assert response.status_code == 400
    assert "Invalid genome release" in response.json()["detail"]


@pytest.mark.asyncio
async def test_predict_strucvar_invalid_variant(client: TestClient):
    """Test predicting an invalid structural variant."""
    # Act
    response = client.get(
        f"{settings.API_V1_STR}/predict/strucvar",
        params={"variant_name": "invalid_variant", "genome_release": "GRCh38"},
    )
    # Assert
    assert response.status_code == 400
    assert "No valid structural variant prediction was made" in response.json()["detail"]


@pytest.mark.asyncio
async def test_predict_strucvar_missing_variant_name(client: TestClient):
    """Test predicting a structural variant without providing a variant name."""
    # Act
    response = client.get(
        f"{settings.API_V1_STR}/predict/strucvar",
        params={"genome_release": "GRCh38"},
    )
    # Assert
    assert response.status_code == 422  # Unprocessable Entity
    assert "variant_name" in response.json()["detail"][0]["loc"]


@pytest.mark.asyncio
async def test_predict_strucvar_sequence_variant(client: TestClient):
    """Test predicting a sequence variant using the structural variant endpoint."""
    # Act
    response = client.get(
        f"{settings.API_V1_STR}/predict/strucvar",
        params={"variant_name": "chr1:228282272:G:A", "genome_release": "GRCh38"},
    )
    # Assert
    assert response.status_code == 400
    assert "No valid structural variant prediction was made" in response.json()["detail"]


# ... existing tests ...
