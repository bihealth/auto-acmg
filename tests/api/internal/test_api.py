import pytest
from fastapi.testclient import TestClient

from src.core.config import settings

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
async def test_predict_seqvar_valid_data(client: TestClient, mocker):
    """Test predicting a sequence variant with valid data."""
    # Arrange:
    mock_predict = mocker.patch("src.auto_acmg.AutoACMG.predict")
    mock_predict.return_value = {
        "seqvar": "chr1:123456:G:A",
        "data": {"some": "data"},
        "criteria": {"PVS1": "Present"},
    }

    # Act:
    response = client.get(
        f"{settings.API_V1_STR}/predict/seqvar",
        params={"variant_name": "chr1:123456:G:A", "genome_release": "GRCh38"},
    )

    # Assert:
    # assert response.status_code == 200
    # assert response.json() == {
    #     "prediction": {
    #         "seqvar": "chr1:123456:G:A",
    #         "data": {"some": "data"},
    #         "criteria": {"PVS1": "Present"},
    #     }
    # }
    mock_predict.assert_called_once_with("chr1:123456:G:A", "GRCh38")


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
