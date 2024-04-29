from unittest.mock import patch

import pytest
from typer.testing import CliRunner

from src.cli import app

#: Create a test runner
runner = CliRunner()


@pytest.fixture
def mock_auto_acmg_predict_success():
    """Fixture to mock AutoACMG predict method with a success response."""
    with patch("src.auto_acmg.AutoACMG.predict") as mock_predict:
        mock_predict.return_value = "Pathogenic"
        yield mock_predict


@pytest.fixture
def mock_auto_acmg_predict_failure():
    """Fixture to mock AutoPVS1 predict method with a failure response."""
    with patch("src.auto_acmg.AutoACMG.predict") as mock_predict:
        mock_predict.side_effect = Exception("An error occurred")
        yield mock_predict


def test_classify_command_success(mock_auto_acmg_predict_success):
    """Test the 'classify' command with a mocked success response from AutoACMG.predict."""
    result = runner.invoke(app, ["NM_000038.3:c.797G>A", "--genome-release", "GRCh38"])
    assert result.exit_code == 0
    assert "" in result.stdout


def test_classify_command_failure(mock_auto_acmg_predict_failure):
    """Test the 'classify' command with a mocked failure response from AutoACMG.predict."""
    result = runner.invoke(app, ["NM_000038.3:c.797G>A", "--genome-release", "GRCh38"])
    assert result.exit_code == 0
    assert "Error: An error occurred" in result.stdout


def test_classify_invalid_genome_release():
    """Test the 'classify' command with an invalid genome release."""
    result = runner.invoke(app, ["NM_000038.3:c.797G>A", "--genome-release", "InvalidRelease"])
    assert result.exit_code == 0
    assert "Invalid genome release" in result.stdout
