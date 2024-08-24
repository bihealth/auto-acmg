import io
from typing import Type
from unittest.mock import MagicMock, patch

import pytest
from loguru import logger
from typer.testing import CliRunner

from src.cli import app
from src.defs.exceptions import AutoAcmgBaseException

#: Create a test runner
runner = CliRunner()


@pytest.fixture
def loguru_capture():
    buffer = io.StringIO()
    # Add a sink that logs to a StringIO object
    logger_id = logger.add(buffer, format="{message}")

    yield buffer

    # Remove the sink after the test to clean up
    logger.remove(logger_id)


@pytest.fixture
def mock_auto_acmg_predict_success():
    """Fixture to mock AutoACMG predict method with a success response."""
    with patch("src.auto_acmg.AutoACMG.predict") as mock_predict:
        mock_predict.return_value = MagicMock(save_to_file=MagicMock())
        yield mock_predict


@pytest.fixture
def mock_auto_acmg_predict_failure(request: pytest.FixtureRequest):
    """Fixture to mock AutoPVS1 predict method with a failure response."""
    exception: Type[Exception] = request.param
    with patch("src.auto_acmg.AutoACMG.predict") as mock_predict:
        mock_predict.side_effect = exception("An error occurred")
        yield mock_predict


def test_classify_success(mock_auto_acmg_predict_success, loguru_capture):
    """Test the 'classify' command with a mocked success response from AutoACMG.predict."""
    result = runner.invoke(app, ["NM_000038.3:c.797G>A", "--genome-release", "GRCh38"])
    log_output = loguru_capture.getvalue()
    assert result.exit_code == 0
    assert "" in log_output


@pytest.mark.parametrize("mock_auto_acmg_predict_failure", [AutoAcmgBaseException], indirect=True)
def test_classify_auto_acmg_failure(mock_auto_acmg_predict_failure, loguru_capture):
    """Test the 'classify' command with a mocked failure response from AutoACMG.predict using AutoAcmgBaseException."""
    result = runner.invoke(app, ["NM_000038.3:c.797G>A", "--genome-release", "GRCh38"])
    log_output = loguru_capture.getvalue()
    assert result.exit_code == 0
    assert "Error occurred" in log_output


@pytest.mark.parametrize("mock_auto_acmg_predict_failure", [ValueError], indirect=True)
def test_classify_generic_exception(mock_auto_acmg_predict_failure, loguru_capture):
    """Test the 'classify' command with a generic exception using ValueError."""
    result = runner.invoke(app, ["NM_000038.3:c.797G>A", "--genome-release", "GRCh38"])
    log_output = loguru_capture.getvalue()
    assert result.exit_code == 1
    assert "Error occurred" not in log_output


def test_classify_invalid_genome_release(loguru_capture):
    """Test the 'classify' command with an invalid genome release."""
    result = runner.invoke(app, ["NM_000038.3:c.797G>A", "--genome-release", "InvalidRelease"])
    log_output = loguru_capture.getvalue()
    assert result.exit_code == 0
    assert "Invalid genome release: InvalidRelease" in log_output
