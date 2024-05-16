"""Helper functions for the tests."""

import csv
import json
import os
from typing import Any, Dict, List, Tuple

from src.defs.auto_pvs1 import PVS1Prediction, PVS1PredictionSeqVarPath
from src.defs.genome_builds import GenomeRelease
from src.defs.mehari import Exon


def get_json_object(file_name: str) -> Dict[str, Any]:
    """
    Reads a JSON file from the assets folder and returns it as a dictionary.

    :param file_name: The name of the JSON file to read.
    :type file_name: str
    :return: A dictionary containing the JSON content.
    :rtype: Dict[str, Any]
    """
    file_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "assets", file_name))

    # Check if the file exists
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file {file_name} does not exist in the assets folder.")

    with open(file_path, "r") as file:
        json_object = json.load(file)

    return json_object


def load_test_data(
    path: str,
) -> List[Tuple[str, GenomeRelease, PVS1Prediction, PVS1PredictionSeqVarPath]]:
    """Load CSV test data/."""
    result = []
    with open(path, "rt") as inputf:
        reader = csv.DictReader(inputf)
        for record in reader:
            if record["section"].startswith("#"):
                continue
            result.append(
                (
                    record["variant_name"],
                    GenomeRelease[record["genome_release"]],
                    PVS1Prediction[record["expected_prediction"]],
                    PVS1PredictionSeqVarPath[record["expected_path"]],
                )
            )
    return result
