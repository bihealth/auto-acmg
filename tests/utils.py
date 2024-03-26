"""Helper functions for the tests."""

import json
import os
from typing import Any, Dict


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
