import csv
import os

import genebe as gnb  # type: ignore

from src.core.config import settings

#: Path to the root directory
path_to_root = settings.PATH_TO_ROOT

#: Genebe API key
genebe_api_key = settings.GENEBE_API_KEY

#: Genebe username
genebe_username = settings.GENEBE_USERNAME


# =========== Prepare the data ===========
# Extract criteria from the csv file
path = os.path.join(path_to_root, "src", "bench", "comparison_criteria_custom.csv")
print(f"Data path: {path}")
variants = []
with open(path, "rt") as inputf:
    reader = csv.DictReader(inputf)
    for record in reader:
        if record["variant"].startswith("#"):
            continue
        variants.append(
            (record["variant"], record["expected_criteria"].split(";"), record["comment"])
        )
print(f"Number of variants: {len(variants)}")
# ========================================

# =========== Run the benchmark ===========
list_of_variants = gnb.parse_hgvs(
    [variant[0] for variant in variants], username=genebe_username, api_key=genebe_api_key
)
df = gnb.annotate_variants_list_to_dataframe(
    list_of_variants, flatten_consequences=True, username=genebe_username, api_key=genebe_api_key
)

# Export the dataframe to a csv file
df.to_csv("output.csv", index=False)
