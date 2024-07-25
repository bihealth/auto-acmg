import csv
import json
import os
import time

import httpx
import pandas as pd

from src.auto_acmg import AutoACMG, AutoACMGPrediction
from src.core.config import settings
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar

#: Path to the root directory
path_to_root = settings.PATH_TO_ROOT

#: List of all criteria
criteria = [
    "PVS1",
    "PS1",
    "PS2",
    "PS3",
    "PS4",
    "PM1",
    "PM2",
    "PM3",
    "PM4",
    "PM5",
    "PM6",
    "PP1",
    "PP2",
    "PP3",
    "PP4",
    "PP5",
    "BA1",
    "BS1",
    "BS2",
    "BS3",
    "BS4",
    "BP1",
    "BP2",
    "BP3",
    "BP4",
    "BP5",
    "BP6",
    "BP7",
]

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


# =========== Helper functions ===========
#: Extract criteria from the line
def extract_criteria(line):
    ret = []
    for c in criteria:
        if c in line:
            ret.append(c)
    return ret


#: Append a row to the DataFrame
def append_row(df, row):
    return pd.concat([df, pd.DataFrame([row], columns=row.index)]).reset_index(drop=True)


def intervar_response(variant: str):
    """
    Implement searching for ACMG classification for SNVs and indels.
    Proxy httpx to the `WinterVar <http://wintervar.wglab.org/>`_ backend.

    :param variant: request
    :return: ACMG classification
    :rtype: dict
    """
    auto_acmg = AutoACMG(variant, GenomeRelease.GRCh37)
    seqvar = auto_acmg.resolve_variant()
    assert isinstance(seqvar, SeqVar)
    chromosome = seqvar.chrom
    position = seqvar.pos
    reference = seqvar.delete
    alternative = seqvar.insert

    if not chromosome or not position or not reference or not alternative:
        return

    url = (
        f"http://wintervar.wglab.org/api_new.php?"
        f"queryType=position&chr={chromosome}&pos={position}"
        f"&ref={reference}&alt={alternative}&build=hg19"
    )
    print("Requesting:", url)
    backend_resp = httpx.get(url)
    backend_resp.raise_for_status()
    return backend_resp.json()


def genebe_response(variant: str):
    """
    Implement searching for ACMG classification for SNVs and indels using GeneBe.

    :param variant: variant string
    :return: GeneBe response
    :rtype: dict
    """
    auto_acmg = AutoACMG(variant, GenomeRelease.GRCh37)
    seqvar = auto_acmg.resolve_variant()
    assert isinstance(seqvar, SeqVar)
    chromosome = seqvar.chrom
    position = seqvar.pos
    reference = seqvar.delete
    alternative = seqvar.insert

    if not chromosome or not position or not reference or not alternative:
        return

    url = (
        f"https://api.genebe.net/cloud/api-public/v1/variant?"
        f"chr={chromosome}&pos={position}"
        f"&ref={reference}&alt={alternative}&genome=hg19"
    )
    print("Requesting:", url)
    backend_resp = httpx.get(url)
    backend_resp.raise_for_status()
    return backend_resp.json()


def eval_autoacmg(pred, expected):
    """
    Evaluate the AutoACMG prediction.

    :param pred: AutoACMG prediction
    :param expected: expected criteria
    :return: criteria met, true positives, false negatives, false positives
    :rtype: tuple
    """
    crit_met = []
    for crit in pred.model_dump().values():
        if crit["prediction"] == AutoACMGPrediction.Met:
            crit_met.append(crit["name"])
    tp = list(set(expected) & set(crit_met))
    fn = list(set(expected) - set(crit_met))
    fp = list(set(crit_met) - set(expected))
    return crit_met, tp, fn, fp


def eval_intervar(pred, expected):
    """
    Evaluate the InterVar prediction.

    :param pred: InterVar prediction
    :param expected: expected criteria
    :return: criteria met, true positives, false negatives, false positives
    :rtype: tuple
    """
    crit_met = []
    for crit in pred:
        if crit in [
            "Intervar",
            "Build",
            "Chromosome",
            "Position",
            "Ref_allele",
            "Alt_allele",
            "Gene",
        ]:
            continue
        if pred[crit] == 1:
            crit_met.append(crit)
    tp = list(set(expected) & set(crit_met))
    fn = list(set(expected) - set(crit_met))
    fp = list(set(crit_met) - set(expected))
    return crit_met, tp, fn, fp


def eval_genebe(resp, expected):
    """
    Evaluate the VarSome prediction.

    :param resp: VarSome prediction
    :param expected: expected criteria
    :return: criteria met, true positives, false negatives, false positives
    :rtype: tuple
    """
    crit_met = []
    try:
        acmg_crit = resp["variants"][0]["acmg_criteria"]
        crit_met = extract_criteria(acmg_crit)
    except Exception as e:
        print(e)
    tp = list(set(expected) & set(crit_met))
    fn = list(set(expected) - set(crit_met))
    fp = list(set(crit_met) - set(expected))
    return crit_met, tp, fn, fp


# ========================================


# =========== Main comparison run ===========
#: Create a pandas DataFrame
stats = pd.DataFrame(
    columns=[
        "Variant",
        "Expected Criteria",
        "AutoACMG Criteria",
        "AutoACMG Prediction time",
        "AutoACMG True Positives",
        "AutoACMG False Negatives",
        "AutoACMG False Positives",
        "Intervar Criteria",
        "Intervar Prediction time",
        "Intervar True Positives",
        "Intervar False Negatives",
        "Intervar False Positives",
        "Genebe Criteria",
        "Genebe Prediction time",
        "Genebe True Positives",
        "Genebe False Negatives",
        "Genebe False Positives",
        "Comment",
        "AutoACMG Full Response",
        "Intervar Full Response",
        "Genebe Full Response",
    ]
)

for i, var in enumerate(variants):
    # Save the stats every 10 variants
    if i % 10 == 0 and i > 0:
        print(f"Processed {i} variants")
        output_path = os.path.join(path_to_root, "src", "bench", "tmp", f"_stats_{i}.csv")
        stats.to_csv(output_path, index=False)

    record = {
        "Variant": var[0],
        "Expected Criteria": ";".join(var[1]),
        "AutoACMG Criteria": "",
        "AutoACMG Prediction time": 0,
        "AutoACMG True Positives": "",
        "AutoACMG False Negatives": "",
        "AutoACMG False Positives": "",
        "Intervar Criteria": "",
        "Intervar Prediction time": 0,
        "Intervar True Positives": "",
        "Intervar False Negatives": "",
        "Intervar False Positives": "",
        "Genebe Criteria": "",
        "Genebe Prediction time": 0,
        "Genebe True Positives": "",
        "Genebe False Negatives": "",
        "Genebe False Positives": "",
        "Comment": var[2],
        "AutoACMG Full Response": "",
        "Intervar Full Response": "",
        "Genebe Full Response": "",
    }

    # AutoACMG
    try:
        start_time = time.time()
        auto_acmg = AutoACMG(var[0], GenomeRelease.GRCh37)
        pred = auto_acmg.predict()
        end_time = time.time()
        assert pred
        # Evaluate the model
        crit_met, tp, fn, fp = eval_autoacmg(pred, var[1])
        record["AutoACMG Criteria"] = ";".join(crit_met)
        record["AutoACMG Prediction time"] = end_time - start_time
        record["AutoACMG True Positives"] = ";".join(tp)
        record["AutoACMG False Negatives"] = ";".join(fn)
        record["AutoACMG False Positives"] = ";".join(fp)
        record["AutoACMG Full Response"] = pred.model_dump()
    except Exception as e:
        print(f"Exception was raised for {var[0]} in AutoACMG:\n{e}")

    # Intervar
    try:
        start_time = time.time()
        resp = intervar_response(var[0])
        end_time = time.time()
        crit_met, tp, fn, fp = eval_intervar(resp, var[1])
        record["Intervar Criteria"] = ";".join(crit_met)
        record["Intervar Prediction time"] = end_time - start_time
        record["Intervar True Positives"] = ";".join(tp)
        record["Intervar False Negatives"] = ";".join(fn)
        record["Intervar False Positives"] = ";".join(fp)
        record["Intervar Full Response"] = resp
    except Exception as e:
        print(f"Exception was raised for {var[0]} in InterVar:\n{e}")

    # Genebe
    try:
        start_time = time.time()
        resp = genebe_response(var[0])
        end_time = time.time()
        crit_met, tp, fn, fp = eval_genebe(resp, var[1])
        record["Genebe Criteria"] = ";".join(crit_met)
        record["Genebe Prediction time"] = end_time - start_time
        record["Genebe True Positives"] = ";".join(tp)
        record["Genebe False Negatives"] = ";".join(fn)
        record["Genebe False Positives"] = ";".join(fp)
        record["Genebe Full Response"] = json.dumps(resp)
    except Exception as e:
        print(f"Exception was raised for {var[0]} in Genebe:\n{e}")

    # Skip if both AutoACMG and Intervar did not return any criteria
    if (
        record["AutoACMG Criteria"] == ""
        and record["Intervar Criteria"] == ""
        and record["Genebe Criteria"] == ""
    ):
        print(f"Skipping {var[0]}")
        continue

    # Add the record to the stats
    print(f"Record for {var[0]}: {record}")
    stats = append_row(stats, pd.Series(record))

# Save the final stats
output_path = os.path.join(path_to_root, "src", "bench", "stats.csv")
stats.to_csv(output_path, index=False)
