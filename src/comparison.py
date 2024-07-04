import csv
import os
import time

import pandas as pd
import requests

from src.auto_acmg import AutoACMG, AutoACMGPrediction
from src.core.config import settings
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar

# Path to the root directory
path_to_root = settings.PATH_TO_ROOT

# Extract criteria from the csv file
path = os.path.join(path_to_root, "tests", "assets", "e2e_variants", "comparison_criteria.csv")
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


def eval_autoacmg(pred, expected):
    crit_met = []
    for crit in pred.model_dump().values():
        if crit["prediction"] == AutoACMGPrediction.Positive:
            crit_met.append(crit["name"])
    tp = list(set(expected) & set(crit_met))
    fn = list(set(expected) - set(crit_met))
    fp = list(set(crit_met) - set(expected))
    return crit_met, tp, fn, fp


def eval_intervar(pred, expected):
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


def intervar_response(variant: str):
    """
    Implement searching for ACMG classification for SNVs and indels.
    Proxy requests to the `WinterVar <http://wintervar.wglab.org/>`_ backend.

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
    backend_resp = requests.get(url)
    backend_resp.raise_for_status()
    return backend_resp.json()


# Run AutoACMG predictions


def append_row(df, row):
    return pd.concat([df, pd.DataFrame([row], columns=row.index)]).reset_index(drop=True)


# Create a pandas DataFrame
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
    ]
)

for var in variants:
    record = {
        "Variant": var[0],
        "Expected Criteria": ";".join(var[1]),
        "Comment": var[2],
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
        "AutoACMG Full Response": "",
        "Intervar Full Response": "",
    }

    # AutoACMG
    try:
        start_time = time.time()
        auto_acmg = AutoACMG(var[0], GenomeRelease.GRCh37)
        pred = auto_acmg.predict()
        end_time = time.time()
        # Evaluate the model
        crit_met, tp, fn, fp = eval_autoacmg(pred, var[1])
        record["AutoACMG Criteria"]: ";".join(crit_met)  # type: ignore
        record["AutoACMG Prediction time"]: end_time - start_time  # type: ignore
        record["AutoACMG True Positives"]: ";".join(tp)  # type: ignore
        record["AutoACMG False Negatives"]: ";".join(fn)  # type: ignore
        record["AutoACMG False Positives"]: ";".join(fp)  # type: ignore
        record["AutoACMG Full Response"]: pred.model_dump()  # type: ignore
    except Exception as e:
        print(f"Exception was raised for {var[0]} in AutoACMG:\n{e}")

    # Intervar
    try:
        start_time = time.time()
        resp = intervar_response(var[0])
        end_time = time.time()
        crit_met, tp, fn, fp = eval_intervar(resp, var[1])
        record["Intervar Criteria"]: ";".join(crit_met)  # type: ignore
        record["Intervar Prediction time"]: end_time - start_time  # type: ignore
        record["Intervar True Positives"]: ";".join(tp)  # type: ignore
        record["Intervar False Negatives"]: ";".join(fn)  # type: ignore
        record["Intervar False Positives"]: ";".join(fp)  # type: ignore
        record["Intervar Full Response"]: resp  # type: ignore
    except Exception as e:
        print(f"Exception was raised for {var[0]} in InterVar:\n{e}")

    stats = append_row(stats, pd.Series(record))

output_path = os.path.join(path_to_root, "src", "jupyter", "stats.csv")
stats.to_csv(output_path, index=False)
