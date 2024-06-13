"""Integration tests for `PS1` and `PM5` critria using upstream server."""

from typing import Tuple

import pytest

from src.auto_acmg import AutoACMG
from src.core.config import Config
from src.criteria.auto_criteria import AutoACMGCriteria
from src.criteria.auto_ps1_pm5 import AutoPS1PM5
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction",
    [
        # ("NM_176795.4(HRAS):c.173C>T", GenomeRelease.GRCh38, (True, False)),
        ("NM_000277.2(PAH):c.865G>C", GenomeRelease.GRCh38, (True, False)),
        ("NM_004333.4(BRAF):c.735A>T", GenomeRelease.GRCh38, (True, False)),
        # ("NM_004333.4(BRAF):c.739T>C", GenomeRelease.GRCh38, (True, False)),
        # ("NM_000314.6(PTEN):c.388C>G", GenomeRelease.GRCh38, (False, True)),
        ("NM_000277.1(PAH):c.782G>C", GenomeRelease.GRCh38, (False, True)),
        # ("NM_000260.3(MYO7A):c.3503G>A", GenomeRelease.GRCh38, (False, True)),
        # ("NM_000152.3(GAA):c.1933G>A", GenomeRelease.GRCh38, (False, True)),
    ],
)
def test_ps1_pm5(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: Tuple[bool, bool],
    config: Config,
):
    # First, resolve variant
    auto_acmg = AutoACMG(variant_name, genome_release, config=config)
    seqvar = auto_acmg.resolve_variant()
    assert isinstance(seqvar, SeqVar)
    # Then, fetch the variant_info from Annonars
    auto_acmg_criteria = AutoACMGCriteria(seqvar, config=config)
    variant_info = auto_acmg_criteria._get_variant_info(seqvar)
    assert isinstance(variant_info, AnnonarsVariantResponse)
    # Then, make prediction
    auto_ps1_pm5 = AutoPS1PM5(seqvar, variant_info.result, config=config)
    prediction, details = auto_ps1_pm5.predict()
    print(details)
    if expected_prediction is None:
        assert prediction is None, f"Failed for {variant_name}"
    else:
        assert prediction
        assert prediction.PS1 == expected_prediction[0]
        assert prediction.PM5 == expected_prediction[1]
