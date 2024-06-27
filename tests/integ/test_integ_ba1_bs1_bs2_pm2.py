"""Integration tests for `BA1`, `BS1`, `BS2`, `PM2` criteria using upstream server."""

from typing import Tuple

import pytest

from src.auto_acmg import AutoACMG
from src.core.config import Config
from src.criteria.auto_ba1_bs1_bs2_pm2 import AutoBA1BS1BS2PM2
from src.criteria.auto_criteria import AutoACMGCriteria
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction",
    [
        # ("NM_005249.5(FOXG1):c.489C>T", GenomeRelease.GRCh37, [True, False, False, False]),
        # ("NM_000257.3(MYH7):c.2358G>T", GenomeRelease.GRCh37, [True, False, False, False]),
        # ("NM_000277.2(PAH):c.1278T>C", GenomeRelease.GRCh37, [True, False, False, False]),
        # ("NM_000277.2(PAH):c.707-7A>T", GenomeRelease.GRCh37, [True, False, False, False]),
        # ("NM_000277.2(PAH):c.1242C>T", GenomeRelease.GRCh37, [False, True, True, False]),
        # ("NM_000277.2(PAH):c.1A>G", GenomeRelease.GRCh37, [False, False, False, True]),
        ("NM_000277.2(PAH):c.503delA", GenomeRelease.GRCh37, [False, False, False, True]),
        ("NM_000277.1(PAH):c.814G>T", GenomeRelease.GRCh37, [False, False, False, True]),
        ("NM_000277.1(PAH):c.1162G>A", GenomeRelease.GRCh37, [False, False, False, True]),
        ("NM_000277.1(PAH):c.722G>A", GenomeRelease.GRCh37, [False, False, False, True]),
    ],
)
def test_ba1_bs1_bs2_pm2(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: Tuple[bool, bool, bool, bool],
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
    # Then, predict BA1, BS1, BS2, PM2
    auto_ba1_bs1_bs2_pm2 = AutoBA1BS1BS2PM2(seqvar, variant_info.result, config=config)
    prediction, details = auto_ba1_bs1_bs2_pm2.predict()
    print(details)
    if expected_prediction is None:
        assert prediction is None, f"Failed for {variant_name}"
    else:
        assert prediction
        assert prediction.BA1 == expected_prediction[0], f"Failed for {variant_name}"
        assert prediction.BS1 == expected_prediction[1], f"Failed for {variant_name}"
        assert prediction.BS2 == expected_prediction[2], f"Failed for {variant_name}"
        assert prediction.PM2 == expected_prediction[3], f"Failed for {variant_name}"
