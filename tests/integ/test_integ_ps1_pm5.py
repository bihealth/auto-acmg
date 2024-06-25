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
        # ("NM_004333.4(BRAF):c.739T>C", GenomeRelease.GRCh38, (True, False)),
        # ("NM_000527.5(LDLR):c.1A>T", GenomeRelease.GRCh38, (True, False)),
        ("NM_000277.2(PAH):c.865G>C", GenomeRelease.GRCh38, (True, False)),
        ("NM_004333.4(BRAF):c.735A>T", GenomeRelease.GRCh38, (True, False)),
        ("NM_005633.3(SOS1):c.1656G>T", GenomeRelease.GRCh38, (True, False)),
        ("NM_000277.3(PAH):c.865G>A", GenomeRelease.GRCh38, (True, False)),
        ("NM_000540.2(RYR1):c.1021G>C", GenomeRelease.GRCh38, (True, False)),
        ("NM_000540.2(RYR1):c.742G>C", GenomeRelease.GRCh38, (True, False)),
        ("NM_001110792.2(MECP2):c.408G>C", GenomeRelease.GRCh38, (True, False)),
        ("NM_001110792.2(MECP2):c.408G>T", GenomeRelease.GRCh38, (True, False)),
        ("NM_000162.5(GCK):c.523G>A", GenomeRelease.GRCh38, (True, False)),
        # ("NM_001110792.2(MECP2):c.507C>G", GenomeRelease.GRCh38, (True, True)),
        # ("NM_000314.6(PTEN):c.388C>G", GenomeRelease.GRCh38, (False, True)),
        # ("NM_000260.3(MYO7A):c.3503G>A", GenomeRelease.GRCh38, (False, True)),
        # ("NM_000152.3(GAA):c.1933G>A", GenomeRelease.GRCh38, (False, True)),
        # ("NM_000277.1(PAH):c.1241A>G", GenomeRelease.GRCh38, (False, True)),
        # ("NM_000277.1(PAH):c.722G>A", GenomeRelease.GRCh38, (False, True)),
        # ("NM_000257.3(MYH7):c.2155C>T", GenomeRelease.GRCh38, (False, True)),
        # ("NM_002834.4(PTPN11):c.1530G>C", GenomeRelease.GRCh38, (False, True)),
        # ("NM_004985.4(KRAS):c.458A>T", GenomeRelease.GRCh38, (False, True)),
        ("NM_000277.1(PAH):c.782G>C", GenomeRelease.GRCh38, (False, True)),
        ("NM_000277.2(PAH):c.734T>C", GenomeRelease.GRCh38, (False, True)),
        ("NM_000257.3(MYH7):c.2221G>T", GenomeRelease.GRCh38, (False, True)),
        ("NM_000257.3(MYH7):c.2221G>C", GenomeRelease.GRCh38, (False, True)),
        ("NM_000257.3(MYH7):c.2167C>T", GenomeRelease.GRCh38, (False, True)),
        ("NM_000257.3(MYH7):c.2167C>G", GenomeRelease.GRCh38, (False, True)),
        ("NM_000257.3(MYH7):c.2156G>A", GenomeRelease.GRCh38, (False, True)),
        ("NM_000257.3(MYH7):c.1208G>A", GenomeRelease.GRCh38, (False, True)),
        ("NM_000257.3(MYH7):c.1207C>T", GenomeRelease.GRCh38, (False, True)),
        ("NM_000257.3(MYH7):c.1358G>A", GenomeRelease.GRCh38, (False, True)),
        ("NM_002880.3(RAF1):c.775T>A", GenomeRelease.GRCh38, (False, True)),
        ("NM_004985.4(KRAS):c.101C>T", GenomeRelease.GRCh38, (False, True)),
        ("NM_002755.3(MAP2K1):c.388T>C", GenomeRelease.GRCh38, (False, True)),
        ("NM_000540.3(RYR1):c.1840C>T", GenomeRelease.GRCh38, (False, True)),
        ("NM_000540.3(RYR1):c.6487C>T", GenomeRelease.GRCh38, (False, True)),
        ("NM_003159.2(CDKL5):c.215T>C", GenomeRelease.GRCh38, (False, True)),
        ("NM_001323289.1:c.59G>T", GenomeRelease.GRCh38, (False, True)),
        ("NM_014297.5(ETHE1):c.488G>A", GenomeRelease.GRCh38, (False, True)),
        ("NM_000257.4(MYH7):c.2710C>T", GenomeRelease.GRCh38, (False, True)),
        ("NM_004958.4(MTOR):c.4448G>A", GenomeRelease.GRCh38, (False, True)),
        ("NM_004958.3:c.5930C>G", GenomeRelease.GRCh38, (False, True)),
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
