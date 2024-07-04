"""Integration tests for `BP7` criteria using upstream server."""

import pytest

from src.auto_acmg import AutoACMG
from src.core.config import Config
from src.criteria.auto_bp7 import AutoBP7
from src.criteria.auto_criteria import AutoACMGCriteria
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction",
    [
        # ("NM_000277.1:c.772C>T", GenomeRelease.GRCh37, True),
        # ("NM_000314.6(PTEN):c.254-39G>T", GenomeRelease.GRCh37, True),
        # ("NM_000314.6(PTEN):c.1104T>C", GenomeRelease.GRCh37, True),
        # ("NM_000314.6(PTEN):c.360A>C", GenomeRelease.GRCh37, True),
        # ("NM_000314.6(PTEN):c.18A>G", GenomeRelease.GRCh37, True),
        ("NM_000257.3(MYH7):c.3036C>T", GenomeRelease.GRCh37, True),
        ("NM_000156.6(GAMT):c.279C>T", GenomeRelease.GRCh37, True),
        ("NM_004700.3(KCNQ4):c.720C>G", GenomeRelease.GRCh37, True),
        ("NM_004985.4(KRAS):c.451-14T>C", GenomeRelease.GRCh37, True),
        ("NM_005343.3(HRAS):c.510G>A", GenomeRelease.GRCh37, True),
        ("NM_004333.4(BRAF):c.111G>A", GenomeRelease.GRCh37, True),
        ("NM_000152.5(GAA):c.1332T>C", GenomeRelease.GRCh37, True),
        ("NM_001754.4(RUNX1):c.1317C>T", GenomeRelease.GRCh37, True),
        ("NM_001754.4(RUNX1):c.36G>A", GenomeRelease.GRCh37, True),
        ("NM_001754.4(RUNX1):c.843C>T", GenomeRelease.GRCh37, True),
        ("NM_005633.3(SOS1):c.1230G>A", GenomeRelease.GRCh37, True),
    ],
)
def test_bp7(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: bool,
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
    # Then, predict BP7
    auto_bp7 = AutoBP7(seqvar, variant_info.result, config=config)
    prediction, details = auto_bp7.predict()
    print(details)
    if expected_prediction is None:
        assert prediction is None, f"Failed for {variant_name}"
    else:
        assert prediction
        assert prediction.BP7 == expected_prediction, f"Failed for {variant_name}"
