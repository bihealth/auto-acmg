"""Integration tests for `BP7` criteria using upstream server."""

import pytest

from src.auto_acmg import AutoACMG
from src.core.config import Config
from src.criteria.auto_criteria import AutoACMGCriteria
from src.criteria.auto_pm1 import AutoPM1
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction",
    [
        # ("NM_000540.3(RYR1):c.1589G>A", GenomeRelease.GRCh37, True),
        ("NM_004700.3(KCNQ4):c.825G>C", GenomeRelease.GRCh37, True),
        ("NM_000257.3(MYH7):c.1157A>G", GenomeRelease.GRCh37, True),
        ("NM_005343.3(HRAS):c.175G>A", GenomeRelease.GRCh37, True),
        ("NM_030662.3(MAP2K2):c.400T>C", GenomeRelease.GRCh37, True),
        ("NM_004333.4(BRAF):c.739T>G", GenomeRelease.GRCh37, True),
        ("NM_001754.4(RUNX1):c.316T>A", GenomeRelease.GRCh37, True),
        ("NM_001754.4(RUNX1):c.314A>C", GenomeRelease.GRCh37, True),
        ("NM_001754.4:c.315C>A", GenomeRelease.GRCh37, True),
        ("NM_002834.4(PTPN11):c.782T>A", GenomeRelease.GRCh37, True),
        ("NM_004333.6(BRAF):c.793G>C", GenomeRelease.GRCh37, True),
        ("NM_005633.3(SOS1):c.1276C>A", GenomeRelease.GRCh37, True),
        ("NM_000546.5(TP53):c.396G>C", GenomeRelease.GRCh37, True),
        ("NM_005343.4(HRAS):c.175_176delinsCT", GenomeRelease.GRCh37, True),
        ("NM_001754.4(RUNX1):c.485G>A", GenomeRelease.GRCh37, True),
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
    # Then, predict PM1
    auto_pm1 = AutoPM1(seqvar, variant_info.result, config=config)
    prediction, details = auto_pm1.predict()
    print(details)
    if expected_prediction is None:
        assert prediction is None, f"Failed for {variant_name}"
    else:
        assert prediction
        assert prediction.PM1 == expected_prediction, f"Failed for {variant_name}"
