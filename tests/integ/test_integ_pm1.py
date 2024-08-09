"""Integration tests for `PM1` criteria using upstream server."""

import pytest

from src.auto_acmg import AutoACMG
from src.core.config import Config
from src.criteria.auto_pm1 import AutoPM1
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGPrediction, AutoACMGResult
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


@pytest.mark.default_cassette("integ_pm1.yaml")
@pytest.mark.vcr
@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction",
    [
        # ("NM_000540.3(RYR1):c.1589G>A", GenomeRelease.GRCh37, True),
        ("NM_004700.3(KCNQ4):c.825G>C", GenomeRelease.GRCh37, AutoACMGPrediction.Met),
        # ("NM_000257.3(MYH7):c.1157A>G", GenomeRelease.GRCh37, True),
        ("NM_005343.3(HRAS):c.175G>A", GenomeRelease.GRCh37, AutoACMGPrediction.Met),
        # ("NM_030662.3(MAP2K2):c.400T>C", GenomeRelease.GRCh37, True),
        ("NM_004333.4(BRAF):c.739T>G", GenomeRelease.GRCh37, AutoACMGPrediction.Met),
        ("NM_001754.4(RUNX1):c.316T>A", GenomeRelease.GRCh37, AutoACMGPrediction.Met),
        ("NM_001754.4(RUNX1):c.314A>C", GenomeRelease.GRCh37, AutoACMGPrediction.Met),
        ("NM_001754.4:c.315C>A", GenomeRelease.GRCh37, AutoACMGPrediction.Met),
        ("NM_002834.4(PTPN11):c.782T>A", GenomeRelease.GRCh37, AutoACMGPrediction.Met),
        # ("NM_004333.6(BRAF):c.793G>C", GenomeRelease.GRCh37, True),
        # ("NM_005633.3(SOS1):c.1276C>A", GenomeRelease.GRCh37, True),
        ("NM_000546.5(TP53):c.396G>C", GenomeRelease.GRCh37, AutoACMGPrediction.Met),
        ("NM_005343.4(HRAS):c.175_176delinsCT", GenomeRelease.GRCh37, AutoACMGPrediction.Met),
        ("NM_001754.4(RUNX1):c.485G>A", GenomeRelease.GRCh37, AutoACMGPrediction.Met),
    ],
)
def test_pm1(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: bool,
    config: Config,
):
    # First, resolve variant
    auto_acmg = AutoACMG(variant_name, genome_release, config=config)
    seqvar = auto_acmg.resolve_variant()
    assert isinstance(seqvar, SeqVar)
    # Then, setup the data
    auto_acmg_result = auto_acmg.parse_data(seqvar)
    assert isinstance(auto_acmg_result, AutoACMGResult)
    # Then, predict PM1
    auto_pm1 = AutoPM1()
    pm1 = auto_pm1.predict_pm1(seqvar, auto_acmg_result.data)
    assert pm1.name == "PM1"
    assert pm1.prediction == expected_prediction
