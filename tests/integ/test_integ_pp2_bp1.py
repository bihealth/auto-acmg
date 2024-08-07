"""Integration tests for `PP2` and `BP1` critria using upstream server."""

from typing import Tuple

import pytest

from src.auto_acmg import AutoACMG
from src.core.config import Config
from src.criteria.auto_pp2_bp1 import AutoPP2BP1
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGResult
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


@pytest.mark.default_cassette("integ_pp2_bp1.yaml")
@pytest.mark.vcr
@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction",
    [
        # ("NM_002880.3(RAF1):c.769T>C", GenomeRelease.GRCh37, (True, False)),
        # ("NM_005343.3(HRAS):c.175G>A", GenomeRelease.GRCh37, (True, True)),
        # ("NM_005633.3(SOS1):c.1018C>T", GenomeRelease.GRCh37, (True, False)),
        # ("NM_002880.3(RAF1):c.94A>G", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000314.6(PTEN):c.520T>A", GenomeRelease.GRCh37, (True, False)),
        # ("NM_004333.4(BRAF):c.1595G>A", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000314.6(PTEN):c.209T>C", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000314.6(PTEN):c.170T>G", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000314.6(PTEN):c.331T>C", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000314.6(PTEN):c.181C>G", GenomeRelease.GRCh37, (True, False)),
        # ("NM_004333.4(BRAF):c.739T>G", GenomeRelease.GRCh37, (True, False)),
        # ("NM_002834.4(PTPN11):c.167T>C", GenomeRelease.GRCh37, (True, False)),
        # ("NM_005633.3(SOS1):c.512T>C", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000314.6(PTEN):c.737C>T", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000314.6(PTEN):c.510T>A", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000314.6(PTEN):c.368A>G", GenomeRelease.GRCh37, (True, False)),
        # ("NM_005633.3(SOS1):c.322G>A", GenomeRelease.GRCh37, (True, False)),
        # ("NM_007373.3(SHOC2):c.4A>G", GenomeRelease.GRCh37, (True, False)),
        # ("NM_002880.3(RAF1):c.775T>A", GenomeRelease.GRCh37, (True, False)),
        # ("NM_002880.3(RAF1):c.1082G>C", GenomeRelease.GRCh37, (True, False)),
        # ("NM_030662.3(MAP2K2):c.619G>A", GenomeRelease.GRCh37, (True, False)),
        # ("NM_004985.4(KRAS):c.15A>T", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000038.6(APC):c.295C>T", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000038.6(APC):c.715G>C", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000038.6(APC):c.995G>A", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000038.6(APC):c.1240C>T", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000038.6(APC):c.2476T>G", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000038.6(APC):c.4399C>T", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000038.6(APC):c.4420G>A", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000038.6(APC):c.5465T>A", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000038.6(APC):c.7862C>G", GenomeRelease.GRCh37, (False, True)),
        # ("NM_024675.3(PALB2):c.721A>G", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000038.6(APC):c.235A>G", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000038.6(APC):c.420G>C", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000038.6(APC):c.1487C>T", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000038.6(APC):c.4906G>T", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000038.6(APC):c.8438C>A", GenomeRelease.GRCh37, (False, True)),
    ],
)
def test_pp2_bp1(
    variant_name: str,
    genome_release: GenomeRelease,
    expected_prediction: Tuple[bool, bool],
    config: Config,
):
    # First, resolve variant
    auto_acmg = AutoACMG(variant_name, genome_release, config=config)
    seqvar = auto_acmg.resolve_variant()
    assert isinstance(seqvar, SeqVar)
    # Then, setup the data
    auto_acmg_result = auto_acmg.parse_data(seqvar)
    assert isinstance(auto_acmg_result, AutoACMGResult)
    # Then, predict PP2 and BP1
    auto_pp2_bp1 = AutoPP2BP1()
    pp2_bp1 = auto_pp2_bp1.predict_pp2bp1(seqvar, auto_acmg_result.data)
    assert pp2_bp1[0].name == "PP2"
    assert pp2_bp1[1].name == "BP1"
