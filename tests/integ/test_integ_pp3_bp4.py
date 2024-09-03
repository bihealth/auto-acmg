"""Integration tests for `PP3` and `BP4` critria using upstream server."""

from typing import Tuple

import pytest

from src.auto_acmg import VCEP_MAPPING, AutoACMG
from src.core.config import Config
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGSeqVarResult
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.auto_pp3_bp4 import AutoPP3BP4


@pytest.mark.default_cassette("integ_pp3_bp4.yaml")
@pytest.mark.vcr
@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction",
    [
        # ("NM_000277.2(PAH):c.472C>T", GenomeRelease.GRCh38, (True, False)),
        # ("NM_000277.2(PAH):c.533A>G", GenomeRelease.GRCh38, (True, False)),
        # ("NM_000277.2(PAH):c.194T>C", GenomeRelease.GRCh38, (True, False)),
        # ("NM_004999.3(MYO6):c.1025C>T", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000441.1(SLC26A4):c.1069G>A", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000527.5(LDLR):c.970G>A", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000441.1(SLC26A4):c.1363A>T", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000545.8(HNF1A):c.1135C>G", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000545.8(HNF1A):c.1136C>A", GenomeRelease.GRCh37, (True, False)),
        # ("NM_001110792.2(MECP2):c.968C>T", GenomeRelease.GRCh37, (True, False)),
        # ("NM_004985.4(KRAS):c.15A>T", GenomeRelease.GRCh37, (True, False)),
        # ("NM_004333.4(BRAF):c.1785T>G", GenomeRelease.GRCh37, (True, False)),
        # ("NM_004985.4(KRAS):c.101C>T", GenomeRelease.GRCh37, (True, False)),
        # ("NM_004985.4(KRAS):c.173C>T", GenomeRelease.GRCh37, (True, False)),
        # ("NM_176795.4(HRAS):c.173C>T", GenomeRelease.GRCh37, (True, False)),
        # ("NM_005343.3(HRAS):c.350A>G", GenomeRelease.GRCh37, (True, False)),
        # ("NM_005633.3(SOS1):c.1642A>C", GenomeRelease.GRCh37, (True, False)),
        # ("NM_005633.3(SOS1):c.806T>C", GenomeRelease.GRCh37, (True, False)),
        # ("NM_002880.3(RAF1):c.1472C>T", GenomeRelease.GRCh37, (True, False)),
        # ("NM_002755.3(MAP2K1):c.389A>G", GenomeRelease.GRCh37, (True, False)),
        # ("NM_002755.3(MAP2K1):c.275T>G", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000277.2(PAH):c.529G>A", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000277.2(PAH):c.848T>A", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000277.2(PAH):c.799C>G", GenomeRelease.GRCh37, (True, False)),
        # ("NM_004086.2(COCH):c.151C>T", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000277.2(PAH):c.1278T>C", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000277.2(PAH):c.707-7A>T", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000277.2(PAH):c.735G>A", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000257.3(MYH7):c.2162+4G>A", GenomeRelease.GRCh37, (False, True)),
        # ("NM_005633.3(SOS1):c.1230G>A", GenomeRelease.GRCh37, (False, True)),
        # ("NM_007373.3(SHOC2):c.1302C>T", GenomeRelease.GRCh37, (False, True)),
        # ("NM_002834.4(PTPN11):c.1658C>T", GenomeRelease.GRCh37, (False, True)),
        # ("NM_004333.4(BRAF):c.1227A>G", GenomeRelease.GRCh37, (False, True)),
        # ("NM_001754.4(RUNX1):c.614-34C>T", GenomeRelease.GRCh37, (False, True)),
        # ("NM_001754.4(RUNX1):c.613+8C>T", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000441.2(SLC26A4):c.416-7T>C", GenomeRelease.GRCh37, (False, True)),
        # ("NM_005633.3(SOS1):c.*4C>T", GenomeRelease.GRCh37, (False, True)),
        # ("NM_004985.5(KRAS):c.-160A>G", GenomeRelease.GRCh37, (False, True)),
        # ("NM_001754.4(RUNX1):c.1086G>C", GenomeRelease.GRCh37, (False, True)),
        # ("NM_001754.4(RUNX1):c.423G>A", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000546.5(TP53):c.704A>G", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000546.5(TP53):c.935C>G", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000546.5(TP53):c.869G>A", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000441.2(SLC26A4):c.1708-18T>A", GenomeRelease.GRCh37, (False, True)),
        # ("NM_206933.3(USH2A):c.12575G>A", GenomeRelease.GRCh37, (False, True)),
        # ("NM_001754.4(RUNX1):c.1396A>T", GenomeRelease.GRCh37, (False, True)),
        # ("NM_001754.4(RUNX1):c.952T>G", GenomeRelease.GRCh37, (False, True)),
        # ("NM_001754.4(RUNX1):c.183G>A", GenomeRelease.GRCh37, (False, True)),
        # ("NM_001754.4(RUNX1):c.144C>T", GenomeRelease.GRCh37, (False, True)),
        # ("NM_000212.2(ITGB3):c.342T>C", GenomeRelease.GRCh37, (False, True)),
    ],
)
def test_pp3_bp4(
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
    auto_acmg_result = auto_acmg.parse_seqvar_data(seqvar)
    assert isinstance(auto_acmg_result, AutoACMGSeqVarResult)
    # Then, predict PP3 and BP4
    if auto_acmg_result.data.hgnc_id in VCEP_MAPPING:
        predictor_class = VCEP_MAPPING[auto_acmg_result.data.hgnc_id]
        predictor = predictor_class(seqvar, auto_acmg_result, config)
        pp3_bp4 = predictor.predict_pp3bp4(seqvar, auto_acmg_result.data)
    else:
        auto_pp3_bp4 = AutoPP3BP4()
        pp3_bp4 = auto_pp3_bp4.predict_pp3bp4(seqvar, auto_acmg_result.data)
    assert pp3_bp4[0].name == "PP3"
    assert pp3_bp4[1].name == "BP4"
