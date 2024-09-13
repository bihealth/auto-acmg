"""Integration tests for `PM2`, `BA1`, `BS1`, `BS2` criteria using upstream server."""

import pytest

from src.auto_acmg import VCEP_MAPPING, AutoACMG
from src.core.config import Config
from src.defs.auto_acmg import AutoACMGSeqVarResult
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.auto_pm2_ba1_bs1_bs2 import AutoPM2BA1BS1BS2


@pytest.mark.default_cassette("integ_pm2_ba1_bs1_bs2.yaml")
@pytest.mark.vcr
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
        # ("NM_000277.2(PAH):c.503delA", GenomeRelease.GRCh37, [True, False, False, False]),
        # ("NM_000277.1(PAH):c.814G>T", GenomeRelease.GRCh37, [True, False, False, False]),
        # ("NM_000277.1(PAH):c.1162G>A", GenomeRelease.GRCh37, [True, False, False, False]),
        # ("NM_000277.1(PAH):c.722G>A", GenomeRelease.GRCh37, [False, False, False, True]),
        # ("NM_000277.2(PAH):c.735G>A", GenomeRelease.GRCh37, (True, False, True, False)),
        # ("NM_000314.6(PTEN):c.-1311T>C", GenomeRelease.GRCh37, (True, False, False, False)),
        # ("NM_000314.6(PTEN):c.-903G>A", GenomeRelease.GRCh37, (True, False, False, False)),
        # ("NM_000314.6(PTEN):c.-9C>G", GenomeRelease.GRCh37, (True, False, False, False)),
        # ("NM_206933.2(USH2A):c.15562A>G", GenomeRelease.GRCh37, (True, False, False, False)),
        # ("NM_004999.3(MYO6):c.475G>A", GenomeRelease.GRCh37, (True, False, False, False)),
        # ("NM_000260.3(MYO7A):c.324C>T", GenomeRelease.GRCh37, (True, False, False, False)),
        # ("NM_000257.3(MYH7):c.3036C>T", GenomeRelease.GRCh37, (True, False, False, False)),
        # ("NM_005633.3(SOS1):c.1018C>T", GenomeRelease.GRCh37, (False, True, False, False)),
        # ("NM_002880.3(RAF1):c.94A>G", GenomeRelease.GRCh37, (False, True, True, False)),
        # ("NM_001754.4(RUNX1):c.423G>A", GenomeRelease.GRCh37, (False, True, False, False)),
        # ("NM_000546.5(TP53):c.704A>G", GenomeRelease.GRCh37, (False, True, False, False)),
        # ("NM_000546.5(TP53):c.935C>G", GenomeRelease.GRCh37, (False, True, False, False)),
        # ("NM_000546.5(TP53):c.869G>A", GenomeRelease.GRCh37, (False, True, False, False)),
        # ("NM_001754.4(RUNX1):c.179C>T", GenomeRelease.GRCh37, (False, True, False, False)),
        # ("NM_001754.4(RUNX1):c.205G>C", GenomeRelease.GRCh37, (False, True, False, False)),
        # ("NM_000540.2(RYR1):c.12884C>T", GenomeRelease.GRCh37, (False, True, False, False)),
        # ("NM_000540.3(RYR1):c.6961A>G", GenomeRelease.GRCh37, (False, True, False, False)),
        # ("NM_000284.4(PDHA1):c.999A>C", GenomeRelease.GRCh37, (False, True, False, False)),
        # ("NM_004992.3(MECP2):c.1140G>A", GenomeRelease.GRCh37, (False, True, True, False)),
        # ("NM_004992.3(MECP2):c.1030C>T", GenomeRelease.GRCh37, (False, True, True, False)),
        # ("NM_000527.4(LDLR):c.1323C>T", GenomeRelease.GRCh37, (False, True, True, False)),
        # ("NM_000277.3(PAH):c.399T>C", GenomeRelease.GRCh37, (False, True, True, False)),
        # ("NM_001110792.2(MECP2):c.1375G>A", GenomeRelease.GRCh37, (False, True, True, False)),
        # ("NM_001110792.2(MECP2):c.1199C>T", GenomeRelease.GRCh37, (False, True, True, False)),
        # ("NM_001110792.2(MECP2):c.1196C>T", GenomeRelease.GRCh37, (False, True, True, False)),
        # ("NM_001110792.2(MECP2):c.915C>G", GenomeRelease.GRCh37, (False, True, True, False)),
        # ("NM_000051.3(ATM):c.1073A>G", GenomeRelease.GRCh37, (False, True, False, False)),
        # ("NM_000051.4(ATM):c.1066-6T>G", GenomeRelease.GRCh37, (False, True, False, False)),
        # ("NM_000051.4(ATM):c.3925G>A", GenomeRelease.GRCh37, (False, True, False, False)),
        # ("NC_012920.1:m.8557G>A", GenomeRelease.GRCh37, (False, True, False, False)),
        # ("NM_000545.8(HNF1A):c.1849G>A", GenomeRelease.GRCh37, (False, True, True, False)),
        # ("NM_004086.2(COCH):c.151C>T", GenomeRelease.GRCh37, (False, False, False, True)),
        # ("NM_004700.3(KCNQ4):c.853G>A", GenomeRelease.GRCh37, (False, False, False, True)),
        # ("NM_206933.2(USH2A):c.4510dupA", GenomeRelease.GRCh37, (False, False, False, True)),
        # ("NM_000441.1(SLC26A4):c.365dupT", GenomeRelease.GRCh37, (False, False, False, True)),
        # ("NM_000277.2(PAH):c.350delC", GenomeRelease.GRCh37, (False, False, False, True)),
        # ("NM_000277.2(PAH):c.58C>T", GenomeRelease.GRCh37, (False, False, False, True)),
        # ("NM_000277.2(PAH):c.498C>A", GenomeRelease.GRCh37, (False, False, False, True)),
        # ("NM_000277.2(PAH):c.521T>C", GenomeRelease.GRCh37, (False, False, False, True)),
        # ("NM_000277.2(PAH):c.504C>A", GenomeRelease.GRCh37, (False, False, False, True)),
    ],
)
def test_pm2ba1bs1bs2(
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
    auto_acmg_result = auto_acmg._parse_seqvar_data(seqvar)
    assert isinstance(auto_acmg_result, AutoACMGSeqVarResult)
    # Then, predict PM2, BA1, BS1, BS2
    if auto_acmg_result.data.hgnc_id in VCEP_MAPPING:
        predictor_class = VCEP_MAPPING[auto_acmg_result.data.hgnc_id]
        predictor = predictor_class(seqvar, auto_acmg_result, config)
        pm2ba1bs1bs2 = predictor.predict_pm2ba1bs1bs2(seqvar, auto_acmg_result.data)
    else:
        auto_pm2ba1bs1bs2 = AutoPM2BA1BS1BS2()
        pm2ba1bs1bs2 = auto_pm2ba1bs1bs2.predict_pm2ba1bs1bs2(seqvar, auto_acmg_result.data)
    assert pm2ba1bs1bs2[0].name == "PM2"
    assert pm2ba1bs1bs2[1].name == "BA1"
    assert pm2ba1bs1bs2[2].name == "BS1"
    assert pm2ba1bs1bs2[3].name == "BS2"
