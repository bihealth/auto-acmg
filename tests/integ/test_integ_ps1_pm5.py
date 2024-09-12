"""Integration tests for `PS1` and `PM5` critria using upstream server."""

from typing import Tuple

import pytest

from src.auto_acmg import VCEP_MAPPING, AutoACMG
from src.core.config import Config
from src.defs.auto_acmg import AutoACMGPrediction, AutoACMGSeqVarResult
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.auto_ps1_pm5 import AutoPS1PM5


@pytest.mark.default_cassette("integ_ps1_pm5.yaml")
@pytest.mark.vcr
@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction",
    [
        # ("NM_176795.4(HRAS):c.173C>T", GenomeRelease.GRCh38, (True, False)),
        # ("NM_004333.4(BRAF):c.739T>C", GenomeRelease.GRCh38, (True, False)),
        # ("NM_000527.5(LDLR):c.1A>T", GenomeRelease.GRCh38, (True, False)),
        (
            "NM_000277.2(PAH):c.865G>C",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.Applicable, AutoACMGPrediction.NotApplicable),
        ),
        (
            "NM_004333.4(BRAF):c.735A>T",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.Applicable, AutoACMGPrediction.NotApplicable),
        ),
        (
            "NM_005633.3(SOS1):c.1656G>T",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.Applicable, AutoACMGPrediction.NotApplicable),
        ),
        (
            "NM_000277.3(PAH):c.865G>A",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.Applicable, AutoACMGPrediction.NotApplicable),
        ),
        (
            "NM_000540.2(RYR1):c.1021G>C",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.Applicable, AutoACMGPrediction.NotApplicable),
        ),
        (
            "NM_000540.2(RYR1):c.742G>C",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.Applicable, AutoACMGPrediction.NotApplicable),
        ),
        (
            "NM_001110792.2(MECP2):c.408G>C",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.Applicable, AutoACMGPrediction.NotApplicable),
        ),
        (
            "NM_001110792.2(MECP2):c.408G>T",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.Applicable, AutoACMGPrediction.NotApplicable),
        ),
        (
            "NM_000162.5(GCK):c.523G>A",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.Applicable, AutoACMGPrediction.NotApplicable),
        ),
        # ("NM_001110792.2(MECP2):c.507C>G", GenomeRelease.GRCh38, (True, True)),
        # ("NM_000314.6(PTEN):c.388C>G", GenomeRelease.GRCh38, (False, True)),
        # ("NM_000260.3(MYO7A):c.3503G>A", GenomeRelease.GRCh38, (False, True)),
        # ("NM_000152.3(GAA):c.1933G>A", GenomeRelease.GRCh38, (False, True)),
        # ("NM_000277.1(PAH):c.1241A>G", GenomeRelease.GRCh38, (False, True)),
        # ("NM_000277.1(PAH):c.722G>A", GenomeRelease.GRCh38, (False, True)),
        # ("NM_000257.3(MYH7):c.2155C>T", GenomeRelease.GRCh38, (False, True)),
        # ("NM_002834.4(PTPN11):c.1530G>C", GenomeRelease.GRCh38, (False, True)),
        # ("NM_004985.4(KRAS):c.458A>T", GenomeRelease.GRCh38, (False, True)),
        (
            "NM_000277.1(PAH):c.782G>C",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_000277.2(PAH):c.734T>C",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_000257.3(MYH7):c.2221G>T",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_000257.3(MYH7):c.2221G>C",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_000257.3(MYH7):c.2167C>T",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_000257.3(MYH7):c.2167C>G",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_000257.3(MYH7):c.2156G>A",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_000257.3(MYH7):c.1208G>A",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_000257.3(MYH7):c.1207C>T",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_000257.3(MYH7):c.1358G>A",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_002880.3(RAF1):c.775T>A",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_004985.4(KRAS):c.101C>T",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_002755.3(MAP2K1):c.388T>C",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_000540.3(RYR1):c.1840C>T",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_000540.3(RYR1):c.6487C>T",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_003159.2(CDKL5):c.215T>C",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_001323289.1:c.59G>T",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_014297.5(ETHE1):c.488G>A",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_000257.4(MYH7):c.2710C>T",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_004958.4(MTOR):c.4448G>A",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
        (
            "NM_004958.3:c.5930C>G",
            GenomeRelease.GRCh38,
            (AutoACMGPrediction.NotApplicable, AutoACMGPrediction.Applicable),
        ),
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
    # Then, setup the data
    auto_acmg_result = auto_acmg._parse_seqvar_data(seqvar)
    assert isinstance(auto_acmg_result, AutoACMGSeqVarResult)
    # Then, predict PS1 and PM5
    if auto_acmg_result.data.hgnc_id in VCEP_MAPPING:
        predictor_class = VCEP_MAPPING[auto_acmg_result.data.hgnc_id]
        predictor = predictor_class(seqvar, auto_acmg_result, config)
        ps1_pm5 = predictor.predict_ps1pm5(seqvar, auto_acmg_result.data)
    else:
        auto_ps1_pm5 = AutoPS1PM5()
        ps1_pm5 = auto_ps1_pm5.predict_ps1pm5(seqvar, auto_acmg_result.data)
    assert ps1_pm5[0].name == "PS1"
    assert ps1_pm5[1].name == "PM5"
