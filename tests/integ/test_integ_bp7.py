"""Integration tests for `BP7` criteria using upstream server."""

import pytest

from src.auto_acmg import VCEP_MAPPING, AutoACMG
from src.core.config import Config
from src.defs.auto_acmg import AutoACMGPrediction, AutoACMGSeqVarResult
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.auto_bp7 import AutoBP7


@pytest.mark.default_cassette("integ_bp7.yaml")
@pytest.mark.vcr
@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction",
    [
        # ("NM_000277.1:c.772C>T", GenomeRelease.GRCh37, True),
        # ("NM_000314.6(PTEN):c.254-39G>T", GenomeRelease.GRCh37, True),
        # ("NM_000314.6(PTEN):c.1104T>C", GenomeRelease.GRCh37, True),
        # ("NM_000314.6(PTEN):c.360A>C", GenomeRelease.GRCh37, True),
        # ("NM_000314.6(PTEN):c.18A>G", GenomeRelease.GRCh37, True),
        (
            "NM_000257.3(MYH7):c.3036C>T",
            GenomeRelease.GRCh37,
            AutoACMGPrediction.Applicable,
        ),
        (
            "NM_000156.6(GAMT):c.279C>T",
            GenomeRelease.GRCh37,
            AutoACMGPrediction.Applicable,
        ),
        (
            "NM_004700.3(KCNQ4):c.720C>G",
            GenomeRelease.GRCh37,
            AutoACMGPrediction.Applicable,
        ),
        (
            "NM_004985.4(KRAS):c.451-14T>C",
            GenomeRelease.GRCh37,
            AutoACMGPrediction.Applicable,
        ),
        (
            "NM_005343.3(HRAS):c.510G>A",
            GenomeRelease.GRCh37,
            AutoACMGPrediction.Applicable,
        ),
        (
            "NM_004333.4(BRAF):c.111G>A",
            GenomeRelease.GRCh37,
            AutoACMGPrediction.Applicable,
        ),
        (
            "NM_000152.5(GAA):c.1332T>C",
            GenomeRelease.GRCh37,
            AutoACMGPrediction.Applicable,
        ),
        (
            "NM_001754.4(RUNX1):c.1317C>T",
            GenomeRelease.GRCh37,
            AutoACMGPrediction.Applicable,
        ),
        (
            "NM_001754.4(RUNX1):c.36G>A",
            GenomeRelease.GRCh37,
            AutoACMGPrediction.Applicable,
        ),
        (
            "NM_001754.4(RUNX1):c.843C>T",
            GenomeRelease.GRCh37,
            AutoACMGPrediction.Applicable,
        ),
        (
            "NM_005633.3(SOS1):c.1230G>A",
            GenomeRelease.GRCh37,
            AutoACMGPrediction.Applicable,
        ),
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
    # Then, setup the data
    auto_acmg_result = auto_acmg._parse_seqvar_data(seqvar)
    assert isinstance(auto_acmg_result, AutoACMGSeqVarResult)
    # Then, predict BP7
    if auto_acmg_result.data.hgnc_id in VCEP_MAPPING:
        predictor_class = VCEP_MAPPING[auto_acmg_result.data.hgnc_id]
        predictor = predictor_class(seqvar, auto_acmg_result, config)
        bp7 = predictor.predict_bp7(seqvar, auto_acmg_result.data)
    else:
        auto_bp7 = AutoBP7()
        bp7 = auto_bp7.predict_bp7(seqvar, auto_acmg_result.data)
    assert bp7.name == "BP7"
    assert bp7.prediction == expected_prediction
