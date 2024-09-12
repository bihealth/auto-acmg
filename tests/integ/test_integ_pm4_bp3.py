"""Integration tests for `PM4` and `BP3` critria using upstream server."""

from typing import Tuple

import pytest

from src.auto_acmg import VCEP_MAPPING, AutoACMG
from src.core.config import Config
from src.defs.auto_acmg import AutoACMGPrediction, AutoACMGSeqVarResult
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.seqvar.auto_pm4_bp3 import AutoPM4BP3


@pytest.mark.default_cassette("integ_pm4_bp3.yaml")
@pytest.mark.vcr
@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction",
    [
        # ("NM_000277.2(PAH):c.169_171delGAG", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000277.3(PAH):c.124_126del", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000419.4:c.480C>G", GenomeRelease.GRCh37, (True, False)),
        # ("NM_001110792.2(MECP2):c.1497A>G", GenomeRelease.GRCh37, (True, False)),
        # ("NM_000277.1(PAH):c.967_969delACA", GenomeRelease.GRCh37, (True, False)),
        # ("NM_206933.2(USH2A):c.8559-2A>G", GenomeRelease.GRCh37, (True, False)),
        (
            "NM_000277.2(PAH):c.1092_1094del",
            GenomeRelease.GRCh37,
            (AutoACMGPrediction.Applicable, AutoACMGPrediction.NotApplicable),
        ),
        (
            "NM_000277.2(PAH):c.1092_1106del",
            GenomeRelease.GRCh37,
            (AutoACMGPrediction.Applicable, AutoACMGPrediction.NotApplicable),
        ),
        (
            "NM_000277.1(PAH):c.208_210delTCT",
            GenomeRelease.GRCh37,
            (AutoACMGPrediction.Applicable, AutoACMGPrediction.NotApplicable),
        ),
        # ("NM_005249.4(FOXG1):c.209_232del24", GenomeRelease.GRCh37, (False, True)),
        # ("NM_005249.5(FOXG1):c.209_235del", GenomeRelease.GRCh37, (False, True)),
        # ("NM_005249.5(FOXG1):c.209_235dup", GenomeRelease.GRCh37, (False, True)),
    ],
)
def test_pm4_bp3(
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
    # Then, predict PM4 and BP3
    if auto_acmg_result.data.hgnc_id in VCEP_MAPPING:
        predictor_class = VCEP_MAPPING[auto_acmg_result.data.hgnc_id]
        predictor = predictor_class(seqvar, auto_acmg_result, config)
        pm4_bp3 = predictor.predict_pm4bp3(seqvar, auto_acmg_result.data)
    else:
        auto_pm4_bp3 = AutoPM4BP3()
        pm4_bp3 = auto_pm4_bp3.predict_pm4bp3(seqvar, auto_acmg_result.data)
    assert pm4_bp3[0].name == "PM4"
    assert pm4_bp3[0].prediction == expected_prediction[0]
    assert pm4_bp3[1].name == "BP3"
    assert pm4_bp3[1].prediction == expected_prediction[1]
