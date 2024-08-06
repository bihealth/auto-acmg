"""Integration tests for `PM4` and `BP3` critria using upstream server."""

from typing import Tuple

import pytest

from src.auto_acmg import AutoACMG
from src.core.config import Config
from src.criteria.auto_criteria import AutoACMGCriteria
from src.criteria.auto_pm4_bp3 import AutoPM4BP3
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


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
        ("NM_000277.2(PAH):c.1092_1094del", GenomeRelease.GRCh37, (True, False)),
        ("NM_000277.2(PAH):c.1092_1106del", GenomeRelease.GRCh37, (True, False)),
        ("NM_000277.1(PAH):c.208_210delTCT", GenomeRelease.GRCh37, (True, False)),
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
    # Then, fetch the variant_info from Annonars
    auto_acmg_criteria = AutoACMGCriteria(seqvar, config=config)
    variant_info = auto_acmg_criteria._get_variant_info(seqvar)
    assert isinstance(variant_info, AnnonarsVariantResponse)
    # Then, make prediction
    auto_pm4_bp3 = AutoPM4BP3(seqvar, variant_info.result, config=config)
    prediction, details = auto_pm4_bp3.predict()
    print(details)
    if expected_prediction is None:
        assert prediction is None, f"Failed for {variant_name}"
    else:
        assert prediction
        assert prediction.PM4 == expected_prediction[0]
        assert prediction.BP3 == expected_prediction[1]
