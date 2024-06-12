"""Integration tests for `PP3` and `BP4` critria using upstream server."""

from typing import Tuple

import pytest

from src.auto_acmg import AutoACMG
from src.core.config import Config
from src.criteria.auto_criteria import AutoACMGCriteria
from src.criteria.auto_pp3_bp4 import AutoPP3BP4
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction",
    [
        ("NM_000277.2(PAH):c.472C>T", GenomeRelease.GRCh38, (True, True)),
        ("NM_004985.4(KRAS):c.173C>T", GenomeRelease.GRCh38, (True, True)),
        ("NM_000277.1(PAH):c.721C>T", GenomeRelease.GRCh38, None),
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
    # Then, fetch the variant_info from Annonars
    auto_acmg_criteria = AutoACMGCriteria(seqvar, config=config)
    variant_info = auto_acmg_criteria._get_variant_info(seqvar)
    assert isinstance(variant_info, AnnonarsVariantResponse)
    # Then, make prediction
    auto_pp3_bp4 = AutoPP3BP4(seqvar, variant_info.result, config=config)
    prediction, details = auto_pp3_bp4.predict()
    print(details)
    if expected_prediction is None:
        assert prediction is None, f"Failed for {variant_name}"
    else:
        assert prediction
        assert prediction.PP3 == expected_prediction[0]
        assert prediction.BP4 == expected_prediction[1]
