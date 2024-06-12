"""Integration tests for `BP7` criteria using upstream server."""

import pytest

from src.auto_acmg import AutoACMG
from src.core.config import Config
from src.criteria.auto_bp7 import AutoBP7
from src.criteria.auto_criteria import AutoACMGCriteria
from src.defs.annonars_variant import AnnonarsVariantResponse
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


@pytest.mark.remote
@pytest.mark.parametrize(
    "variant_name, genome_release, expected_prediction",
    [
        (
            "NM_000257.3(MYH7):c.3036C>T",
            GenomeRelease.GRCh37,
            True,
        ),
        (
            "NM_000156.6(GAMT):c.279C>T",
            GenomeRelease.GRCh37,
            True,
        ),
        (
            "NM_000277.2(PAH):c.963C>T",
            GenomeRelease.GRCh37,
            False,
        ), ### Found pathogenic variant in 2 bp proximity
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
    auto_acmg_criteria = AutoACMGCriteria(seqvar, genome_release, config=config)
    variant_info = auto_acmg_criteria._get_variant_info(seqvar)
    assert isinstance(variant_info, AnnonarsVariantResponse)
    # Then, predict BP7
    auto_bp7 = AutoBP7(seqvar, genome_release, variant_info.result, config=config)
    prediction, details = auto_bp7.predict()
    print(details)
    if expected_prediction is None:
        assert prediction is None, f"Failed for {variant_name}"
    else:
        assert prediction
        assert prediction.BP7 == expected_prediction, f"Failed for {variant_name}"
