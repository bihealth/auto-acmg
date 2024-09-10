"""Entry point for the command line interface."""

import typer
from loguru import logger
from typing_extensions import Annotated

from src.auto_acmg import AutoACMG
from src.core.config import settings
from src.defs.exceptions import AutoAcmgBaseException, InvalidGenomeBuild
from src.defs.genome_builds import GenomeRelease

app = typer.Typer()

#: Allowed genome releases
ALLOWED_GENOME_RELEASES = ["GRCh37", "GRCh38", "hg19", "hg38", "grch37", "grch38"]
#: Allowed sequence variant formats
ALLOWED_SEQVAR_FORMATS = ["Canonical SPDI", "gnomAD", "relaxed SPDI", "dbSNP", "ClinVar"]
#: Allowed structural variant formats
ALLOWED_STRUCVAR_FORMATS = ["Colon-separated", "Hyphen-separated"]


@app.command()
def classify(
    variant: Annotated[
        str,
        typer.Argument(
            help=(
                f"Variant to be classified, e.g., 'NM_000038.3:c.797G>A'. "
                f"Accepted sequence variants formats: {', '.join(ALLOWED_SEQVAR_FORMATS)}. "
                f"Accepted structural variants formats: {', '.join(ALLOWED_STRUCVAR_FORMATS)}."
            )
        ),
    ],
    genome_release: Annotated[
        str,
        typer.Option(
            "--genome-release",
            "-g",
            help=f"Accepted genome Releases: {', '.join(ALLOWED_GENOME_RELEASES)}",
        ),
    ] = "GRCh38",
    duplication_tandem: Annotated[
        bool,
        typer.Option(
            "--duplication-tandem",
            "-dt",
            help=(
                "Flag to indicate if the duplication is in tandem AND disrupts reading frame AND "
                "undergoes NMD."
            ),
        ),
    ] = False,
):
    """
    Classify sequence variant on the ACMG guidelines.
    """
    try:
        genome_release_enum = GenomeRelease.from_string(genome_release)
        if not genome_release_enum:
            logger.error(
                (
                    f"Invalid genome release: {genome_release}. "
                    f"Please use one of {', '.join(ALLOWED_GENOME_RELEASES)}."
                )
            )
            raise InvalidGenomeBuild("Invalid genome release")
        # Temporary save the duplication tandem flag in the settings
        settings.DUPLICATION_TANDEM = duplication_tandem
        auto_acmg = AutoACMG(variant, genome_release_enum)
        prediction = auto_acmg.predict()
        prediction.save_to_file() if prediction else logger.error("No prediction was made.")
    except AutoAcmgBaseException as e:
        logger.error("Error occurred: {}", e)


if __name__ == "__main__":
    app()
