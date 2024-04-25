"""Entry point for the command line interface."""

import typer
from typing_extensions import Annotated

from src.auto_acmg import AutoACMG
from src.defs.genome_builds import GenomeRelease
from src.pvs1.auto_pvs1 import AutoPVS1

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
):
    """
    Classify sequence variant on the ACMG guidelines.
    """
    try:
        genome_release_enum = GenomeRelease.from_string(genome_release)
        if not genome_release_enum:
            raise ValueError(
                (
                    f"Invalid genome release: {genome_release}. "
                    f"Please use one of {', '.join(ALLOWED_GENOME_RELEASES)}."
                )
            )

        auto_acmg = AutoACMG(variant, genome_release_enum)
        auto_acmg.predict()
    except Exception as e:
        typer.secho(f"Error: {e}", err=True, fg=typer.colors.RED)


if __name__ == "__main__":
    app()
