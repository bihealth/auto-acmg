"""Implementations of the PVS1 algorithm."""

import typer

from src.genome_builds import GenomeRelease
from src.seqvar import SeqVar, SeqVarResolver
from src.seqvar_pvs1 import SeqVarPVS1
from src.types.autopvs1 import PVS1Prediction, PVS1PredictionSeqVarPath


class AutoPVS1:
    """AutoPVS1 algorithm for PVS1 criteria prediction."""

    def __init__(self, variant_name: str, genome_release: GenomeRelease = GenomeRelease.GRCh38):
        self.variant_name = variant_name
        self.genome_release = genome_release

    def resolve_variant(self) -> SeqVar | None:
        """Resolve the variant."""
        # TODO: Add resolve for Structure variants
        try:
            seqvar_resolver = SeqVarResolver()
            seqvar: SeqVar = seqvar_resolver.resolve_seqvar(self.variant_name, self.genome_release)
            typer.secho(f"Resolved variant: {seqvar}.", fg=typer.colors.BLUE)
            return seqvar
        except Exception as e:
            typer.secho(e, err=True, fg=typer.colors.RED)
            return None

    def predict(self):
        """Run the AutoPVS1 algorithm."""
        typer.secho(f"Running AutoPVS1 for variant {self.variant_name}.", fg=typer.colors.BLUE)
        variant = self.resolve_variant()

        if isinstance(variant, SeqVar):
            self.seqvar: SeqVar = variant
            self.prediction: PVS1Prediction = PVS1Prediction.NotPVS1
            self.prediction_path: PVS1PredictionSeqVarPath = PVS1PredictionSeqVarPath.NotSet

            try:
                typer.secho(
                    f"Predicting PVS1 for variant {self.seqvar.user_representation}, genome release: {self.genome_release.name}.",
                    fg=typer.colors.BLUE,
                )
                seqvar_pvs1 = SeqVarPVS1(self.seqvar)
                seqvar_pvs1.initialize()
                seqvar_pvs1.verify_PVS1()
                self.prediction, self.prediction_path = seqvar_pvs1.get_prediction()
                typer.secho(
                    f"PVS1 prediction for {self.seqvar.user_representation}: {self.prediction.name}",
                    fg=typer.colors.GREEN,
                )
            except Exception as e:
                typer.secho(
                    f"Failed to predict PVS1 for variant {self.seqvar.user_representation}.",
                    err=True,
                    fg=typer.colors.RED,
                )
                typer.secho(e, err=True)
                return

        elif isinstance(variant, str):
            # TODO: Add Structure variants PVS1 prediction
            pass
        else:
            typer.secho(
                f"Failed to resolve variant {self.variant_name}.", err=True, fg=typer.colors.RED
            )
            return
