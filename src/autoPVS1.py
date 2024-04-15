"""Implementations of the PVS1 algorithm."""

import typer

from src.defs.autopvs1 import (
    PVS1Prediction,
    PVS1PredictionSeqVarPath,
    PVS1PredictionStrucVarPath,
)
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar, SeqVarResolver
from src.defs.strucvar import StrucVar, StrucVarResolver
from src.seqvar_pvs1 import SeqVarPVS1
from src.strucvar_pvs1 import StrucVarPVS1


class AutoPVS1:
    """AutoPVS1 algorithm for PVS1 criteria prediction."""

    def __init__(self, variant_name: str, genome_release: GenomeRelease = GenomeRelease.GRCh38):
        self.variant_name = variant_name
        self.genome_release = genome_release

    def resolve_variant(self) -> SeqVar | StrucVar | None:
        """Resolve the variant."""
        try:
            try:
                seqvar_resolver = SeqVarResolver()
                seqvar: SeqVar = seqvar_resolver.resolve_seqvar(
                    self.variant_name, self.genome_release
                )
                typer.secho(f"Resolved variant: {seqvar}.", fg=typer.colors.BLUE)
                return seqvar
            except Exception as e:
                strucvar_resolver = StrucVarResolver()
                strucvar: StrucVar = strucvar_resolver.resolve_strucvar(
                    self.variant_name, self.genome_release
                )
                typer.secho(f"Resolved structural variant: {strucvar}.", fg=typer.colors.BLUE)
                return strucvar
        except Exception as e:
            typer.secho(e, err=True, fg=typer.colors.RED)
            return None

    def predict(self):
        """Run the AutoPVS1 algorithm."""
        typer.secho(f"Running AutoPVS1 for variant {self.variant_name}.", fg=typer.colors.BLUE)
        variant = self.resolve_variant()

        if isinstance(variant, SeqVar):
            self.seqvar: SeqVar = variant
            self.seqvar_prediction: PVS1Prediction = PVS1Prediction.NotPVS1
            self.seqvar_prediction_path: PVS1PredictionSeqVarPath = PVS1PredictionSeqVarPath.NotSet

            try:
                typer.secho(
                    f"Predicting PVS1 for variant {self.seqvar.user_repr}, genome release: {self.genome_release.name}.",
                    fg=typer.colors.BLUE,
                )
                seqvar_pvs1 = SeqVarPVS1(self.seqvar)
                seqvar_pvs1.initialize()
                seqvar_pvs1.verify_PVS1()
                self.seqvar_prediction, self.seqvar_prediction_path = seqvar_pvs1.get_prediction()
                typer.secho(
                    f"PVS1 prediction for {self.seqvar.user_repr}: {self.seqvar_prediction.name}",
                    fg=typer.colors.GREEN,
                )
            except Exception as e:
                typer.secho(
                    f"Failed to predict PVS1 for variant {self.seqvar.user_repr}.",
                    err=True,
                    fg=typer.colors.RED,
                )
                typer.secho(e, err=True)
                return

        elif isinstance(variant, StrucVar):
            self.strucvar: StrucVar = variant
            self.strucvar_prediction: PVS1Prediction = PVS1Prediction.NotPVS1  # type: ignore
            self.strucvar_prediction_path: PVS1PredictionStrucVarPath = PVS1PredictionStrucVarPath.NotSet  # type: ignore

            try:
                typer.secho(
                    f"Predicting PVS1 for structural variant {self.strucvar.user_repr}, genome release: {self.genome_release.name}.",
                    fg=typer.colors.BLUE,
                )
                strucvar_pvs1 = StrucVarPVS1(self.strucvar)
                strucvar_pvs1.initialize()
                strucvar_pvs1.verify_PVS1()
                self.strucvar_prediction, self.strucvar_prediction_path = (
                    strucvar_pvs1.get_prediction()
                )
                typer.secho(
                    f"PVS1 prediction for {self.strucvar.user_repr}: {self.strucvar_prediction.name}",
                    fg=typer.colors.GREEN,
                )
            except Exception as e:
                typer.secho(
                    f"Failed to predict PVS1 for structural variant {self.strucvar.user_repr}.",
                    err=True,
                    fg=typer.colors.RED,
                )
                typer.secho(e, err=True)
                return
        else:
            typer.secho(
                f"Failed to resolve variant {self.variant_name}.", err=True, fg=typer.colors.RED
            )
            return
