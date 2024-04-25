"""Implementations of the PVS1 algorithm."""

import typer

from src.core.config import settings
from src.defs.autopvs1 import (
    PVS1Prediction,
    PVS1PredictionPathMapping,
    PVS1PredictionSeqVarPath,
    PVS1PredictionStrucVarPath,
)
from src.defs.exceptions import ParseError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar, SeqVarResolver
from src.defs.strucvar import StrucVar, StrucVarResolver
from src.pvs1.seqvar_pvs1 import SeqVarPVS1
from src.pvs1.strucvar_pvs1 import StrucVarPVS1


class AutoACMG:
    """Implements the AutoPVS1 algorithm for predicting PVS1 criteria based on genomic variants.

    This class handles both sequence variants and structural variants to determine their potential
    impact under the PVS1 criteria of the ACMG guidelines for variant classification.

    Attributes:
        variant_name (str): The name or identifier of the variant being analyzed.
        genome_release (GenomeRelease): The genome release version, defaults to GRCh38.
    """

    def __init__(self, variant_name: str, genome_release: GenomeRelease = GenomeRelease.GRCh38):
        """Initializes the AutoPVS1 with the specified variant and genome release.

        Args:
            variant_name: The name or identifier of the variant.
            genome_release (Optional): The genome release version, such as GRCh38 or GRCh37.
        """
        self.variant_name = variant_name
        self.genome_release = genome_release

    def resolve_variant(self) -> SeqVar | StrucVar | None:
        """Attempts to resolve the specified variant as either a sequence or structural variant.

        This method first tries to resolve the variant as a sequence variant. If it fails, it then
        attempts to resolve it as a structural variant.

        Returns:
            SeqVar, StrucVar, or None: The resolved variant object or None if resolution fails.

        Raises:
            Exception: Specific exceptions are caught and logged, but generic exceptions may be raised if both resolutions fail.
        """
        try:
            try:
                seqvar_resolver = SeqVarResolver()
                seqvar: SeqVar = seqvar_resolver.resolve_seqvar(
                    self.variant_name, self.genome_release
                )
                typer.secho(f"Resolved variant: {seqvar}.", fg=typer.colors.BLUE)
                return seqvar
            except ParseError:
                strucvar_resolver = StrucVarResolver()
                strucvar: StrucVar = strucvar_resolver.resolve_strucvar(
                    self.variant_name, self.genome_release
                )
                typer.secho(f"Resolved structural variant: {strucvar}.", fg=typer.colors.BLUE)
                return strucvar
        except ParseError:
            return None
        except Exception as e:
            typer.secho(e, err=True, fg=typer.colors.RED)
            return None

    def predict(self):
        """Runs the prediction algorithm to assess the PVS1 criteria for the resolved variant.

        This method resolves the variant and then, based on the type of variant, predicts its
        classification according to the PVS1 criteria. It handles both sequence and structural variants.

        Raises:
            Exception: Handles general exceptions that may occur during prediction and logs them.
        """
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
                    (
                        f"PVS1 prediction for {self.seqvar.user_repr}: "
                        f"{self.seqvar_prediction.name}.\n"
                        f"The prediction path is:\n"
                        f"{PVS1PredictionPathMapping[self.seqvar_prediction_path]}."
                    ),
                    fg=typer.colors.GREEN,
                )
                if settings.DEBUG:
                    typer.secho(
                        (
                            f"\nConsequence: {seqvar_pvs1._consequence.name},\n"
                            f"pHGVS: {seqvar_pvs1.pHGVS},\n"
                            f"Transcript tags: {seqvar_pvs1.transcript_tags},\n"
                            f"Exons: {seqvar_pvs1.exons}."
                        ),
                        fg=typer.colors.YELLOW,
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
