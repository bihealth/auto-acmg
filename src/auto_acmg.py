"""Implementations of the PVS1 algorithm."""

import typer

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
from src.pvs1.auto_pvs1 import AutoPVS1


class AutoACMG:
    """Class for predicting ACMG criteria.

    This class handles both sequence variants and structural variants to determine their potential
    impact under the various criteria of the ACMG guidelines for variant classification. Currently
    it only implements the PVS1 criterion (not finished yet).

    Attributes:
        variant_name (str): The name or identifier of the variant being analyzed.
        genome_release (GenomeRelease): The genome release version, defaults to GRCh38.
    """

    def __init__(self, variant_name: str, genome_release: GenomeRelease = GenomeRelease.GRCh38):
        """Initializes the AutoACMG with the specified variant and genome release.

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
            Exception: Specific exceptions are caught and logged, but generic exceptions may be
            raised if both resolutions fail.
        """
        try:
            try:
                seqvar_resolver = SeqVarResolver()
                seqvar: SeqVar = seqvar_resolver.resolve_seqvar(
                    self.variant_name, self.genome_release
                )
                return seqvar
            except ParseError:
                strucvar_resolver = StrucVarResolver()
                strucvar: StrucVar = strucvar_resolver.resolve_strucvar(
                    self.variant_name, self.genome_release
                )
                return strucvar
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
        typer.secho(f"Resolving variant: {self.variant_name}.", fg=typer.colors.BLUE)
        variant = self.resolve_variant()
        if not variant:
            typer.secho(
                f"Failed to resolve variant {self.variant_name}.", err=True, fg=typer.colors.RED
            )
            return
        else:
            typer.secho(f"Variant resolved: {variant}.", fg=typer.colors.BLUE)

        if isinstance(variant, SeqVar):
            typer.secho(f"Classifying ACMG criteria for sequence variant.", fg=typer.colors.BLUE)
            # PVS1
            try:
                typer.secho(
                    (
                        f"Predicting PVS1 for variant {variant.user_repr}, genome release: "
                        f"{self.genome_release.name}."
                    ),
                    fg=typer.colors.BLUE,
                )
                self.seqvar: SeqVar = variant
                self.seqvar_prediction: PVS1Prediction = PVS1Prediction.NotPVS1
                self.seqvar_prediction_path: PVS1PredictionSeqVarPath = (
                    PVS1PredictionSeqVarPath.NotSet
                )

                pvs1 = AutoPVS1(self.seqvar, self.genome_release)
                seqvar_prediction, seqvar_prediction_path = pvs1.predict()
                if seqvar_prediction is None or seqvar_prediction_path is None:
                    typer.secho(
                        f"Failed to predict PVS1 for variant {self.seqvar.user_repr}.",
                        err=True,
                        fg=typer.colors.RED,
                    )
                    return
                else:
                    # Double check if the prediction path is indeed for sequence variant
                    assert isinstance(seqvar_prediction_path, PVS1PredictionSeqVarPath)
                    self.seqvar_prediction = seqvar_prediction
                    self.seqvar_prediction_path = seqvar_prediction_path
                    typer.secho(
                        (
                            f"PVS1 prediction for {self.seqvar.user_repr}: "
                            f"{self.seqvar_prediction.name}.\n"
                            f"The prediction path is:\n"
                            f"{PVS1PredictionPathMapping[self.seqvar_prediction_path]}."
                        ),
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
            typer.secho(f"Classifying ACMG criteria for structural variant.", fg=typer.colors.BLUE)
            # PVS1
            try:
                typer.secho(
                    (
                        f"Predicting PVS1 for structural variant {variant.user_repr}, genome "
                        f"release: {self.genome_release.name}."
                    ),
                    fg=typer.colors.BLUE,
                )
                self.strucvar: StrucVar = variant
                self.strucvar_prediction: PVS1Prediction = PVS1Prediction.NotPVS1  # type: ignore
                self.strucvar_prediction_path: PVS1PredictionStrucVarPath = PVS1PredictionStrucVarPath.NotSet  # type: ignore

                pvs1 = AutoPVS1(self.strucvar, self.genome_release)
                strucvar_prediction, strucvar_prediction_path = pvs1.predict()
                if strucvar_prediction is None or strucvar_prediction_path is None:
                    typer.secho(
                        f"Failed to predict PVS1 for structural variant {self.strucvar.user_repr}.",
                        err=True,
                        fg=typer.colors.RED,
                    )
                    return
                else:
                    # Double check if the prediction path is indeed for structural variant
                    assert isinstance(strucvar_prediction_path, PVS1PredictionStrucVarPath)
                    self.strucvar_prediction = strucvar_prediction
                    self.strucvar_prediction_path = strucvar_prediction_path
                    typer.secho(
                        f"PVS1 prediction for {self.strucvar.user_repr}: {self.strucvar_prediction.name}",
                        fg=typer.colors.GREEN,
                    )
            except Exception as e:
                typer.secho(
                    (
                        f"Failed to predict PVS1 for structural variant {self.strucvar.user_repr}."
                        f"The prediction path is:\n"
                        f"{PVS1PredictionPathMapping[self.strucvar_prediction_path]}."
                    ),
                    err=True,
                    fg=typer.colors.RED,
                )
                typer.secho(e, err=True)
                return
        typer.secho(f"ACMG criteria classification completed.", fg=typer.colors.BLUE)
