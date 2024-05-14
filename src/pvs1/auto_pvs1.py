"""Implementations of the PVS1 algorithm."""

from typing import Optional, Tuple, Union

from loguru import logger

from src.core.config import Config
from src.defs.auto_pvs1 import PVS1Prediction, PVS1PredictionSeqVarPath, PVS1PredictionStrucVarPath
from src.defs.exceptions import AutoAcmgBaseException, AutoPVS1Error
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.defs.strucvar import StrucVar
from src.pvs1.seqvar_pvs1 import SeqVarPVS1
from src.pvs1.strucvar_pvs1 import StrucVarPVS1


class AutoPVS1:
    """Implements the AutoPVS1 algorithm for predicting PVS1 criteria based on genomic variants.

    This class handles both sequence variants and structural variants to determine their potential
    impact under the PVS1 criteria of the ACMG guidelines for variant classification.

    Attributes:
        variant_name (str): The name or identifier of the variant being analyzed.
        genome_release (GenomeRelease): The genome release version, defaults to GRCh38.
    """

    def __init__(
        self,
        variant: SeqVar | StrucVar,
        genome_release: GenomeRelease = GenomeRelease.GRCh38,
        *,
        config: Optional[Config] = None,
    ):
        """Initializes the AutoPVS1 with the specified variant and genome release.

        Args:
            variant_name: The name or identifier of the variant.
            genome_release (Optional): The genome release version, such as GRCh38 or GRCh37.
        """
        #: Configuration to use.
        self.config = config or Config()
        self.variant = variant
        self.genome_release = genome_release

    def predict(
        self,
    ) -> Tuple[PVS1Prediction, Union[PVS1PredictionSeqVarPath, PVS1PredictionStrucVarPath]]:
        """Runs the prediction algorithm to assess the PVS1 criteria for the resolved variant.

        This method resolves the variant and then, based on the type of variant, predicts its
        classification according to the PVS1 criteria. It handles both sequence and structural
        variants.

        Returns:
            Tuple[PVS1Prediction, Union[PVS1PredictionSeqVarPath, PVS1PredictionStrucVarPath]]: The
            prediction result and the path taken to reach it. If the prediction fails, returns None.
        """
        if isinstance(self.variant, SeqVar):
            self.seqvar: SeqVar = self.variant
            self.seqvar_prediction: PVS1Prediction = PVS1Prediction.NotSet
            self.seqvar_prediction_path: PVS1PredictionSeqVarPath = PVS1PredictionSeqVarPath.NotSet

            try:
                seqvar_pvs1 = SeqVarPVS1(self.seqvar, config=self.config)
                seqvar_pvs1.initialize()
                seqvar_pvs1.verify_PVS1()
                self.seqvar_prediction, self.seqvar_prediction_path = seqvar_pvs1.get_prediction()
                return self.seqvar_prediction, self.seqvar_prediction_path
            except AutoAcmgBaseException as e:
                logger.exception("Error occurred: {}", e)
                raise AutoPVS1Error("Error occurred while predicting PVS1.") from e

        elif isinstance(self.variant, StrucVar):
            self.strucvar: StrucVar = self.variant
            self.strucvar_prediction: PVS1Prediction = PVS1Prediction.NotSet  # type: ignore
            self.strucvar_prediction_path: PVS1PredictionStrucVarPath = PVS1PredictionStrucVarPath.NotSet  # type: ignore

            try:
                strucvar_pvs1 = StrucVarPVS1(self.strucvar)
                strucvar_pvs1.initialize()
                strucvar_pvs1.verify_PVS1()
                self.strucvar_prediction, self.strucvar_prediction_path = strucvar_pvs1.get_prediction()
                return self.strucvar_prediction, self.strucvar_prediction_path
            except AutoAcmgBaseException as e:
                logger.exception("Error occurred: {}", e)
                raise AutoPVS1Error("Error occurred while predicting PVS1.") from e

        else:
            raise AutoPVS1Error("Invalid variant type provided for PVS1 prediction.")
