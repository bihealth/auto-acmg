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
    """Implements the AutoPVS1 algorithm for predicting PVS1 criteria based on genomic variants."""

    def __init__(
        self,
        variant: Union[SeqVar, StrucVar],
        *,
        config: Optional[Config] = None,
    ):
        """Initializes the AutoPVS1 with the specified variant and genome release.

        Args:
            variant_name: The name or identifier of the variant.
        """
        #: Configuration to use.
        self.config: Config = config or Config()
        self.variant: Union[SeqVar, StrucVar] = variant

    def predict(
        self,
    ) -> Tuple[PVS1Prediction, Union[PVS1PredictionSeqVarPath, PVS1PredictionStrucVarPath], str]:
        """Runs the prediction algorithm to assess the PVS1 criteria for the resolved variant.

        This method resolves the variant and then, based on the type of variant, predicts its
        classification according to the PVS1 criteria. It handles both sequence and structural
        variants.

        Returns:
            The prediction result, the path taken to reach it and comment with details of
            prediction. If the prediction fails, returns None.
        """
        if isinstance(self.variant, SeqVar):
            try:
                seqvar_pvs1 = SeqVarPVS1(self.variant, config=self.config)
                seqvar_pvs1.initialize()
                return seqvar_pvs1.verify_PVS1()
            except AutoAcmgBaseException as e:
                logger.exception("Error occurred: {}", e)
                raise AutoPVS1Error("Error occurred while predicting PVS1.") from e

        elif isinstance(self.variant, StrucVar):
            try:
                strucvar_pvs1 = StrucVarPVS1(self.variant, config=self.config)
                strucvar_pvs1.initialize()
                return strucvar_pvs1.verify_PVS1()
            except AutoAcmgBaseException as e:
                logger.exception("Error occurred: {}", e)
                raise AutoPVS1Error("Error occurred while predicting PVS1.") from e

        else:
            raise AutoPVS1Error("Invalid variant type provided for PVS1 prediction.")
