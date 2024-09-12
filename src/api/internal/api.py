from fastapi import APIRouter, HTTPException, Query

from src.auto_acmg import AutoACMG
from src.core.config import settings
from src.defs.api import (
    ApiAutoACMGSeqVarData,
    ApiAutoACMGSeqVarResult,
    SeqVarPredictionResponse,
    StrucVarPredictionResponse,
)
from src.defs.auto_acmg import AutoACMGSeqVarResult, AutoACMGStrucVarResult
from src.defs.exceptions import AutoAcmgBaseException
from src.defs.genome_builds import GenomeRelease

router = APIRouter()


@router.post("/predict/seqvar", response_model=SeqVarPredictionResponse)
async def predict_seqvar(
    variant_name: str = Query(..., description="The name or identifier of the sequence variant"),
    genome_release: str = Query(default="GRCh38", description="The genome release version"),
):
    try:
        genome_release_enum = GenomeRelease.from_string(genome_release)
        if not genome_release_enum:
            raise HTTPException(status_code=400, detail="Invalid genome release")

        auto_acmg = AutoACMG(variant_name, genome_release_enum)
        prediction = auto_acmg.predict()

        if (
            prediction is None
            or not isinstance(prediction, AutoACMGSeqVarResult)
            or prediction.seqvar is None
        ):
            raise HTTPException(
                status_code=400, detail="No valid sequence variant prediction was made"
            )

        # Convert AutoACMGSeqVarResult to ApiAutoACMGSeqVarResult
        api_prediction = ApiAutoACMGSeqVarResult(
            seqvar=prediction.seqvar,
            data=ApiAutoACMGSeqVarData(**prediction.data.model_dump()),
            criteria=prediction.criteria,
        )

        return SeqVarPredictionResponse(prediction=api_prediction)
    except AutoAcmgBaseException as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/predict/strucvar", response_model=StrucVarPredictionResponse)
async def predict_strucvar(
    variant_name: str = Query(..., description="The name or identifier of the structural variant"),
    genome_release: str = Query(default="GRCh38", description="The genome release version"),
    duplication_tandem: bool = Query(
        default=False,
        description="The duplication is in tandem and disrupts reading frame and undergoes NMD",
    ),
):
    try:
        # Set default duplication tandem if provided
        settings.DUPLICATION_TANDEM = duplication_tandem

        genome_release_enum = GenomeRelease.from_string(genome_release)
        if not genome_release_enum:
            raise HTTPException(status_code=400, detail="Invalid genome release")

        auto_acmg = AutoACMG(variant_name, genome_release_enum)
        prediction = auto_acmg.predict()

        if prediction is None or not isinstance(prediction, AutoACMGStrucVarResult):
            raise HTTPException(
                status_code=400, detail="No valid structural variant prediction was made"
            )

        return StrucVarPredictionResponse(prediction=prediction)
    except AutoAcmgBaseException as e:
        raise HTTPException(status_code=400, detail=str(e))
