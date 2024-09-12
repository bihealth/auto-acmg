"""API response models for AutoACMG."""

from typing import List, Union

from pydantic import BaseModel, Field

from src.defs.auto_acmg import (
    AutoACMGConsequence,
    AutoACMGCriteriaPred,
    AutoACMGSeqVarScores,
    AutoACMGSeqVarTresholds,
    AutoACMGStrucVarResult,
    GenomicStrand,
)
from src.defs.mehari import Exon
from src.defs.seqvar import SeqVar
from src.defs.strucvar import StrucVar


class ApiAutoACMGSeqVarData(BaseModel):
    consequence: AutoACMGConsequence = Field(
        ..., description="The consequence of the sequence variant"
    )
    gene_symbol: str = Field(..., description="The gene symbol of the sequence variant")
    hgnc_id: str = Field(..., description="The HGNC ID of the sequence variant")
    transcript_id: str = Field(..., description="The transcript ID of the sequence variant")
    transcript_tags: List[str] = Field(
        ..., description="The transcript tags of the sequence variant"
    )
    tx_pos_utr: int = Field(
        ..., description="The transcript position in the UTR of the sequence variant"
    )
    cds_pos: int = Field(..., description="The CDS position of the sequence variant")
    prot_pos: int = Field(..., description="The protein position of the sequence variant")
    prot_length: int = Field(..., description="The protein length of the sequence variant")
    pHGVS: str = Field(..., description="The pHGVS notation of the sequence variant")
    cds_start: int = Field(..., description="The CDS start position of the sequence variant")
    cds_end: int = Field(..., description="The CDS end position of the sequence variant")
    strand: GenomicStrand = Field(..., description="The strand of the sequence variant")
    exons: List[Exon] = Field(..., description="The exons of the sequence variant")
    scores: AutoACMGSeqVarScores = Field(..., description="The scores of the sequence variant")
    thresholds: AutoACMGSeqVarTresholds = Field(
        ..., description="The thresholds of the sequence variant"
    )

    class Config:
        extra = "ignore"  # This will ignore any extra fields not explicitly defined


class ApiAutoACMGSeqVarResult(BaseModel):
    seqvar: SeqVar
    data: ApiAutoACMGSeqVarData
    criteria: AutoACMGCriteriaPred


class SeqVarPredictionResponse(BaseModel):
    prediction: ApiAutoACMGSeqVarResult = Field(
        ..., description="The prediction result for the sequence variant"
    )


class StrucVarPredictionResponse(BaseModel):
    prediction: AutoACMGStrucVarResult = Field(
        ..., description="The prediction result for the structural variant"
    )


class VariantResolveResponse(BaseModel):
    variant_type: str = Field(..., description="The type of the variant (sequence or structural)")
    resolved_variant: Union[SeqVar, StrucVar] = Field(..., description="The resolved variant")
