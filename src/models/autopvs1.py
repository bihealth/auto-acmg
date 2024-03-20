from typing import List, Optional

from pydantic import BaseModel

from src.models.mehari_gene import TranscriptGene
from src.models.mehari_seqvar import TranscriptSeqvar


class TranscriptInfo(BaseModel):
    seqvar: Optional[TranscriptSeqvar]
    gene: Optional[TranscriptGene]
