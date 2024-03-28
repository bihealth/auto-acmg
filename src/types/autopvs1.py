from typing import List, Optional

from pydantic import BaseModel

from src.types.mehari import TranscriptGene, TranscriptSeqvar


class TranscriptInfo(BaseModel):
    seqvar: Optional[TranscriptSeqvar]
    gene: Optional[TranscriptGene]
