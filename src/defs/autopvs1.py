from enum import Enum, auto
from typing import List, Optional

from pydantic import BaseModel

from src.defs.mehari import TranscriptGene, TranscriptSeqvar


class TranscriptInfo(BaseModel):
    seqvar: Optional[TranscriptSeqvar]
    gene: Optional[TranscriptGene]


#: Enumeration for sequence variant consequence
class SeqVarConsequence(Enum):
    """Consequence of a sequence variant."""

    NonsenseFrameshift = auto()
    SpliceSites = auto()
    InitiationCodon = auto()
    NotSet = auto()


#: Enumeration for PVS1 prediction status
class PVS1Prediction(Enum):
    """PVS1 prediction."""

    PVS1 = auto()
    PVS1_Strong = auto()
    PVS1_Moderate = auto()
    PVS1_Supporting = auto()
    NotPVS1 = auto()


#: Enumeration for PVS1 prediction path for sequence variant
class PVS1PredictionSeqVarPath(Enum):
    """PVS1 prediction path."""

    NotSet = auto()
    NF1 = auto()
    NF2 = auto()
    NF3 = auto()
    NF4 = auto()
    NF5 = auto()
    NF6 = auto()
    SS1 = auto()
    SS2 = auto()
    SS3 = auto()
    SS4 = auto()
    SS5 = auto()
    SS6 = auto()
    SS7 = auto()
    SS8 = auto()
    SS9 = auto()
    SS10 = auto()
    IC1 = auto()
    IC2 = auto()
    IC3 = auto()


#: Enumeration for PVS1 prediction path for structure variants
class PVS1PredictionsStrucVarPath(Enum):
    """PVS1 prediction path for structure variants."""

    DEL1 = auto()
    DEL2 = auto()
    DEL3 = auto()
    DEL4 = auto()
    DEL5 = auto()
    DEL6 = auto()
    DEL7 = auto()
    DEL8 = auto()
    DUP1 = auto()
    DUP2 = auto()
    DUP3 = auto()
    DUP4 = auto()


#: Mapping of consequence from transcript info to SeqVarConsequence
SeqvarConsequenceMapping = {
    "intergenic_variant": SeqVarConsequence.NotSet,
    "intron_variant": SeqVarConsequence.NotSet,
    "upstream_gene_variant": SeqVarConsequence.NotSet,
    "downstream_gene_variant": SeqVarConsequence.NotSet,
    "5_prime_utr_variant": SeqVarConsequence.NotSet,
    "5_prime_UTR_variant": SeqVarConsequence.NotSet,
    "3_prime_utr_variant": SeqVarConsequence.NonsenseFrameshift,
    "3_prime_UTR_variant": SeqVarConsequence.NonsenseFrameshift,
    "splice_region_variant": SeqVarConsequence.SpliceSites,  # Can affect splicing
    "splice_donor_variant": SeqVarConsequence.SpliceSites,  # Canonical splice site
    "splice_acceptor_variant": SeqVarConsequence.SpliceSites,  # Canonical splice site
    "frameshift_variant": SeqVarConsequence.NonsenseFrameshift,  # Loss of function
    "transcript_ablation": SeqVarConsequence.NotSet,  # Severe effect, but not specifically classified here
    "transcript_amplification": SeqVarConsequence.NotSet,
    "inframe_insertion": SeqVarConsequence.NotSet,  # Not necessarily loss of function
    "inframe_deletion": SeqVarConsequence.NotSet,  # Not necessarily loss of function
    "synonymous_variant": SeqVarConsequence.NotSet,  # Usually benign, but exceptions exist
    "stop_retained_variant": SeqVarConsequence.NotSet,
    "missense_variant": SeqVarConsequence.NotSet,  # Not loss of function in a direct way
    "initiator_codon_variant": SeqVarConsequence.InitiationCodon,  # Affects start codon
    "stop_gained": SeqVarConsequence.NonsenseFrameshift,  # Nonsense variant
    "stop_lost": SeqVarConsequence.NotSet,  # Could be significant, but not classified here as Nonsense/Frameshift
    "mature_mirna_variant": SeqVarConsequence.NotSet,
    "non_coding_exon_variant": SeqVarConsequence.NotSet,  # Impact unclear
    "nc_transcript_variant": SeqVarConsequence.NotSet,
    "incomplete_terminal_codon_variant": SeqVarConsequence.NotSet,  # Rarely significant
    "nmd_transcript_variant": SeqVarConsequence.NotSet,  # Effect on NMD, not directly LOF
    "coding_sequence_variant": SeqVarConsequence.NotSet,  # Ambiguous
    "tfbs_ablation": SeqVarConsequence.NotSet,  # Regulatory, not LOF
    "tfbs_amplification": SeqVarConsequence.NotSet,
    "tf_binding_site_variant": SeqVarConsequence.NotSet,
    "regulatory_region_ablation": SeqVarConsequence.NotSet,  # Regulatory, not LOF
    "regulatory_region_variant": SeqVarConsequence.NotSet,
    "regulatory_region_amplification": SeqVarConsequence.NotSet,
}
