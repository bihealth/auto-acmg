from enum import auto
from typing import Dict, List, Optional, Union

from pydantic import BaseModel

from src.defs.core import AutoAcmgBaseEnum
from src.defs.mehari import Exon, TranscriptGene, TranscriptSeqvar


class GenomicStrand(AutoAcmgBaseEnum):
    """Enumeration for genomic strand."""

    Plus = auto()
    Minus = auto()

    @staticmethod
    def from_string(value: str):
        """Converts string to enum member if possible, otherwise returns None."""
        strand_mapping = {
            "STRAND_PLUS": "Plus",
            "STRAND_MINUS": "Minus",
        }
        value_mapped = strand_mapping.get(value, value)
        for member in GenomicStrand:
            if member.name == value_mapped:
                return member
        return None


class TranscriptInfo(BaseModel):
    """Information about a transcript."""

    seqvar: Optional[TranscriptSeqvar]
    gene: Optional[TranscriptGene]


class CdsInfo(BaseModel):
    """Information about the coding sequence."""

    start_codon: int
    stop_codon: int
    cds_start: int
    cds_end: int
    cds_strand: GenomicStrand
    exons: List[Exon]


#: Enumeration for sequence variant consequence
class SeqVarPVS1Consequence(AutoAcmgBaseEnum):
    """Consequence of a sequence variant."""

    Missense = auto()
    NonsenseFrameshift = auto()
    SpliceSites = auto()
    InitiationCodon = auto()
    NotSet = auto()


#: Enumeration for PVS1 prediction status
class PVS1Prediction(AutoAcmgBaseEnum):
    """PVS1 prediction."""

    NotSet = auto()
    PVS1 = auto()
    PVS1_Strong = auto()
    PVS1_Moderate = auto()
    PVS1_Supporting = auto()
    NotPVS1 = auto()
    UnsupportedConsequence = auto()


#: Enumeration for PVS1 prediction path for sequence variant
class PVS1PredictionSeqVarPath(AutoAcmgBaseEnum):
    """PVS1 prediction path."""

    NotSet = auto()
    PTEN = auto()
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


#: Enumeration for PVS1 prediction path for structural variant
class PVS1PredictionStrucVarPath(AutoAcmgBaseEnum):
    """PVS1 prediction path for structure variants."""

    NotSet = auto()
    DEL1 = auto()
    DEL2 = auto()
    DEL3 = auto()
    DEL4 = auto()
    DEL5_1 = auto()
    DEL6_1 = auto()
    DEL7_1 = auto()
    DEL5_2 = auto()
    DEL6_2 = auto()
    DEL7_2 = auto()
    DEL8 = auto()
    DUP1 = auto()
    DUP2_1 = auto()
    DUP2_2 = auto()
    DUP3 = auto()
    DUP4 = auto()


#: Mapping of consequence from transcript info to SeqVarConsequence
SeqvarConsequenceMapping: Dict[str, SeqVarPVS1Consequence] = {
    "intergenic_variant": SeqVarPVS1Consequence.NotSet,
    "intron_variant": SeqVarPVS1Consequence.NotSet,
    "upstream_gene_variant": SeqVarPVS1Consequence.InitiationCodon,
    "downstream_gene_variant": SeqVarPVS1Consequence.InitiationCodon,
    "start_lost": SeqVarPVS1Consequence.InitiationCodon,
    "5_prime_utr_variant": SeqVarPVS1Consequence.NotSet,
    "5_prime_UTR_variant": SeqVarPVS1Consequence.NotSet,
    "3_prime_utr_variant": SeqVarPVS1Consequence.NonsenseFrameshift,
    "3_prime_UTR_variant": SeqVarPVS1Consequence.NonsenseFrameshift,
    "splice_region_variant": SeqVarPVS1Consequence.SpliceSites,  # Can affect splicing
    "splice_donor_variant": SeqVarPVS1Consequence.SpliceSites,  # Canonical splice site
    "splice_donor_5th_base_variant": SeqVarPVS1Consequence.SpliceSites,  # Non-canonical splice site
    "splice_donor_region_variant": SeqVarPVS1Consequence.SpliceSites,  # Non-canonical splice site
    "splice_polypyrimidine_tract_variant": SeqVarPVS1Consequence.SpliceSites,  # Non-canonical splice site
    "splice_acceptor_variant": SeqVarPVS1Consequence.SpliceSites,  # Canonical splice site
    "frameshift_variant": SeqVarPVS1Consequence.NonsenseFrameshift,  # Loss of function
    "transcript_ablation": SeqVarPVS1Consequence.NotSet,  # Severe effect, but not specifically classified here
    "transcript_amplification": SeqVarPVS1Consequence.NotSet,
    "inframe_insertion": SeqVarPVS1Consequence.NotSet,  # Not necessarily loss of function
    "inframe_deletion": SeqVarPVS1Consequence.NotSet,  # Not necessarily loss of function
    "synonymous_variant": SeqVarPVS1Consequence.NotSet,  # Usually benign, but exceptions exist
    "stop_retained_variant": SeqVarPVS1Consequence.NotSet,
    "missense_variant": SeqVarPVS1Consequence.Missense,  # Not loss of function in a direct way
    "initiator_codon_variant": SeqVarPVS1Consequence.InitiationCodon,  # Affects start codon
    "start_retained_variant": SeqVarPVS1Consequence.InitiationCodon,  # Affects start codon
    "stop_gained": SeqVarPVS1Consequence.NonsenseFrameshift,  # Nonsense variant
    "stop_lost": SeqVarPVS1Consequence.NotSet,  # Could be significant, but not classified here as Nonsense/Frameshift
    "mature_mirna_variant": SeqVarPVS1Consequence.NotSet,
    "mature_miRNA_variant": SeqVarPVS1Consequence.NotSet,
    "non_coding_exon_variant": SeqVarPVS1Consequence.NotSet,  # Impact unclear
    "nc_transcript_variant": SeqVarPVS1Consequence.NotSet,
    "incomplete_terminal_codon_variant": SeqVarPVS1Consequence.NotSet,  # Rarely significant
    "NMD_transcript_variant": SeqVarPVS1Consequence.NotSet,  # Effect on NMD, not directly LOF
    "nmd_transcript_variant": SeqVarPVS1Consequence.NotSet,  # Effect on NMD, not directly LOF
    "coding_sequence_variant": SeqVarPVS1Consequence.NotSet,  # Ambiguous
    "sequence_variant": SeqVarPVS1Consequence.NotSet,  # Ambiguous
    "tfbs_ablation": SeqVarPVS1Consequence.NotSet,  # Regulatory, not LOF
    "tfbs_amplification": SeqVarPVS1Consequence.NotSet,
    "tf_binding_site_variant": SeqVarPVS1Consequence.NotSet,
    "regulatory_region_ablation": SeqVarPVS1Consequence.NotSet,  # Regulatory, not LOF
    "regulatory_region_variant": SeqVarPVS1Consequence.NotSet,
    "regulatory_region_amplification": SeqVarPVS1Consequence.NotSet,
    "feature_elongation": SeqVarPVS1Consequence.NotSet,  # Ambiguous
    "feature_truncation": SeqVarPVS1Consequence.NotSet,  # Ambiguous
    "protein_altering_variant": SeqVarPVS1Consequence.NotSet,  # Ambiguous
    "non_coding_transcript_exon_variant": SeqVarPVS1Consequence.NotSet,  # Impact unclear
    "non_coding_transcript_variant": SeqVarPVS1Consequence.NotSet,  # Impact unclear
    "coding_transcript_variant": SeqVarPVS1Consequence.NotSet,  # Ambiguous
}


#: Mapping of PVS1 prediction path to description for sequence variant
PVS1PredictionPathMapping: Dict[
    Union[PVS1PredictionSeqVarPath, PVS1PredictionStrucVarPath], str
] = {
    PVS1PredictionSeqVarPath.NotSet: "Not Set",
    PVS1PredictionSeqVarPath.PTEN: "Special guideline for PTEN -> Predicted to undergo NMD",
    PVS1PredictionSeqVarPath.NF1: (
        "Predicted to undergo NMD -> Exon is present in biologically-relevant transcript(s)"
    ),
    PVS1PredictionSeqVarPath.NF2: (
        "Predicted to undergo NMD -> Exon is absent from biologically-relevant transcript(s)"
    ),
    PVS1PredictionSeqVarPath.NF3: (
        "Not predicted to undergo NMD -> "
        "Truncated/altered region is critical to protein function"
    ),
    PVS1PredictionSeqVarPath.NF4: (
        "Not predicted to undergo NMD -> "
        "Role of region in protein function is unknown -> "
        "LoF variants in this exon are frequent in the general population and/or "
        "exon is absent from biologically-relevant transcript(s)"
    ),
    PVS1PredictionSeqVarPath.NF5: (
        "Not predicted to undergo NMD -> "
        "Role of region in protein function is unknown -> "
        "LoF variants in this exon are not frequent in the general population and "
        "exon is present in biologically-relevant transcript(s) -> "
        "Variant removes >10% of protein"
    ),
    PVS1PredictionSeqVarPath.NF6: (
        "Not predicted to undergo NMD -> "
        "Role of region in protein function is unknown -> "
        "LoF variants in this exon are not frequent in the general population and "
        "exon is present in biologically-relevant transcript(s) -> "
        "Variant removes <10% of protein"
    ),
    PVS1PredictionSeqVarPath.SS1: (
        "Exon skipping or use of a cryptic slice site disrupts reading frame and "
        "is predicted to undergo NMD -> "
        "Exon is present in biologically-relevant transcript(s)"
    ),
    PVS1PredictionSeqVarPath.SS2: (
        "Exon skipping or use of a cryptic slice site disrupts reading frame and "
        "is predicted to undergo NMD -> "
        "Exon is absent from biologically-relevant transcript(s)"
    ),
    PVS1PredictionSeqVarPath.SS3: (
        "Exon skipping or use of a cryptic slice site disrupts reading frame and "
        "is not predicted to undergo NMD -> "
        "Truncated/altered region is critical to protein function"
    ),
    PVS1PredictionSeqVarPath.SS4: (
        "Exon skipping or use of a cryptic slice site disrupts reading frame and "
        "is not predicted to undergo NMD -> "
        "Role of region in protein function is unknown -> "
        "LoF variants in this exon are frequent in the general population and/or "
        "exon is absent from biologically-relevant transcript(s)"
    ),
    PVS1PredictionSeqVarPath.SS5: (
        "Exon skipping or use of a cryptic slice site disrupts reading frame and "
        "is not predicted to undergo NMD -> "
        "Role of region in protein function is unknown -> "
        "LoF variants in this exon are not frequent in the general population and "
        "exon is present in biologically-relevant transcript(s) -> "
        "Variant removes >10% of protein"
    ),
    PVS1PredictionSeqVarPath.SS6: (
        "Exon skipping or use of a cryptic slice site disrupts reading frame and "
        "is not predicted to undergo NMD -> "
        "Role of region in protein function is unknown -> "
        "LoF variants in this exon are not frequent in the general population and "
        "exon is present in biologically-relevant transcript(s) -> "
        "Variant removes <10% of protein"
    ),
    PVS1PredictionSeqVarPath.SS7: (
        "Exon skipping or use of a cryptic slice site preserves reading frame -> "
        "Role of region in protein function is unknown -> "
        "LoF variants in this exon are frequent in the general population and/or "
        "exon is absent from biologically-relevant transcript(s)"
    ),
    PVS1PredictionSeqVarPath.SS8: (
        "Exon skipping or use of a cryptic slice site preserves reading frame -> "
        "Role of region in protein function is unknown -> "
        "LoF variants in this exon are not frequent in the general population and "
        "exon is present in biologically-relevant transcript(s) -> "
        "Variant removes >10% of protein"
    ),
    PVS1PredictionSeqVarPath.SS9: (
        "Exon skipping or use of a cryptic slice site preserves reading frame -> "
        "Role of region in protein function is unknown -> "
        "LoF variants in this exon are not frequent in the general population and "
        "exon is present in biologically-relevant transcript(s) -> "
        "Variant removes <10% of protein"
    ),
    PVS1PredictionSeqVarPath.SS10: (
        "Exon skipping or use of a cryptic slice site preserves reading frame -> "
        "Truncated/altered region is critical to protein function"
    ),
    PVS1PredictionSeqVarPath.IC1: (
        "No known alternative start codon in other transcripts -> "
        ">=1 pathogenic variant(s) upstream of closest potential in-frame start codon"
    ),
    PVS1PredictionSeqVarPath.IC2: (
        "No known alternative start codon in other transcripts -> "
        "No pathogenic variant(s) upstream of closest potential in-frame start codon"
    ),
    PVS1PredictionSeqVarPath.IC3: "Different functional transcript uses alternative start codon",
    PVS1PredictionStrucVarPath.NotSet: "Not Set",
    PVS1PredictionStrucVarPath.DEL1: "Full gene deletion",
    PVS1PredictionStrucVarPath.DEL2: (
        "Single to multi exon deletion disrupts reading frame and "
        "is predicted to undergo NMD -> "
        "Exon is present in biologically-relevant transcript(s)"
    ),
    PVS1PredictionStrucVarPath.DEL3: (
        "Single to multi exon deletion disrupts reading frame and "
        "is predicted to undergo NMD -> "
        "Exon is absent from biologically-relevant transcript(s)"
    ),
    PVS1PredictionStrucVarPath.DEL4: (
        "Single to multi exon deletion disrupts reading frame and "
        "is not predicted to undergo NMD -> "
        "Truncated/altered region is critical to protein function"
    ),
    PVS1PredictionStrucVarPath.DEL5_1: (
        "Single to multi exon deletion disrupts reading frame and "
        "is not predicted to undergo NMD -> "
        "Role of region in protein function is unknown -> "
        "LoF variants in this exon are frequent in the general population and/or "
        "exon is absent from biologically-relevant transcript(s)"
    ),
    PVS1PredictionStrucVarPath.DEL6_1: (
        "Single to multi exon deletion disrupts reading frame and "
        "is not predicted to undergo NMD -> "
        "Role of region in protein function is unknown -> "
        "LoF variants in this exon are not frequent in the general population and "
        "exon is present in biologically-relevant transcript(s) -> "
        "Variant removes >10% of protein"
    ),
    PVS1PredictionStrucVarPath.DEL7_1: (
        "Single to multi exon deletion disrupts reading frame and "
        "is not predicted to undergo NMD -> "
        "Role of region in protein function is unknown -> "
        "LoF variants in this exon are not frequent in the general population and "
        "exon is present in biologically-relevant transcript(s) -> "
        "Variant removes <10% of protein"
    ),
    PVS1PredictionStrucVarPath.DEL5_2: (
        "Single to multi exon deletion preserves reading frame -> "
        "Role of region in protein function is unknown -> "
        "LoF variants in this exon are frequent in the general population and/or "
        "exon is absent from biologically-relevant transcript(s)"
    ),
    PVS1PredictionStrucVarPath.DEL6_2: (
        "Single to multi exon deletion preserves reading frame -> "
        "Role of region in protein function is unknown -> "
        "LoF variants in this exon are not frequent in the general population and "
        "exon is present in biologically-relevant transcript(s) -> "
        "Variant removes >10% of protein"
    ),
    PVS1PredictionStrucVarPath.DEL7_2: (
        "Single to multi exon deletion preserves reading frame -> "
        "Role of region in protein function is unknown -> "
        "LoF variants in this exon are not frequent in the general population and "
        "exon is present in biologically-relevant transcript(s) -> "
        "Variant removes <10% of protein"
    ),
    PVS1PredictionStrucVarPath.DEL8: (
        "Single to multi exon deletion preserves reading frame -> "
        "Truncated/altered region is critical to protein function"
    ),
    PVS1PredictionStrucVarPath.DUP1: (
        "Proven in tandem -> " "Reading frame disrupted and NMD predicted to occur"
    ),
    PVS1PredictionStrucVarPath.DUP2_1: (
        "Proven in tandem -> " "No or unknown impact on reading frame and NMD"
    ),
    PVS1PredictionStrucVarPath.DUP2_2: (
        "Presumed in tandem -> " "No or unknown impact on reading frame and NMD"
    ),
    PVS1PredictionStrucVarPath.DUP3: (
        "Proven in tandem -> " "Reading frame presumed disrupted and NMD predicted to occur"
    ),
    PVS1PredictionStrucVarPath.DUP4: "Proven not in tandem",
}
