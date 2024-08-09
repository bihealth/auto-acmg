"""PVS1 related definitions."""

from typing import Dict, Union

from src.defs.core import AutoAcmgBaseEnum


class SeqVarPVS1Consequence(AutoAcmgBaseEnum):
    """Consequence of a sequence variant specifically for PVS1."""

    Missense = "missense"
    NonsenseFrameshift = "nonsense_frameshift"
    SpliceSites = "splice_sites"
    InitiationCodon = "initiation_codon"
    NotSet = "not_set"


class PVS1Prediction(AutoAcmgBaseEnum):
    """PVS1 prediction."""

    NotSet = "not_set"
    PVS1 = "pvs1"
    PVS1_Strong = "pvs1_strong"
    PVS1_Moderate = "pvs1_moderate"
    PVS1_Supporting = "pvs1_supporting"
    NotPVS1 = "not_pvs1"
    UnsupportedConsequence = "unsupported_consequence"


#: Enumeration for PVS1 prediction path for sequence variant
class PVS1PredictionSeqVarPath(AutoAcmgBaseEnum):
    """PVS1 prediction path."""

    NotSet = "not_set"
    PTEN = "pten"
    NF1 = "nf1"
    NF2 = "nf2"
    NF3 = "nf3"
    NF4 = "nf4"
    NF5 = "nf5"
    NF6 = "nf6"
    SS1 = "ss1"
    SS2 = "ss2"
    SS3 = "ss3"
    SS4 = "ss4"
    SS5 = "ss5"
    SS6 = "ss6"
    SS7 = "ss7"
    SS8 = "ss8"
    SS9 = "ss9"
    SS10 = "ss10"
    IC1 = "ic1"
    IC2 = "ic2"
    IC3 = "ic3"


#: Enumeration for PVS1 prediction path for structural variant
class PVS1PredictionStrucVarPath(AutoAcmgBaseEnum):
    """PVS1 prediction path for structure variants."""

    NotSet = "not_set"
    DEL1 = "del1"
    DEL2 = "del2"
    DEL3 = "del3"
    DEL4 = "del4"
    DEL5_1 = "del5_1"
    DEL6_1 = "del6_1"
    DEL7_1 = "del7_1"
    DEL5_2 = "del5_2"
    DEL6_2 = "del6_2"
    DEL7_2 = "del7_2"
    DEL8 = "del8"
    DUP1 = "dup1"
    DUP2_1 = "dup2_1"
    DUP2_2 = "dup2_2"
    DUP3 = "dup3"
    DUP4 = "dup4"


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
