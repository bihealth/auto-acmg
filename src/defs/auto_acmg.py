from enum import auto
from typing import Dict, List, Optional

from pydantic import BaseModel, ConfigDict

from src.defs.annonars_variant import GnomadExomes, GnomadMtDna
from src.defs.core import AutoAcmgBaseEnum, AutoAcmgBaseModel
from src.defs.genome_builds import GenomeRelease
from src.defs.mehari import Exon, TranscriptGene, TranscriptSeqvar
from src.defs.seqvar import SeqVar
from src.defs.strucvar import StrucVar

# ============ General ACMG Definitions ============


class SpliceType(AutoAcmgBaseEnum):
    """Splicing type enumeration."""

    Donor = "donor"
    Acceptor = "acceptor"
    # Region = "region"
    Unknown = "unknown"


class GenomicStrand(AutoAcmgBaseEnum):
    """Genomic strand enumeration."""

    Plus = "plus"
    Minus = "minus"
    NotSet = "not_set"

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


class AlleleCondition(AutoAcmgBaseEnum):
    """Allele condition enumeration. Used for BS2 criteria."""

    Dominant = "dominant"
    Recessive = "recessive"
    Unknown = "unknown"


#: Mapping of Clingen Dosage Map to Allele Condition
ClingenDosageMap: dict[str, AlleleCondition] = {
    "CLINGEN_DOSAGE_SCORE_UNKNOWN": AlleleCondition.Unknown,
    "CLINGEN_DOSAGE_SCORE_SUFFICIENT_EVIDENCE_AVAILABLE": AlleleCondition.Dominant,
    "CLINGEN_DOSAGE_SCORE_SOME_EVIDENCE_AVAILABLE": AlleleCondition.Dominant,
    "CLINGEN_DOSAGE_SCORE_LITTLE_EVIDENCE": AlleleCondition.Unknown,
    "CLINGEN_DOSAGE_SCORE_NO_EVIDENCE_AVAILABLE": AlleleCondition.Unknown,
    "CLINGEN_DOSAGE_SCORE_RECESSIVE": AlleleCondition.Recessive,
    "CLINGEN_DOSAGE_SCORE_UNLIKELY": AlleleCondition.Unknown,
}


class AminoAcid(AutoAcmgBaseEnum):
    """Amino acids enumeration."""

    Ala = "A"
    Cys = "C"
    Asp = "D"
    Glu = "E"
    Phe = "F"
    Gly = "G"
    His = "H"
    Ile = "I"
    Lys = "K"
    Leu = "L"
    Met = "M"
    Asn = "N"
    Pro = "P"
    Gln = "Q"
    Arg = "R"
    Ser = "S"
    Thr = "T"
    Val = "V"
    Trp = "W"
    Tyr = "Y"
    Stop = "*"


#: Exception list for BA1 criteria.
BA1_ESCEPTION_LIST = [
    SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="chr3",
        pos=128598490,
        delete="C",
        insert="CTAAG",
        user_repr="NM_014049.4:c.-44_-41dupTAAG",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh38,
        chrom="chr3",
        pos=128879647,
        delete="C",
        insert="CTAAG",
        user_repr="NM_014049.4:c.-44_-41dupTAAG",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="chr13",
        pos=20763612,
        delete="C",
        insert="T",
        user_repr="NM_004004.5:c.109G>A",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh38,
        chrom="chr13",
        pos=20189473,
        delete="C",
        insert="T",
        user_repr="NM_004004.5:c.109G>A",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="chr6",
        pos=26091179,
        delete="C",
        insert="G",
        user_repr="NM_000410.3:c.187C>G",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh38,
        chrom="chr6",
        pos=26090951,
        delete="C",
        insert="G",
        user_repr="NM_000410.3:c.187C>G",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="chr6",
        pos=26093141,
        delete="G",
        insert="A",
        user_repr="NM_000410.3:c.845G>A",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh38,
        chrom="chr6",
        pos=26092913,
        delete="G",
        insert="A",
        user_repr="NM_000410.3:c.845G>A",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="chr16",
        pos=3299586,
        delete="G",
        insert="A",
        user_repr="NM_000243.2:c.1105C>T",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh38,
        chrom="chr16",
        pos=3249586,
        delete="G",
        insert="A",
        user_repr="NM_000243.2:c.1105C>T",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="chr16",
        pos=3299468,
        delete="C",
        insert="T",
        user_repr="NM_000243.2:c.1223G>A",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh38,
        chrom="chr16",
        pos=3249468,
        delete="C",
        insert="T",
        user_repr="NM_000243.2:c.1223G>A",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="chr13",
        pos=73409497,
        delete="G",
        insert="A",
        user_repr="NM_006346.2:c.1214G>A",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh38,
        chrom="chr13",
        pos=72835359,
        delete="G",
        insert="A",
        user_repr="NM_006346.2:c.1214G>A",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="chr12",
        pos=121175678,
        delete="C",
        insert="T",
        user_repr="NM_000017.3:c.511C>T",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh38,
        chrom="chr12",
        pos=120737875,
        delete="C",
        insert="T",
        user_repr="NM_000017.3:c.511C>T",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="chr3",
        pos=15686693,
        delete="G",
        insert="C",
        user_repr="NM_000060.4:c.1330G>C",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh38,
        chrom="chr3",
        pos=15645186,
        delete="G",
        insert="C",
        user_repr="NM_000060.4:c.1330G>C",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="chr19",
        pos=39071022,
        delete="G",
        insert="A",
        user_repr="NM_000540.3:c.14524G>A",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh38,
        chrom="chr19",
        pos=38580382,
        delete="G",
        insert="A",
        user_repr="NM_000540.3:c.14524G>A",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="chr19",
        pos=38934252,
        delete="C",
        insert="T",
        user_repr="NM_000540.3:c.325C>T",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh38,
        chrom="chr19",
        pos=38443612,
        delete="C",
        insert="T",
        user_repr="NM_000540.3:c.325C>T",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="chr2",
        pos=152417726,
        delete="C",
        insert="A",
        user_repr="NM_001271208.2:c.19097G>T",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh38,
        chrom="chr2",
        pos=151561212,
        delete="C",
        insert="A",
        user_repr="NM_001271208.2:c.19097G>T",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh37,
        chrom="chr2",
        pos=152382489,
        delete="T",
        insert="G",
        user_repr="NM_001271208.2:c.22249A>C",
    ),
    SeqVar(
        genome_release=GenomeRelease.GRCh38,
        chrom="chr2",
        pos=151525975,
        delete="T",
        insert="G",
        user_repr="NM_001271208.2:c.22249A>C",
    ),
]

# ============ Definitions for algorithms ============


class AutoACMGPrediction(AutoAcmgBaseEnum):
    """ACMG criteria prediction enumeration."""

    NotSet = "not_set"
    NotApplicable = "not_applicable"
    NotAutomated = "not_automated"
    Depricated = "depricated"
    Met = "met"
    NotMet = "not_met"
    Failed = "failed"


class AutoACMGStrength(AutoAcmgBaseEnum):
    """ACMG criteria strength enumeration."""

    NotSet = "not_set"
    PathogenicVeryStrong = "pathogenic_very_strong"
    PathogenicStrong = "pathogenic_strong"
    PathogenicModerate = "pathogenic_moderate"
    PathogenicSupporting = "pathogenic_supporting"
    BenignStandAlone = "benign_stand_alone"
    BenignStrong = "benign_strong"
    BenignSupporting = "benign_supporting"


# ============ Data Structures ============


class CdsInfo(AutoAcmgBaseModel):
    """Information about the coding sequence."""

    start_codon: int
    stop_codon: int
    cds_start: int
    cds_end: int
    cds_strand: GenomicStrand
    exons: List[Exon]


class TranscriptInfo(AutoAcmgBaseModel):
    """Information about a transcript. Used in `SeqVarTranscriptsHelper`."""

    seqvar: Optional[TranscriptSeqvar]
    gene: Optional[TranscriptGene]


class PS1PM5(AutoAcmgBaseModel):
    """PS1 and PM5 criteria prediction."""

    PS1: bool = False
    PS1_strength: AutoACMGStrength = AutoACMGStrength.PathogenicStrong
    PM5: bool = False
    PM5_strength: AutoACMGStrength = AutoACMGStrength.PathogenicModerate


class PM4BP3(AutoAcmgBaseModel):
    """PM4 and BP3 criteria prediction."""

    PM4: bool = False
    PM4_strength: AutoACMGStrength = AutoACMGStrength.PathogenicModerate
    BP3: bool = False
    BP3_strength: AutoACMGStrength = AutoACMGStrength.BenignSupporting


class PM2BA1BS1BS2(AutoAcmgBaseModel):
    """BA1, BS1, BS2, and PM2 criteria prediction."""

    PM2: bool = False
    PM2_strength: AutoACMGStrength = AutoACMGStrength.PathogenicModerate
    BA1: bool = False
    BA1_strength: AutoACMGStrength = AutoACMGStrength.BenignStandAlone
    BS1: bool = False
    BS1_strength: AutoACMGStrength = AutoACMGStrength.BenignStrong
    BS2: bool = False
    BS2_strength: AutoACMGStrength = AutoACMGStrength.BenignStrong


class PM1(AutoAcmgBaseModel):
    """PM1 criteria prediction."""

    PM1: bool = False
    PM1_strength: AutoACMGStrength = AutoACMGStrength.PathogenicModerate


class PP2BP1(AutoAcmgBaseModel):
    """PP2 and BP6 criteria prediction."""

    PP2: bool = False
    PP2_strength: AutoACMGStrength = AutoACMGStrength.PathogenicSupporting
    BP1: bool = False
    BP1_strength: AutoACMGStrength = AutoACMGStrength.BenignSupporting


class BP7(AutoAcmgBaseModel):
    """BP7 criteria prediction."""

    BP7: bool = False
    BP7_strength: AutoACMGStrength = AutoACMGStrength.BenignSupporting


class PP3BP4(AutoAcmgBaseModel):
    """PP3 and BP4 criteria prediction."""

    PP3: bool = False
    PP3_strength: AutoACMGStrength = AutoACMGStrength.PathogenicSupporting
    BP4: bool = False
    BP4_strength: AutoACMGStrength = AutoACMGStrength.BenignSupporting


class AutoACMGCriteria(AutoAcmgBaseModel):
    """Criteria prediction."""

    name: str
    prediction: AutoACMGPrediction = AutoACMGPrediction.NotSet
    strength: AutoACMGStrength = AutoACMGStrength.NotSet
    summary: str = ""
    description: str = ""


class AutoACMGCriteriaPred(AutoAcmgBaseModel):
    """ACMG criteria prediction."""

    pvs1: AutoACMGCriteria = AutoACMGCriteria(
        name="PVS1", strength=AutoACMGStrength.PathogenicVeryStrong
    )
    ps1: AutoACMGCriteria = AutoACMGCriteria(name="PS1", strength=AutoACMGStrength.PathogenicStrong)
    ps2: AutoACMGCriteria = AutoACMGCriteria(
        name="PS2",
        prediction=AutoACMGPrediction.NotAutomated,
        strength=AutoACMGStrength.PathogenicStrong,
    )
    ps3: AutoACMGCriteria = AutoACMGCriteria(
        name="PS3",
        prediction=AutoACMGPrediction.NotAutomated,
        strength=AutoACMGStrength.PathogenicStrong,
    )
    ps4: AutoACMGCriteria = AutoACMGCriteria(
        name="PS4",
        prediction=AutoACMGPrediction.NotAutomated,
        strength=AutoACMGStrength.PathogenicStrong,
    )
    pm1: AutoACMGCriteria = AutoACMGCriteria(
        name="PM1", strength=AutoACMGStrength.PathogenicModerate
    )
    pm2: AutoACMGCriteria = AutoACMGCriteria(
        name="PM2", strength=AutoACMGStrength.PathogenicModerate
    )
    pm3: AutoACMGCriteria = AutoACMGCriteria(
        name="PM3",
        prediction=AutoACMGPrediction.NotAutomated,
        strength=AutoACMGStrength.PathogenicModerate,
    )
    pm4: AutoACMGCriteria = AutoACMGCriteria(
        name="PM4", strength=AutoACMGStrength.PathogenicModerate
    )
    pm5: AutoACMGCriteria = AutoACMGCriteria(
        name="PM5", strength=AutoACMGStrength.PathogenicModerate
    )
    pm6: AutoACMGCriteria = AutoACMGCriteria(
        name="PM6",
        prediction=AutoACMGPrediction.NotAutomated,
        strength=AutoACMGStrength.PathogenicModerate,
    )
    pp1: AutoACMGCriteria = AutoACMGCriteria(
        name="PP1",
        prediction=AutoACMGPrediction.NotAutomated,
        strength=AutoACMGStrength.PathogenicSupporting,
    )
    pp2: AutoACMGCriteria = AutoACMGCriteria(
        name="PP2",
        strength=AutoACMGStrength.PathogenicSupporting,
    )
    pp3: AutoACMGCriteria = AutoACMGCriteria(
        name="PP3",
        strength=AutoACMGStrength.PathogenicSupporting,
    )
    pp4: AutoACMGCriteria = AutoACMGCriteria(
        name="PP4",
        prediction=AutoACMGPrediction.NotAutomated,
        strength=AutoACMGStrength.PathogenicSupporting,
    )
    pp5: AutoACMGCriteria = AutoACMGCriteria(
        name="PP5",
        prediction=AutoACMGPrediction.Depricated,
        strength=AutoACMGStrength.PathogenicSupporting,
    )
    ba1: AutoACMGCriteria = AutoACMGCriteria(name="BA1", strength=AutoACMGStrength.BenignStandAlone)
    bs1: AutoACMGCriteria = AutoACMGCriteria(name="BS1", strength=AutoACMGStrength.BenignStrong)
    bs2: AutoACMGCriteria = AutoACMGCriteria(name="BS2", strength=AutoACMGStrength.BenignStrong)
    bs3: AutoACMGCriteria = AutoACMGCriteria(
        name="BS3",
        prediction=AutoACMGPrediction.NotAutomated,
        strength=AutoACMGStrength.BenignStrong,
    )
    bs4: AutoACMGCriteria = AutoACMGCriteria(
        name="BS4",
        prediction=AutoACMGPrediction.NotAutomated,
        strength=AutoACMGStrength.BenignStrong,
    )
    bp1: AutoACMGCriteria = AutoACMGCriteria(name="BP1", strength=AutoACMGStrength.BenignSupporting)
    bp2: AutoACMGCriteria = AutoACMGCriteria(
        name="BP2",
        prediction=AutoACMGPrediction.NotAutomated,
        strength=AutoACMGStrength.BenignSupporting,
    )
    bp3: AutoACMGCriteria = AutoACMGCriteria(name="BP3", strength=AutoACMGStrength.BenignSupporting)
    bp4: AutoACMGCriteria = AutoACMGCriteria(name="BP4", strength=AutoACMGStrength.BenignSupporting)
    bp5: AutoACMGCriteria = AutoACMGCriteria(
        name="BP5",
        prediction=AutoACMGPrediction.NotAutomated,
        strength=AutoACMGStrength.BenignSupporting,
    )
    bp6: AutoACMGCriteria = AutoACMGCriteria(
        name="BP6",
        prediction=AutoACMGPrediction.Depricated,
        strength=AutoACMGStrength.BenignSupporting,
    )
    bp7: AutoACMGCriteria = AutoACMGCriteria(name="BP7", strength=AutoACMGStrength.BenignSupporting)


class AutoACMGConsequence(AutoAcmgBaseModel):

    mehari: List[str] = []
    cadd: str = ""
    cadd_consequence: str = ""


class AutoACMGCADD(AutoAcmgBaseModel):

    phyloP100: Optional[float] = None
    gerp: Optional[float] = None
    spliceAI_acceptor_gain: Optional[float] = None
    spliceAI_acceptor_loss: Optional[float] = None
    spliceAI_donor_gain: Optional[float] = None
    spliceAI_donor_loss: Optional[float] = None
    ada: Optional[float] = None
    rf: Optional[float] = None


class AutoACMGDbnsfp(AutoAcmgBaseModel):

    alpha_missense: Optional[float] = None
    metaRNN: Optional[float] = None
    bayesDel_noAF: Optional[float] = None
    revel: Optional[float] = None
    phyloP100: Optional[float] = None


class AutoACMGDbscsnv(AutoAcmgBaseModel):

    ada: Optional[float] = None
    rf: Optional[float] = None


class AutoACMGSeqVarScores(AutoAcmgBaseModel):
    """ACMG scores."""

    cadd: AutoACMGCADD = AutoACMGCADD()
    dbnsfp: AutoACMGDbnsfp = AutoACMGDbnsfp()
    dbscsnv: AutoACMGDbscsnv = AutoACMGDbscsnv()


class AutoACMGSeqVarTresholds(AutoAcmgBaseModel):
    """ACMG thresholds."""

    #: Conservation threshold from VarSome
    phyloP100: float = 3.58
    #: Conservarion threshold from Glaucoma VCEP
    gerp: float = 0.0
    #: SpliceAI acceptor gain threshold
    spliceAI_acceptor_gain: float = 0.1
    #: SpliceAI acceptor loss threshold
    spliceAI_acceptor_loss: float = 0.1
    #: SpliceAI donor gain threshold
    spliceAI_donor_gain: float = 0.1
    #: SpliceAI donor loss threshold
    spliceAI_donor_loss: float = 0.1
    #: Ada threshold
    ada: float = 0.957813
    #: RF threshold
    rf: float = 0.584
    #: MetaRNN pathogenic threshold
    metaRNN_pathogenic: float = 0.841
    #: BayesDel_noAF pathogenic threshold
    bayesDel_noAF_pathogenic: float = 0.521
    #: MetaRNN benign threshold
    metaRNN_benign: float = 0.267
    #: BayesDel_noAF benign threshold
    bayesDel_noAF_benign: float = -0.476
    #: PP2 and BP1 pathogenic threshold
    pp2bp1_pathogenic: float = 0.808
    #: PP2 and BP1 benign threshold
    pp2bp1_benign: float = 0.569
    #: PM1 variants count threshold
    pm1_pathogenic: float = 8
    #: BA1 variants count threshold
    ba1_benign: float = 0.05
    #: BS1 variants count threshold
    bs1_benign: float = 0.00015
    #: PM2 variants count threshold
    pm2_pathogenic: float = 0.0001
    #: Minimum number of alleles
    an_min: int = 2000
    #: BP7 donor position
    bp7_donor: int = 1
    #: BP7 acceptor position
    bp7_acceptor: int = 2


class AutoACMGSeqVarData(AutoAcmgBaseModel):

    consequence: AutoACMGConsequence = AutoACMGConsequence()
    gene_symbol: str = ""
    hgnc_id: str = ""
    transcript_id: str = ""
    transcript_tags: List[str] = []
    tx_pos_utr: int = -1
    cds_pos: int = 0
    prot_pos: int = -1
    prot_length: int = -1
    cds_info: Dict[str, CdsInfo] = {}
    pHGVS: str = ""
    cds_start: int = 0
    cds_end: int = 0
    strand: GenomicStrand = GenomicStrand.NotSet
    exons: List[Exon] = []
    scores: AutoACMGSeqVarScores = AutoACMGSeqVarScores()
    thresholds: AutoACMGSeqVarTresholds = AutoACMGSeqVarTresholds()
    gnomad_exomes: Optional[GnomadExomes] = None
    gnomad_mtdna: Optional[GnomadMtDna] = None


class AutoACMGSeqVarResult(AutoAcmgBaseModel):
    """Response of the ACMG criteria prediction for sequence variants."""

    #: Sequence variant for which the ACMG criteria are predicted
    seqvar: Optional[SeqVar] = None
    #: Data, which was used for the prediction
    data: AutoACMGSeqVarData = AutoACMGSeqVarData()
    # ; ACMG criteria prediction
    criteria: AutoACMGCriteriaPred = AutoACMGCriteriaPred()


class AutoACMGStrucVarData(AutoAcmgBaseModel):

    gene_symbol: str = ""
    hgnc_id: str = ""
    transcript_id: str = ""
    transcript_tags: List[str] = []
    prot_pos: int = -1
    prot_length: int = -1
    strand: GenomicStrand = GenomicStrand.NotSet
    exons: List[Exon] = []


class AutoACMGStrucVarPred(AutoAcmgBaseModel):

    pvs1: AutoACMGCriteria = AutoACMGCriteria(
        name="PVS1", strength=AutoACMGStrength.PathogenicVeryStrong
    )


class AutoACMGStrucVarResult(AutoAcmgBaseModel):
    """Response of the ACMG criteria prediction for structural variants."""

    #: Sequence variant for which the ACMG criteria are predicted
    strucvar: Optional[StrucVar] = None
    #: Data, which was used for the prediction
    data: AutoACMGStrucVarData = AutoACMGStrucVarData()
    #: ACMG criteria prediction
    criteria: AutoACMGStrucVarPred = AutoACMGStrucVarPred()


class VcepSpec(BaseModel):
    """VCEP specification for specific gene."""

    #: Identifier, e.g., "GN002"
    identifier: str
    #: Version, e.g., "2.0.0"
    version: str
    #: Title of the VCEP specification, e.g., "ClinGen Cardiomyopathy Expert Panel Specifications
    #: to the ACMG/AMP Variant Interpretation Guidelines for MYH7".
    title: Optional[str] = None

    model_config = ConfigDict(
        frozen=True,
    )
