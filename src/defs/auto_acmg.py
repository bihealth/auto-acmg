from enum import auto
from typing import Dict, List, Optional

from pydantic import BaseModel, ConfigDict

from src.defs.annonars_variant import GnomadExomes, GnomadMtDna
from src.defs.auto_pvs1 import CdsInfo, GenomicStrand
from src.defs.core import AutoAcmgBaseEnum
from src.defs.mehari import Exon
from src.defs.seqvar import SeqVar


class SpliceType(AutoAcmgBaseEnum):
    """Splice type enumeration."""

    Donor = "donor"
    Acceptor = "acceptor"
    # Region = "region"
    Unknown = "unknown"


class AlleleCondition(AutoAcmgBaseEnum):
    """Allele condition enumeration."""

    Dominant = "dominant"
    Recessive = "recessive"
    Unknown = "unknown"


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
    """Amino acid enumeration."""

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


class MissenseScore(BaseModel):
    """Missense score."""

    name: str
    benign_threshold: float = 0.0
    pathogenic_threshold: float = 0.0


# MissenseScores: List[MissenseScore] = [
#     MissenseScore(name="BayesDel_noAF_score", benign_threshold=-0.476, pathogenic_threshold=0.521),
#     MissenseScore(name="REVEL_score", benign_threshold=0.133, pathogenic_threshold=0.946),
#     MissenseScore(name="CADD_raw", benign_threshold=16.1, pathogenic_threshold=33),
#     MissenseScore(name="PrimateAI_score", benign_threshold=0.362, pathogenic_threshold=0.895),
#     # MissenseScore(name="FATHMM", benign_threshold=4.69, pathogenic_threshold=-5.04),  # is <=
#     # MissenseScore(name="MutPred2", benign_threshold=0.197, pathogenic_threshold=0.829),  # We don't have the right score for this
#     MissenseScore(name="Polyphen2_HVAR_score", benign_threshold=0.001, pathogenic_threshold=1),
#     # MissenseScore(name="SIFT_score", benign_threshold=0.327, pathogenic_threshold=0.001),   # is <=
#     # MissenseScore(name="VEST4_score", benign_threshold=0.302, pathogenic_threshold=0.861),
#     MissenseScore(
#         name="phyloP100way_vertebrate", benign_threshold=-1.04, pathogenic_threshold=9.88
#     ),  # Not sure about 100/470/17 way
# ]

MissenseScores: List[MissenseScore] = [
    MissenseScore(name="MetaRNN_score", benign_threshold=0.267, pathogenic_threshold=0.841),
    MissenseScore(name="BayesDel_noAF_score", benign_threshold=-0.476, pathogenic_threshold=0.521),
]


class ACMGPrediction(AutoAcmgBaseEnum):
    """ACMG prediction enumeration."""

    NotSet = auto()
    Met = auto()
    Unmet = auto()


class ACMGCriteria(BaseModel):

    prediction: ACMGPrediction = ACMGPrediction.NotSet
    comment: str = ""


class ACMGResult(BaseModel):
    """ACMG criteria prediction. Note: without PVS1."""

    PS1: ACMGCriteria = ACMGCriteria()
    PM1: ACMGCriteria = ACMGCriteria()
    PM2: ACMGCriteria = ACMGCriteria()
    PM4: ACMGCriteria = ACMGCriteria()
    PM5: ACMGCriteria = ACMGCriteria()
    PP2: ACMGCriteria = ACMGCriteria()
    PP3: ACMGCriteria = ACMGCriteria()
    BA1: ACMGCriteria = ACMGCriteria()
    BS1: ACMGCriteria = ACMGCriteria()
    BS2: ACMGCriteria = ACMGCriteria()
    BP1: ACMGCriteria = ACMGCriteria()
    BP3: ACMGCriteria = ACMGCriteria()
    BP4: ACMGCriteria = ACMGCriteria()
    BP7: ACMGCriteria = ACMGCriteria()


class AutoACMGPrediction(AutoAcmgBaseEnum):
    """ACMG criteria prediction enumeration."""

    NotSet = auto()
    NotApplicable = auto()
    NotAutomated = auto()
    Depricated = auto()
    Met = auto()
    NotMet = auto()
    Failed = auto()


class AutoACMGStrength(AutoAcmgBaseEnum):
    """ACMG criteria strength enumeration."""

    NotSet = auto()
    PathogenicVeryStrong = auto()
    PathogenicStrong = auto()
    PathogenicModerate = auto()
    PathogenicSupporting = auto()
    BenignStandAlone = auto()
    BenignStrong = auto()
    BenignSupporting = auto()


class PS1PM5(BaseModel):
    """PS1 and PM5 criteria prediction."""

    PS1: bool = False
    PS1_strength: AutoACMGStrength = AutoACMGStrength.PathogenicStrong
    PM5: bool = False
    PM5_strength: AutoACMGStrength = AutoACMGStrength.PathogenicModerate


class PM4BP3(BaseModel):
    """PM4 and BP3 criteria prediction."""

    PM4: bool = False
    PM4_strength: AutoACMGStrength = AutoACMGStrength.PathogenicModerate
    BP3: bool = False
    BP3_strength: AutoACMGStrength = AutoACMGStrength.BenignSupporting


class PM2BA1BS1BS2(BaseModel):
    """BA1, BS1, BS2, and PM2 criteria prediction."""

    PM2: bool = False
    PM2_strength: AutoACMGStrength = AutoACMGStrength.PathogenicModerate
    BA1: bool = False
    BA1_strength: AutoACMGStrength = AutoACMGStrength.BenignStandAlone
    BS1: bool = False
    BS1_strength: AutoACMGStrength = AutoACMGStrength.BenignStrong
    BS2: bool = False
    BS2_strength: AutoACMGStrength = AutoACMGStrength.BenignStrong


class PM1(BaseModel):
    """PM1 criteria prediction."""

    PM1: bool = False
    PM1_strength: AutoACMGStrength = AutoACMGStrength.PathogenicModerate


class PP2BP1(BaseModel):
    """PP2 and BP6 criteria prediction."""

    PP2: bool = False
    PP2_strength: AutoACMGStrength = AutoACMGStrength.PathogenicSupporting
    BP1: bool = False
    BP1_strength: AutoACMGStrength = AutoACMGStrength.BenignSupporting


class BP7(BaseModel):
    """BP7 criteria prediction."""

    BP7: bool = False
    BP7_strength: AutoACMGStrength = AutoACMGStrength.BenignSupporting


class PP3BP4(BaseModel):
    """PP3 and BP4 criteria prediction."""

    PP3: bool = False
    PP3_strength: AutoACMGStrength = AutoACMGStrength.PathogenicSupporting
    BP4: bool = False
    BP4_strength: AutoACMGStrength = AutoACMGStrength.BenignSupporting


class AutoACMGCriteria(BaseModel):
    """Criteria prediction."""

    name: str
    prediction: AutoACMGPrediction = AutoACMGPrediction.NotSet
    strength: AutoACMGStrength = AutoACMGStrength.NotSet
    summary: str = ""
    description: str = ""


class AutoACMGCriteriaResult(BaseModel):

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


class AutoACMGConsequence(BaseModel):

    mehari: List[str] = []
    cadd: str = ""
    cadd_consequence: str = ""


class AutoACMGCADD(BaseModel):

    phyloP100: Optional[float] = None
    spliceAI_acceptor_gain: Optional[float] = None
    spliceAI_acceptor_loss: Optional[float] = None
    spliceAI_donor_gain: Optional[float] = None
    spliceAI_donor_loss: Optional[float] = None
    ada: Optional[float] = None
    rf: Optional[float] = None


class AutoACMGDbnsfp(BaseModel):

    alpha_missense: Optional[float] = None
    metaRNN: Optional[float] = None
    bayesDel_noAF: Optional[float] = None
    revel: Optional[float] = None
    phyloP100: Optional[float] = None


class AutoACMGDbscsnv(BaseModel):

    ada: Optional[float] = None
    rf: Optional[float] = None


class AutoACMGScores(BaseModel):
    """ACMG scores."""

    cadd: AutoACMGCADD = AutoACMGCADD()
    dbnsfp: AutoACMGDbnsfp = AutoACMGDbnsfp()
    dbscsnv: AutoACMGDbscsnv = AutoACMGDbscsnv()


class AutoACMGTresholds(BaseModel):
    """ACMG thresholds."""

    #: Conservation threshold from VarSome
    phyloP100: float = 3.58
    #: SpliceAI acceptor gain threshold
    spliceAI_acceptor_gain: float = 0.2
    #: SpliceAI acceptor loss threshold
    spliceAI_acceptor_loss: float = 0.2
    #: SpliceAI donor gain threshold
    spliceAI_donor_gain: float = 0.2
    #: SpliceAI donor loss threshold
    spliceAI_donor_loss: float = 0.2
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


class AutoACMGData(BaseModel):
    """Response of the ACMG criteria prediction."""

    consequence: AutoACMGConsequence = AutoACMGConsequence()
    gene_symbol: str = ""
    hgnc_id: str = ""
    transcript_id: str = ""
    transcript_tags: List[str] = []
    tx_pos_utr: int = -1
    prot_pos: int = -1
    prot_length: int = -1
    cds_info: Dict[str, CdsInfo] = {}
    pHGVS: str = ""
    cds_start: int = 0
    cds_end: int = 0
    strand: Optional[GenomicStrand] = None
    exons: List[Exon] = []
    scores: AutoACMGScores = AutoACMGScores()
    thresholds: AutoACMGTresholds = AutoACMGTresholds()
    gnomad_exomes: Optional[GnomadExomes] = None
    gnomad_mtdna: Optional[GnomadMtDna] = None


class AutoACMGResult(BaseModel):
    """Response of the ACMG criteria prediction."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    #: Sequence variant for which the ACMG criteria are predicted
    seqvar: Optional[SeqVar] = None
    #: Data, which was used for the prediction
    data: AutoACMGData = AutoACMGData()
    # ; ACMG criteria prediction
    criteria: AutoACMGCriteriaResult = AutoACMGCriteriaResult()
