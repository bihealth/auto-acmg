from enum import Enum, EnumMeta, auto
from typing import Any, List, Optional

from pydantic import BaseModel

from src.defs.exceptions import AutoAcmgBaseException


class AutoAcmgBaseEnumMeta(EnumMeta):
    def __getitem__(cls, name: Any) -> Any:
        """Override __getitem__ to raise an AutoAcmgBaseException if the KeyError is raised."""
        try:
            return super().__getitem__(name)
        except KeyError:
            raise AutoAcmgBaseException(f"Invalid {cls.__name__} value: {name}")


class AutoAcmgBaseEnum(Enum, metaclass=AutoAcmgBaseEnumMeta):
    """Base enumeration for ACMG criteria."""

    @classmethod
    def _missing_(cls, value: Any) -> Any:
        """Override _missing_ to raise an AutoAcmgBaseException if the ValueError is raised."""
        raise AutoAcmgBaseException(f"Invalid {cls.__name__} value: {value}")


class SpliceType(AutoAcmgBaseEnum):
    """Splice type enumeration."""

    Donor = "donor"
    Acceptor = "acceptor"
    # Region = "region"
    Unknown = "unknown"


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
    benign_threshold: Optional[float] = None
    pathogenic_threshold: Optional[float] = None


MissenseScores: List[MissenseScore] = [
    MissenseScore(name="BayesDel_noAF_score", benign_threshold=-0.36, pathogenic_threshold=0.27),
    MissenseScore(name="REVEL_score", benign_threshold=0.183, pathogenic_threshold=0.773),
    MissenseScore(name="CADD_raw", benign_threshold=0.15, pathogenic_threshold=28.1),
    MissenseScore(name="PrimateAI_score", benign_threshold=0.362, pathogenic_threshold=0.867),
    # MissenseScore(name="FATHMM", benign_threshold=4.69, pathogenic_threshold=-5.04),  # is <=
    # MissenseScore(name="MutPred2", benign_threshold=0.197, pathogenic_threshold=0.829),  # We don't have the right score for this
    MissenseScore(name="Polyphen2_HVAR_score", benign_threshold=0.009, pathogenic_threshold=0.999),
    # MissenseScore(name="SIFT_score", benign_threshold=0.327, pathogenic_threshold=0.001),   # is <=
    MissenseScore(name="VEST4_score", benign_threshold=0.302, pathogenic_threshold=0.861),
    MissenseScore(
        name="phyloP100way_vertebrate", benign_threshold=0.021, pathogenic_threshold=9.741
    ),  # Not sure about 100/470/17 way
]


class ACMGCriteria(BaseModel):
    """ACMG criteria prediction. Note: without PVS1."""

    PS1: bool = False
    PS2: bool = False
    PS3: bool = False
    PS4: bool = False
    PM1: bool = False
    PM2: bool = False
    PM3: bool = False
    PM4: bool = False
    PM5: bool = False
    PM6: bool = False
    PP1: bool = False
    PP2: bool = False
    PP3: bool = False
    PP4: bool = False
    PP5: bool = False
    BA1: bool = False
    BS1: bool = False
    BS2: bool = False
    BS3: bool = False
    BS4: bool = False
    BP1: bool = False
    BP2: bool = False
    BP3: bool = False
    BP4: bool = False
    BP5: bool = False
    BP6: bool = False
    BP7: bool = False


class PS1PM5(BaseModel):
    """PS1 and PM5 criteria prediction."""

    PS1: bool = False
    PM5: bool = False


class PM4BP3(BaseModel):
    """PM4 and BP3 criteria prediction."""

    PM4: bool = False
    BP3: bool = False


class BA1BS1BS2PM2(BaseModel):
    """BA1, BS1, BS2, and PM2 criteria prediction."""

    BA1: bool = False
    BS1: bool = False
    BS2: bool = False
    PM2: bool = False


class PM1(BaseModel):
    """PM1 criteria prediction."""

    PM1: bool = False


class PP2BP1(BaseModel):
    """PP2 and BP6 criteria prediction."""

    PP2: bool = False
    BP1: bool = False


class BP7(BaseModel):
    """BP7 criteria prediction."""

    BP7: bool = False


class PP3BP4(BaseModel):
    """PP3 and BP4 criteria prediction."""

    PP3: bool = False
    BP4: bool = False


class ACMGPrediction(AutoAcmgBaseEnum):
    """ACMG prediction enumeration."""

    NotSet = auto()
    NotApplicable = auto()
    NotAutomated = auto()
    Depricated = auto()
    Positive = auto()
    Negative = auto()


class CriteriaPrediction(BaseModel):
    """Criteria prediction."""

    name: str
    prediction: ACMGPrediction = ACMGPrediction.NotSet
    summary: str = ""
    description: str = ""


class AutoACMGResult(BaseModel):
    """Response of the ACMG criteria prediction."""

    pvs1: CriteriaPrediction = CriteriaPrediction(name="PVS1")
    ps1: CriteriaPrediction = CriteriaPrediction(name="PS1")
    ps2: CriteriaPrediction = CriteriaPrediction(name="PS2", prediction=ACMGPrediction.NotAutomated)
    ps3: CriteriaPrediction = CriteriaPrediction(name="PS3", prediction=ACMGPrediction.NotAutomated)
    ps4: CriteriaPrediction = CriteriaPrediction(name="PS4", prediction=ACMGPrediction.NotAutomated)
    pm1: CriteriaPrediction = CriteriaPrediction(name="PM1")
    pm2: CriteriaPrediction = CriteriaPrediction(name="PM2")
    pm3: CriteriaPrediction = CriteriaPrediction(name="PM3", prediction=ACMGPrediction.NotAutomated)
    pm4: CriteriaPrediction = CriteriaPrediction(name="PM4")
    pm5: CriteriaPrediction = CriteriaPrediction(name="PM5")
    pm6: CriteriaPrediction = CriteriaPrediction(name="PM6", prediction=ACMGPrediction.NotAutomated)
    pp1: CriteriaPrediction = CriteriaPrediction(name="PP1", prediction=ACMGPrediction.NotAutomated)
    pp2: CriteriaPrediction = CriteriaPrediction(name="PP2")
    pp3: CriteriaPrediction = CriteriaPrediction(name="PP3")
    pp4: CriteriaPrediction = CriteriaPrediction(name="PP4", prediction=ACMGPrediction.NotAutomated)
    pp5: CriteriaPrediction = CriteriaPrediction(name="PP5", prediction=ACMGPrediction.Depricated)
    ba1: CriteriaPrediction = CriteriaPrediction(name="BA1")
    bs1: CriteriaPrediction = CriteriaPrediction(name="BS1")
    bs2: CriteriaPrediction = CriteriaPrediction(name="BS2")
    bs3: CriteriaPrediction = CriteriaPrediction(name="BS3", prediction=ACMGPrediction.NotAutomated)
    bs4: CriteriaPrediction = CriteriaPrediction(name="BS4", prediction=ACMGPrediction.NotAutomated)
    bp1: CriteriaPrediction = CriteriaPrediction(name="BP1")
    bp2: CriteriaPrediction = CriteriaPrediction(name="BP2", prediction=ACMGPrediction.NotAutomated)
    bp3: CriteriaPrediction = CriteriaPrediction(name="BP3")
    bp4: CriteriaPrediction = CriteriaPrediction(name="BP4")
    bp5: CriteriaPrediction = CriteriaPrediction(name="BP5", prediction=ACMGPrediction.NotAutomated)
    bp6: CriteriaPrediction = CriteriaPrediction(name="BP6", prediction=ACMGPrediction.Depricated)
    bp7: CriteriaPrediction = CriteriaPrediction(name="BP7")
