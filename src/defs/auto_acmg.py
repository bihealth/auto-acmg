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
    """ACMG prediction enumeration."""

    NotSet = auto()
    NotApplicable = auto()
    NotAutomated = auto()
    Depricated = auto()
    Met = auto()
    Unmet = auto()


class AutoACMGCriteria(BaseModel):
    """Criteria prediction."""

    name: str
    prediction: AutoACMGPrediction = AutoACMGPrediction.NotSet
    summary: str = ""
    description: str = ""


class AutoACMGResult(BaseModel):
    """Response of the ACMG criteria prediction."""

    pvs1: AutoACMGCriteria = AutoACMGCriteria(name="PVS1")
    ps1: AutoACMGCriteria = AutoACMGCriteria(name="PS1")
    ps2: AutoACMGCriteria = AutoACMGCriteria(name="PS2", prediction=AutoACMGPrediction.NotAutomated)
    ps3: AutoACMGCriteria = AutoACMGCriteria(name="PS3", prediction=AutoACMGPrediction.NotAutomated)
    ps4: AutoACMGCriteria = AutoACMGCriteria(name="PS4", prediction=AutoACMGPrediction.NotAutomated)
    pm1: AutoACMGCriteria = AutoACMGCriteria(name="PM1")
    pm2: AutoACMGCriteria = AutoACMGCriteria(name="PM2")
    pm3: AutoACMGCriteria = AutoACMGCriteria(name="PM3", prediction=AutoACMGPrediction.NotAutomated)
    pm4: AutoACMGCriteria = AutoACMGCriteria(name="PM4")
    pm5: AutoACMGCriteria = AutoACMGCriteria(name="PM5")
    pm6: AutoACMGCriteria = AutoACMGCriteria(name="PM6", prediction=AutoACMGPrediction.NotAutomated)
    pp1: AutoACMGCriteria = AutoACMGCriteria(name="PP1", prediction=AutoACMGPrediction.NotAutomated)
    pp2: AutoACMGCriteria = AutoACMGCriteria(name="PP2")
    pp3: AutoACMGCriteria = AutoACMGCriteria(name="PP3")
    pp4: AutoACMGCriteria = AutoACMGCriteria(name="PP4", prediction=AutoACMGPrediction.NotAutomated)
    pp5: AutoACMGCriteria = AutoACMGCriteria(name="PP5", prediction=AutoACMGPrediction.Depricated)
    ba1: AutoACMGCriteria = AutoACMGCriteria(name="BA1")
    bs1: AutoACMGCriteria = AutoACMGCriteria(name="BS1")
    bs2: AutoACMGCriteria = AutoACMGCriteria(name="BS2")
    bs3: AutoACMGCriteria = AutoACMGCriteria(name="BS3", prediction=AutoACMGPrediction.NotAutomated)
    bs4: AutoACMGCriteria = AutoACMGCriteria(name="BS4", prediction=AutoACMGPrediction.NotAutomated)
    bp1: AutoACMGCriteria = AutoACMGCriteria(name="BP1")
    bp2: AutoACMGCriteria = AutoACMGCriteria(name="BP2", prediction=AutoACMGPrediction.NotAutomated)
    bp3: AutoACMGCriteria = AutoACMGCriteria(name="BP3")
    bp4: AutoACMGCriteria = AutoACMGCriteria(name="BP4")
    bp5: AutoACMGCriteria = AutoACMGCriteria(name="BP5", prediction=AutoACMGPrediction.NotAutomated)
    bp6: AutoACMGCriteria = AutoACMGCriteria(name="BP6", prediction=AutoACMGPrediction.Depricated)
    bp7: AutoACMGCriteria = AutoACMGCriteria(name="BP7")

    def export(self, filename: Optional[str] = None) -> None:
        """Export the ACMG criteria prediction to a JSON file."""
        data = self.model_dump_json(indent=4)
        if not filename:
            filename = f"acmg_criteria_prediction.json"
        with open(filename, "w") as file:
            for line in data.splitlines():
                file.write(line)
