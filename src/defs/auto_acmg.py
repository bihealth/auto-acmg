from enum import Enum, EnumMeta
from typing import Any

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
