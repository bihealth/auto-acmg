"""Enum classes for PVS1."""

from enum import Enum, auto


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
