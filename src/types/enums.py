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
