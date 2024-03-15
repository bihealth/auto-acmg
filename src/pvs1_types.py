"""Enum classes for PVS1."""

from enum import Enum, auto


class SeqVarConsequence(Enum):
    """Consequence of a sequence variant."""

    NonsenseFrameshift = auto()
    SpliceSites = auto()
    InitiationCodon = auto()


class PVS1Prediction(Enum):
    """PVS1 prediction."""

    PVS1 = auto()
    PVS1_Strong = auto()
    PVS1_Moderate = auto()
    PVS1_Supporting = auto()
    NotPVS1 = auto()
