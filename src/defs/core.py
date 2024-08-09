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


class AutoAcmgBaseModel(BaseModel):
    """Base model for ACMG criteria."""

    # class Config:
    #     use_enum_values = True
    #     arbitrary_types_allowed = True
