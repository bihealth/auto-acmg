"""Custom exceptions for the project."""


class AutoAcmgBaseException(Exception):
    """Base class for exceptions."""


class ApiCallException(AutoAcmgBaseException):
    pass


class AnnonarsException(ApiCallException):
    pass


class MehariException(ApiCallException):
    pass


class ParseError(AutoAcmgBaseException):
    """Exception for parsing errors."""

    pass


class InvalidGenomeBuild(AutoAcmgBaseException):
    """Exception for invalid genome build."""

    pass


class InvalidPos(AutoAcmgBaseException):
    """Exception for invalid position."""

    pass


class MappingError(AutoAcmgBaseException):
    """Exception for errors in mapping chromosomes or genome builds."""

    pass


class InvalidAPIResposeError(AutoAcmgBaseException):
    """Exception for API errors."""

    pass


class AlgorithmError(AutoAcmgBaseException):
    """Exception for algorithm errors."""

    pass


class AutoPVS1Error(AutoAcmgBaseException):
    """Exception for AutoPVS1 errors."""
