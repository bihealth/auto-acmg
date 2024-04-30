"""Custom exceptions for the project."""


class ParseError(Exception):
    """Exception for parsing errors."""

    pass


class InvalidGenomeBuild(Exception):
    """Exception for invalid genome build."""

    pass


class InvalidPos(Exception):
    """Exception for invalid position."""

    pass


class MappingError(Exception):
    """Exception for errors in mapping chromosomes or genome builds."""

    pass


class InvalidAPIResposeError(Exception):
    """Exception for API errors."""

    pass


class AlgorithmError(Exception):
    """Exception for algorithm errors."""

    pass
