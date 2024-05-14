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
    pass


class InvalidGenomeBuild(AutoAcmgBaseException):
    pass


class InvalidPos(AutoAcmgBaseException):
    pass


class MappingError(AutoAcmgBaseException):
    pass


class InvalidAPIResposeError(AutoAcmgBaseException):
    pass


class AlgorithmError(AutoAcmgBaseException):
    pass


class AutoPVS1Error(AutoAcmgBaseException):
    pass


class MissingDataError(AutoAcmgBaseException):
    pass
