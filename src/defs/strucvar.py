"""Implementation of structural variant class."""

import re
from enum import Enum, auto
from typing import Optional, Union

from src.defs.exceptions import InvalidPos, ParseError
from src.defs.genome_builds import CHROM_LENGTHS_37, CHROM_LENGTHS_38, GenomeRelease

#: Regular expression for colon-separated SV representation
REGEX_CNV_COLON = re.compile(
    r"^(?P<sv_type>DEL|DUP):(?:(?:(?P<genome_build>\w+):)?(?P<chrom>(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT))):(?P<start>\d+):(?P<stop>\d+)$",
    re.IGNORECASE,
)

#: Regular expression for hyphen-separated SV representation
REGEX_CNV_HYPHEN = re.compile(
    r"^(?P<sv_type>DEL|DUP)-(?:(?:(?P<genome_build>\w+)-)?(?P<chrom>(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)))-(?P<start>\d+)-(?P<stop>\d+)$",
    re.IGNORECASE,
)


class StrucVarType(Enum):
    """Enumeration for structural variant type."""

    DEL = auto()
    DUP = auto()
    # INV = auto()
    # INS = auto()
    # TRA = auto()
    # BND = auto()


class StrucVar:
    """A class to represent a structural variant."""

    def __init__(
        self,
        sv_type: StrucVarType,
        genome_release: GenomeRelease,
        chrom: str,
        start: int,
        stop: int,
        user_repr: Optional[str] = None,
    ):
        self.sv_type = sv_type
        self.genome_release = genome_release
        self.chrom = self._normalize_chromosome(chrom)
        self.start = start
        self.stop = stop
        self.user_repr = (
            user_repr if user_repr else f"{sv_type}-{genome_release.name}-{chrom}-{start}-{stop}"
        )

    def _normalize_chromosome(self, chrom: str) -> str:
        """Normalize the chromosome name."""
        return chrom.lower().replace("chr", "").replace("m", "mt").upper()

    def __repr__(self):
        return self.user_repr


class StrucVarResolver:
    """The class to resolve structural variant representations."""

    def __init__(self):
        pass

    def _validate_strucvar(self, variant: StrucVar) -> StrucVar:
        """
        Validate the structural variant position.

        :param variant: Structural variant
        :type variant: StrucVar
        :return: Structural variant
        :rtype: StrucVar
        :raises InvalidPos: If the structural variant position is invalid
        """
        if variant.start > variant.stop or variant.start < 1:
            raise InvalidPos(f"Invalid positions: start={variant.start}, stop={variant.stop}")

        chrom_lengths = (
            CHROM_LENGTHS_37 if variant.genome_release == GenomeRelease.GRCh37 else CHROM_LENGTHS_38
        )
        if variant.stop > chrom_lengths.get(variant.chrom, 0):
            raise InvalidPos(
                f"Invalid stop position: {variant.stop}. Chromosome length: {chrom_lengths.get(variant.chrom, 0)}"
            )

        return variant

    def _parse_separated_strucvar(
        self, value: str, default_genome_release: GenomeRelease = GenomeRelease.GRCh38
    ) -> StrucVar:
        """
        Parse a separated structural variant representation.

        :param value: Structural variant representation
        :type value: str
        :param default_genome_release: Default genome release
        :type default_genome_release: GenomeRelease
        :return: Structural variant
        :rtype: StrucVar
        :raises ParseError: If the structural variant representation is invalid
        """
        match_colon = REGEX_CNV_COLON.match(value)
        match_hyphen = REGEX_CNV_HYPHEN.match(value)
        match = match_colon if match_colon else match_hyphen
        if not match:
            raise ParseError(f"Unable to parse colon/hyphen separated strucvar: {value}")

        sv_type = StrucVarType[match.group("sv_type").upper()]
        genome_release_value = match.group("genome_build")
        genome_release = (
            GenomeRelease[genome_release_value] if genome_release_value else default_genome_release
        )
        chrom = match.group("chrom")
        start = int(match.group("start"))
        stop = int(match.group("stop"))

        variant = StrucVar(sv_type, genome_release, chrom, start, stop, user_repr=value)
        return self._validate_strucvar(variant)

    def resolve_strucvar(self, value: str, genome_release: GenomeRelease) -> StrucVar:
        """
        Resolve the structural variant representation.

        :param value: Structural variant representation
        :type value: str
        :param genome_release: Genome release
        :type genome_release: GenomeRelease
        :return: Structural variant
        :rtype: StrucVar
        :raises ParseError: If the structural variant representation is invalid
        """
        try:
            return self._parse_separated_strucvar(value, genome_release)
        except Exception as e:
            raise ParseError(f"Unable to parse structural variant: {value}") from e