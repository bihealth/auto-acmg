"""Implementation of sequence variant class."""

import re
from typing import Optional

from pydantic import BaseModel

from src.api.dotty import DottyClient
from src.core.config import Config
from src.defs.exceptions import AutoAcmgBaseException, InvalidPos, ParseError
from src.defs.genome_builds import (
    CHROM_LENGTHS_37,
    CHROM_LENGTHS_38,
    REFSEQ_CHROM_37,
    REFSEQ_CHROM_38,
    GenomeRelease,
)

#: Regular expression for gnomAD-style variant representation
REGEX_GNOMAD_VARIANT = re.compile(
    r"^(?:(?P<genome_build>(?:\w+))-)?(?P<chrom>(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT))-(?P<pos>\d+)-(?P<delete>[ACGT]+)-(?P<insert>[ACGT]+)$",
    re.IGNORECASE,
)
#: Regular expression for canonical SPDI variant representation
REGEX_CANONICAL_SPDI = re.compile(
    r"^(?P<sequence>NC_(?:\d{6}\.\d+)):(?P<pos>\d+):(?P<delete>[ACGT]+):(?P<insert>[ACGT]+)$", re.I
)
#: Regular expression for "relaxed" SPDI variant representation
REGEX_RELAXED_SPDI = re.compile(
    r"^(?:(?P<genome_build>\w+):)?(?P<chrom>(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)):(?P<pos>\d+):(?P<delete>[ACGT]+):(?P<insert>[ACGT]+)$",
    re.IGNORECASE,
)
#: Regular expression for dbSNP
REGEX_DBSNP_ID = re.compile(r"^rs\d+$", re.IGNORECASE)
#: Regular expression for ClinVar
REGEX_CLINVAR_ID = re.compile(
    r"^(?P<accession>(?:RCV|VCV)\d{9})(?:\.(?P<version>\d+))?$", re.IGNORECASE
)


class SeqVar:
    """A class to represent a sequence variant."""

    def __init__(
        self,
        genome_release: GenomeRelease,
        chrom: str,
        pos: int,
        delete: str,
        insert: str,
        user_repr: Optional[str] = None,
    ):
        self.genome_release = genome_release
        self.chrom = self._normalize_chromosome(chrom)
        self.pos = pos
        self.delete = delete.upper()
        self.insert = insert.upper()
        self.user_repr = (
            user_repr
            if user_repr
            else f"{genome_release.name}-{self.chrom}-{pos}-{delete}-{insert}"
        )

    @staticmethod
    def _normalize_chromosome(chrom: str) -> str:
        """Normalize the chromosome name."""
        return chrom.lower().replace("chr", "").replace("m", "mt").replace("mtt", "mt").upper()

    def __repr__(self):
        """Return a user-friendly representation of the variant."""
        return self.user_repr

    def _as_dict(self):
        """Return a dictionary representation of the variant."""
        return {
            "genome_release": self.genome_release,
            "chrom": self.chrom,
            "pos": self.pos,
            "delete": self.delete,
            "insert": self.insert,
        }

    def __eq__(self, other):
        """Return True if the two objects are equal."""
        if not isinstance(other, SeqVar):
            return False
        return self._as_dict() == other._as_dict()


class SeqVarResolver:
    """The class to resolve sequence variants."""

    def __init__(self, *, config: Optional[Config] = None):
        self.config = config or Config()
        self.dotty_client = DottyClient(api_base_url=self.config.api_base_url_dotty)

    def _validate_seqvar(self, variant: SeqVar) -> SeqVar:
        """Validate the sequence variant position.

        Args:
            variant (SeqVar): Sequence variant

        Returns:
            SeqVar: Sequence variant

        Raises:
            InvalidPos: If the position is invalid
        """
        if variant.pos < 1:
            raise InvalidPos(f"Invalid position: {variant.pos}")

        stop_pos = variant.pos + len(variant.delete) - 1
        chrom_lengths = (
            CHROM_LENGTHS_37 if variant.genome_release == GenomeRelease.GRCh37 else CHROM_LENGTHS_38
        )
        if stop_pos > chrom_lengths.get(variant.chrom, 0):
            raise InvalidPos(f"Invalid position: {variant.pos}")

        return variant

    def _normalize_chrom(self, value: str) -> str:
        """Normalize the chromosome name.

        Replaces 'chr' with '' and 'm' with 'mt'.
        """
        return value.lower().replace("chr", "").replace("m", "mt").replace("mtt", "mt").upper()

    def _parse_separated_seqvar(
        self, value: str, default_genome_release: GenomeRelease = GenomeRelease.GRCh38
    ) -> SeqVar:
        """Parse a colon/hyphen separated sequence variant representation.

        Args:
            value (str): Sequence variant representation
            default_genome_release (GenomeRelease): Default genome release

        Returns:
            SeqVar: Sequence variant

        Raises:
            ParseError: If the variant representation is invalid
        """
        match = REGEX_GNOMAD_VARIANT.match(value) or REGEX_RELAXED_SPDI.match(value)
        if not match:
            raise ParseError(f"Unable to parse colon/hyphen separated seqvar: {value}")

        genome_build_value = match.group("genome_build")
        genome_build = (
            GenomeRelease[genome_build_value] if genome_build_value else default_genome_release
        )
        chrom = self._normalize_chrom(match.group("chrom"))
        pos = int(match.group("pos"))
        delete = match.group("delete").upper()
        insert = match.group("insert").upper()

        variant = SeqVar(
            genome_release=genome_build,
            chrom=chrom,
            pos=pos,
            delete=delete,
            insert=insert,
            user_repr=value,
        )
        return self._validate_seqvar(variant)

    def _parse_canonical_spdi_seqvar(self, value: str) -> SeqVar:
        """Parse a canonical SPDI sequence variant representation.

        Args:
            value (str): Sequence variant representation

        Returns:
            SeqVar: Sequence variant

        Raises:
            ParseError: If the variant representation is invalid
        """
        match = REGEX_CANONICAL_SPDI.match(value)
        if not match:
            raise ParseError(f"Unable to parse canonical SPDI variant: {value}")

        sequence = match.group("sequence").upper()
        pos = int(match.group("pos")) + 1  # SPDI is 0-based
        delete = match.group("delete").upper()
        insert = match.group("insert").upper()

        if sequence in REFSEQ_CHROM_37:
            genome_build = GenomeRelease.GRCh37
            chrom = REFSEQ_CHROM_37[sequence]
        elif sequence in REFSEQ_CHROM_38:
            genome_build = GenomeRelease.GRCh38
            chrom = REFSEQ_CHROM_38[sequence]
        else:
            raise ParseError(f"Unknown sequence: {sequence}")

        variant = SeqVar(
            genome_release=genome_build,
            chrom=chrom,
            pos=pos,
            delete=delete,
            insert=insert,
            user_repr=value,
        )
        return self._validate_seqvar(variant)

    def resolve_seqvar(self, value: str, genome_release: GenomeRelease) -> SeqVar:
        """Resolve a sequence variant.

        Supports gnomAD-style, SPDI and dbSNP representations.
        ClinVar IDs are not supported at the moment.

        Args:
            value (str): Sequence variant representation
            genome_release (GenomeRelease): Genome release version

        Returns:
            SeqVar: Sequence variant

        Raises:
            ParseError: If the variant representation is invalid or cannot be resolved
        """
        try:
            return self._parse_separated_seqvar(value, default_genome_release=genome_release)
        except ParseError:
            pass
        except InvalidPos as e:
            raise ParseError(f"Invalid position: {e}")

        try:
            return self._parse_canonical_spdi_seqvar(value)
        except ParseError:
            pass
        except InvalidPos as e:
            raise ParseError(f"Invalid position: {e}")

        try:
            spdi = self.dotty_client.to_spdi(value, assembly=genome_release)
            if spdi is not None and spdi.success and spdi.value is not None:
                return SeqVar(
                    genome_release=GenomeRelease.from_string(spdi.value.assembly),
                    chrom=spdi.value.contig,
                    pos=spdi.value.pos,
                    delete=spdi.value.reference_deleted,
                    insert=spdi.value.alternate_inserted,
                    user_repr=value,
                )
            else:
                raise ParseError(f"Unable to resolve seqvar: {value}")
        except AutoAcmgBaseException as e:
            raise ParseError(f"Unable to resolve seqvar. The error was: {e}")
