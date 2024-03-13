"""Genome build information and operations."""

from enum import Enum, auto
from typing import Literal

from src.core.exceptions import MappingError


#: GenomeRelease enume
class GenomeRelease(Enum):
    """Enumeration for allowed genome release values."""

    GRCh37 = auto()
    GRCh38 = auto()

    @staticmethod
    def from_string(value: str):
        """Converts string to enum member if possible, otherwise returns None."""
        genome_mapping = {"hg19": "GRCh37", "hg38": "GRCh38"}
        value_mapped = genome_mapping.get(value, value)
        for member in GenomeRelease:
            if member.name == value_mapped:
                return member
        return None

    @staticmethod
    def list():
        """Returns list of enum member names."""
        return list(map(lambda c: c.name, GenomeRelease))


#: Mapping from chromosome number to RefSeq accession (GRCh37)
CHROM_REFSEQ_37: dict[str, str] = {
    "1": "NC_000001.10",
    "2": "NC_000002.11",
    "3": "NC_000003.11",
    "4": "NC_000004.11",
    "5": "NC_000005.9",
    "6": "NC_000006.11",
    "7": "NC_000007.13",
    "8": "NC_000008.10",
    "9": "NC_000009.11",
    "10": "NC_000010.10",
    "11": "NC_000011.9",
    "12": "NC_000012.11",
    "13": "NC_000013.10",
    "14": "NC_000014.8",
    "15": "NC_000015.9",
    "16": "NC_000016.9",
    "17": "NC_000017.10",
    "18": "NC_000018.9",
    "19": "NC_000019.9",
    "20": "NC_000020.10",
    "21": "NC_000021.8",
    "22": "NC_000022.10",
    "X": "NC_000023.10",
    "Y": "NC_000024.9",
    "MT": "NC_012920.1",
}

#: Mapping from RefSeq accession to chromosome number (GRCh37)
REFSEQ_CHROM_37: dict[str, str] = {v: k for k, v in CHROM_REFSEQ_37.items()}

#: Mapping from chromosome number to chromosome length (GRCh37)
CHROM_LENGTHS_37: dict[str, int] = {
    "1": 249250621,
    "2": 243199373,
    "3": 198022430,
    "4": 191154276,
    "5": 180915260,
    "6": 171115067,
    "7": 159138663,
    "8": 146364022,
    "9": 141213431,
    "10": 135534747,
    "11": 135006516,
    "12": 133851895,
    "13": 115169878,
    "14": 107349540,
    "15": 102531392,
    "16": 90354753,
    "17": 81195210,
    "18": 78077248,
    "19": 59128983,
    "20": 63025520,
    "21": 48129895,
    "22": 51304566,
    "X": 155270560,
    "Y": 59373566,
    "MT": 16569,
}

#: Mapping from chromosome number to RefSeq accession (GRCh38)
CHROM_REFSEQ_38: dict[str, str] = {
    "1": "NC_000001.11",
    "2": "NC_000002.12",
    "3": "NC_000003.12",
    "4": "NC_000004.12",
    "5": "NC_000005.10",
    "6": "NC_000006.12",
    "7": "NC_000007.14",
    "8": "NC_000008.11",
    "9": "NC_000009.12",
    "10": "NC_000010.11",
    "11": "NC_000011.10",
    "12": "NC_000012.12",
    "13": "NC_000013.11",
    "14": "NC_000014.9",
    "15": "NC_000015.10",
    "16": "NC_000016.10",
    "17": "NC_000017.11",
    "18": "NC_000018.10",
    "19": "NC_000019.10",
    "20": "NC_000020.11",
    "21": "NC_000021.9",
    "22": "NC_000022.11",
    "X": "NC_000023.11",
    "Y": "NC_000024.10",
    "MT": "NC_012920.1",
}

#: Mapping from RefSeq accession to chromosome number (GRCh38)
REFSEQ_CHROM_38: dict[str, str] = {v: k for k, v in CHROM_REFSEQ_38.items()}

#: Mapping from chromosome number to chromosome length (GRCh38)
CHROM_LENGTHS_38: dict[str, int] = {
    "1": 248956422,
    "2": 242193529,
    "3": 198295559,
    "4": 190214555,
    "5": 181538259,
    "6": 170805979,
    "7": 159345973,
    "8": 145138636,
    "9": 138394717,
    "10": 133797422,
    "11": 135086622,
    "12": 133275309,
    "13": 114364328,
    "14": 107043718,
    "15": 101991189,
    "16": 90338345,
    "17": 83257441,
    "18": 80373285,
    "19": 58617616,
    "20": 64444167,
    "21": 46709983,
    "22": 50818468,
    "X": 156040895,
    "Y": 57227415,
    "MT": 16569,
}


def refseq_to_genome_build(refseq_acc: str) -> GenomeRelease:
    """
    Determines the genome build based on RefSeq accession.

    :param refseq_acc: RefSeq accession
    :type refseq_acc: str
    :return: Genome build
    :rtype: GenomeRelease
    :raises MappingError: If RefSeq accession is unknown
    """
    refseq_acc_upper = refseq_acc.upper()
    if refseq_acc_upper in CHROM_REFSEQ_37.values():
        return GenomeRelease.GRCh37
    elif refseq_acc_upper in CHROM_REFSEQ_38.values():
        return GenomeRelease.GRCh38
    else:
        raise MappingError(f"Unknown RefSeq identifier: {refseq_acc}")
