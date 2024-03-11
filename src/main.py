"""Main entry point for the autopvs1 CLI."""

import argparse
import logging
import sys
from enum import Enum, auto
from typing import Optional


class GenomeRelease(Enum):
    """Enumeration for allowed genome release values."""

    hg19 = auto()
    hg38 = auto()
    GRCh37 = auto()
    GRCh38 = auto()

    @staticmethod
    def from_string(value: str):
        """Converts string to enum member if possible, otherwise returns None."""
        for member in GenomeRelease:
            if member.name == value:
                return member
        return None

    @staticmethod
    def list():
        """Returns list of enum member names."""
        return list(map(lambda c: c.name, GenomeRelease))


class GenomeReleaseAction(argparse.Action):
    """Custom action to convert genome release string to enum member."""

    def __call__(self, parser, namespace, values, option_string=None):
        genome_release = GenomeRelease.from_string(values)
        if genome_release is None:
            logging.warning(
                f"Invalid genome release value: '{values}'. It must be one of {', '.join(GenomeRelease.list())}."
            )
            setattr(namespace, self.dest, None)
        else:
            setattr(namespace, self.dest, genome_release)


def create_parser():
    """Creates and returns the argparse parser with all the defined arguments."""
    parser = argparse.ArgumentParser(description="Entry point for the autopvs1 CLI.")
    parser.add_argument(
        "variant", type=str, help='Variant to be classified, e.g., "NM_000038.3:c.797G>A".'
    )
    parser.add_argument(
        "--genome_release",
        action=GenomeReleaseAction,
        help=f'Genome release, e.g., {", ".join(GenomeRelease.list())}.',
    )

    return parser


def main(args: Optional[list[str]] = None):
    """Entry point for the CLI."""
    if args is None:
        args = sys.argv[1:]

    # Setup logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    parser = create_parser()
    parsed_args = parser.parse_args(args)

    logging.info(f"Variant to be classified: {parsed_args.variant}")
    if parsed_args.genome_release:
        logging.info(f"Genome release: {parsed_args.genome_release.name}")
    else:
        logging.info("No valid genome release specified or no genome release provided.")


if __name__ == "__main__":
    main()
