"""Main entry point for the autopvs1 CLI."""

import argparse
import asyncio
import logging
import sys
from typing import Optional

from src.genome_builds import GenomeRelease
from src.seqvar import SeqVar, SeqVarResolver


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


async def main(args: Optional[list[str]] = None):
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

    try:
        seqvar_resolver = SeqVarResolver()
        seqvar: SeqVar = await seqvar_resolver.resolve_seqvar(
            parsed_args.variant, parsed_args.genome_release
        )
        logging.info(f"Resolved variant: {seqvar}. Dictionary representation: {seqvar.__dict__}")
    except Exception as e:
        logging.error(e)


if __name__ == "__main__":
    # Run the main function
    asyncio.run(main())
