"""Main entry point for the autopvs1 CLI."""

import argparse
import asyncio
import logging
import sys
from typing import Optional

from src.genome_builds import GenomeRelease
from src.pvs1 import PVS1
from src.seqvar import SeqVar, SeqVarResolver


class GenomeReleaseAction(argparse.Action):
    """Custom action to convert genome release string to enum member."""

    def __call__(self, _parser, namespace, values, _option_string=None):
        genome_release = GenomeRelease.from_string(values)
        if genome_release is None:
            logging.warning(
                f"Invalid genome release value: '{values}'. It must be one of {', '.join(GenomeRelease.list())}."
            )
            setattr(namespace, self.dest, None)
        else:
            setattr(namespace, self.dest, genome_release)


def create_parser():
    """Helper function for creating the argparse parser with all the defined arguments."""
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
    logger = logging.getLogger(__name__)

    parser = create_parser()
    parsed_args = parser.parse_args(args)

    logger.info(f"Variant to be classified: {parsed_args.variant}")
    if parsed_args.genome_release:
        logger.info(f"Genome release: {parsed_args.genome_release.name}")
    else:
        logger.info("No valid genome release specified or no genome release provided.")

    try:
        seqvar_resolver = SeqVarResolver()
        seqvar: SeqVar = await seqvar_resolver.resolve_seqvar(
            parsed_args.variant, parsed_args.genome_release
        )
        logger.info(f"Resolved variant: {seqvar}. Dictionary representation: {seqvar.__dict__}")
        pvs1 = PVS1(seqvar)
        await pvs1.run()
        logger.info(f"PVS1 transcripts: {pvs1.seqvar_transcripts}\n\n\n{pvs1.gene_transcripts}")
    except Exception as e:
        logger.error(e)


if __name__ == "__main__":
    # Run the main function
    asyncio.run(main())
