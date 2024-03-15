"""Main entry point for the autopvs1 CLI."""

import argparse
import asyncio
import logging
import sys
from typing import Optional

from src.autoPVS1 import AutoPVS1
from src.core.config import settings
from src.genome_builds import GenomeRelease

# Setup logging
logging_level = logging.DEBUG if settings.DEBUG else logging.INFO
logging.basicConfig(level=logging_level, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


class GenomeReleaseAction(argparse.Action):
    """Custom action to convert genome release string to enum member."""

    def __call__(self, _parser, namespace, values, _option_string=None):
        genome_release = GenomeRelease.from_string(values)
        if genome_release is None:
            logger.warning(
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

    parser = create_parser()
    parsed_args = parser.parse_args(args)

    logger.info(f"Variant to be classified: {parsed_args.variant}")
    if parsed_args.genome_release:
        logger.info(f"Genome release: {parsed_args.genome_release.name}")
    else:
        logger.info("No valid genome release specified or no genome release provided.")

    try:
        auto_pvs1 = AutoPVS1(parsed_args.variant, parsed_args.genome_release)
        await auto_pvs1.predict()
    except Exception as e:
        logger.error(e)


if __name__ == "__main__":
    # Run the main function
    asyncio.run(main())
