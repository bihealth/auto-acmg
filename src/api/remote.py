"""API client for remote services."""

from typing import Any

import requests
from loguru import logger

from src.defs.exceptions import ApiCallException
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


class RemoteClient:
    """API client for remote services."""

    def get_repeat_masked_regions(self, seqvar: SeqVar) -> Any:
        """Get repeat masked regions for a sequence variant.

        Args:
            seqvar (SeqVar): The sequence variant.

        Returns:
            Any: The repeat masked regions.
        """
        url = "https://genome.ucsc.edu/cgi-bin/hgTables"
        genome_release = "hg38" if seqvar.genome_release == GenomeRelease.GRCh38 else "hg19"
        position = f"chr{seqvar.chrom}:{seqvar.pos - 50}-{seqvar.pos + 50}"
        params = {
            "clade": "mammal",
            "org": "Human",
            "db": genome_release,
            "hgta_group": "rep",
            "hgta_track": "rmsk",
            "hgta_table": "rmsk",
            "hgta_regionType": "genome",
            "position": position,
            "hgta_outputType": "bed",
            # "hgta_outFileName": "repeat_masked_regions.csv",
            "hgta_outSep": "tab",
            "hgta_compressType": "none",
        }
        logger.debug("POST request to: {}", url)
        response = requests.post(url, data=params)
        try:
            response.raise_for_status()
            return response.text
        except requests.RequestException as e:
            logger.exception("Request failed: {}", e)
            raise ApiCallException("Failed to get repeat masked regions.") from e


#         hgsid: 2298522188_l1mV38Ng7jIGaQrtKXq02A1TEEvP
# jsh_pageVertPos: 0
# clade: mammal
# org: Human
# db: hg38
# hgta_group: rep
# hgta_track: rmsk
# hgta_table: rmsk
# hgta_regionType: genome
# position: chr7:155,799,529-155,799,999
# hgta_outputType: bed
# boolshad.sendToGalaxy: 0
# boolshad.sendToGreat: 0
# hgta_outFileName: ddd
# hgta_outSep: tab
# hgta_compressType: none
# hgta_doTopSubmit: get output
