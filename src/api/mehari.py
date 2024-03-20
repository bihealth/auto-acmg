"""Mehari API client."""

import aiohttp

from src.core.config import settings
from src.genome_builds import GenomeRelease
from src.models.mehari_gene import GeneTranscripts
from src.models.mehari_seqvar import TranscriptsSeqVar
from src.seqvar import SeqVar

#: Mehari API base URL
MEHARI_API_BASE_URL = f"{settings.API_REEV_URL}/mehari"


class MehariClient:
    def __init__(self, api_base_url: str = MEHARI_API_BASE_URL):
        self.api_base_url = api_base_url

    async def get_seqvar_transcripts(self, seqvar: SeqVar) -> TranscriptsSeqVar | None:
        """
        Get transcripts for a sequence variant.

        :param seqvar: Sequence variant
        :type seqvar: SeqVar
        :return: Transcripts
        :rtype: TranscriptsSeqVar | None
        """
        url = (
            f"{self.api_base_url}/seqvars/csq?"
            f"genome_release={seqvar.genome_release.name.lower()}"
            f"&chromosome={seqvar.chrom}"
            f"&position={seqvar.pos}"
            f"&reference={seqvar.insert}"
            f"&alternative={seqvar.delete}"
        )
        async with aiohttp.ClientSession() as session:
            async with session.get(url) as response:
                if response.status == 200:
                    return await response.json()
                else:
                    return None

    async def get_gene_transcripts(
        self, hgnc_id: str, genome_build: GenomeRelease
    ) -> GeneTranscripts | None:
        """ "
        Get transcripts for a gene.

        :param hgnc_id: HGNC gene ID
        :type hgnc_id: str
        :param genome_build: Genome build
        :type genome_build: GenomeRelease
        :return: Transcripts
        :rtype: GeneTranscripts | None
        """
        genome_build_mapping = {
            GenomeRelease.GRCh37: "GENOME_BUILD_GRCH37",
            GenomeRelease.GRCh38: "GENOME_BUILD_GRCH38",
        }
        url = (
            f"{self.api_base_url}/genes/txs?"
            f"hgncId={hgnc_id}"
            f"&genomeBuild={genome_build_mapping[genome_build]}"
        )
        async with aiohttp.ClientSession() as session:
            async with session.get(url) as response:
                if response.status == 200:
                    return await response.json()
                else:
                    return None
