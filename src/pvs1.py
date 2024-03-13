"""Implementations of the PVS1 algorithm."""

import logging
import random
import string

from src.api.mehari import MehariClient
from src.seqvar import SeqVar

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PVS1:
    """PVS1 algorithm implementation."""

    def __init__(self, seqvar: SeqVar):
        self.seqvar = seqvar
        self.seqvar_transcripts = None
        self.gene_transcripts = None
        self.gene_hgnc_id = None
        self.id = self.id_generator()
        # self.vep_input = f"/tmp/vep.{self.id}.vcf"
        # self.vep_output = f"/tmp/vep.{self.id}.tab"

        # self.vep_variation = "na"
        # self.vep_symbol = "na"
        # self.vep_trans = "na"
        # self.vep_canonical = "na"
        # self.vep_pick = "na"
        # self.vep_consequence = "na"
        # self.hgvs_c = "na"
        # self.hgvs_p = "na"
        # self.hgvs_g = "na"
        # self.vep_exon = "na"
        # self.vep_intron = "na"

    @staticmethod
    def id_generator(size: int = 6, chars: str = string.ascii_uppercase + string.digits) -> str:
        """Generates a unique identifier for the VEP."""
        return "".join(random.choice(chars) for _ in range(size))

    async def _get_transcripts(self):
        """Get all transcripts for the given sequence variant."""
        try:
            mehari_client = MehariClient()
            response = await mehari_client.get_seqvar_transcripts(self.seqvar)
            if response["result"]:
                self.seqvar_transcripts = response["result"]
            else:
                self.seqvar_transcripts = None

            if self.seqvar_transcripts and len(self.seqvar_transcripts) > 0:
                self.gene_hgnc_id = self.seqvar_transcripts[0]["gene_id"]

                result = await mehari_client.get_gene_transcripts(
                    self.gene_hgnc_id, self.seqvar.genome_release
                )
                if result["transcripts"]:
                    self.gene_transcripts = result["transcripts"]
                else:
                    self.gene_transcripts = None

        except Exception as e:
            logger.error(
                f"Failed to get transcripts for variant {self.seqvar.user_representation}."
            )
            logger.error(e)

    def _is_nmd_target(self) -> bool:
        """
        Nonsense-mediated decay (NMD) classification. Return if the variant
        undergoes NMD prediction.
        See more at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6185798/#:~:text=Generally%2C%20NMD%20is%20not%20predicted%20to%20occur%20if%20the%20premature%20termination%20codon%20occurs%20in%20the%203%E2%80%99%20most%20exon%20or%20within%20the%203%E2%80%99%2Dmost%2050%20nucleotides%20of%20the%20penultimate%20exon

        Rule: If the variant is located in the last exon or in the last 50 nucleotides
        of the penultimate exon, it is NOT predicted to undergo NMD.
        :return: NMD prediction
        :rtype: bool
        """
        return False

    async def run(self):
        """Make the PVS1 prediction."""
        await self._get_transcripts()

        # self.vep
