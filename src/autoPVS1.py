"""Implementations of the PVS1 algorithm."""

import logging
from typing import Dict, List

from src.api.mehari import MehariClient
from src.core.config import settings
from src.genome_builds import GenomeRelease
from src.pvs1_types import PVS1Prediction
from src.seqvar import SeqVar, SeqVarResolver
from src.seqvar_pvs1 import SeqVarPVS1

# Setup logging
logging_level = logging.DEBUG if settings.DEBUG else logging.INFO
logging.basicConfig(level=logging_level)
logger = logging.getLogger(__name__)


class AutoPVS1:
    """AutoPVS1 algorithm for PVS1 criteria prediction."""

    def __init__(self, variant_name: str, genome_release: GenomeRelease = GenomeRelease.GRCh38):
        self.variant_name = variant_name
        self.genome_release = genome_release

    async def resolve_variant(self) -> SeqVar | None:
        """
        Resolve the variant.

        :return: Sequence variant
        :rtype: SeqVar | None
        """
        # TODO: Add resolve for Structure variants
        try:
            # Try to resolve as Sequence variant
            seqvar_resolver = SeqVarResolver()
            seqvar: SeqVar = await seqvar_resolver.resolve_seqvar(
                self.variant_name, self.genome_release
            )
            logger.debug(f"Resolved variant: {seqvar}.")
            return seqvar
        except Exception as e:
            logger.error(e)
            return None

    async def predict(self):
        """Run the AutoPVS1 algorithm."""
        logger.info(f"Running AutoPVS1 for variant {self.variant_name}.")
        variant = await self.resolve_variant()
        if isinstance(variant, SeqVar):
            self.seqvar = variant
            self.HGVSs: List[str] = []
            self.seqvar_transcripts = None
            self.gene_transcripts = None
            self.gene_hgnc_id = None
            self.predictions: Dict[str, PVS1Prediction] = {}

            logger.debug(f"Retrieving transcripts.")
            await self._get_transcripts_seqvar()
            if self.HGVSs:
                for hgvs in self.HGVSs:
                    logger.info(f"Analyzing transcript {hgvs}.")
                    # Choose transcripts
                    seqvar_t = None
                    gene_t = None
                    for transcript_info in self.seqvar_transcripts:
                        if transcript_info["feature_id"] == hgvs:
                            seqvar_t = transcript_info
                            break
                    for transcript_info in self.gene_transcripts:
                        if transcript_info["id"] == hgvs:
                            gene_t = transcript_info
                            break

                    if seqvar_t and gene_t:
                        # TODO: Add support for SpliceSites and InitiationCodon consequences
                        pvs1 = SeqVarPVS1(self.seqvar, seqvar_t, gene_t)
                        pvs1.verify_PVS1()
                        self.predictions[hgvs] = pvs1.prediction
                        logger.debug(f"PVS1 prediction for {hgvs}: {pvs1.prediction.name}")

            else:
                logger.debug("No transcripts found for the variant.")
        elif isinstance(variant, str):
            # TODO: Add Structure variants PVS1 prediction
            pass
        else:
            logger.error(f"Failed to resolve variant {self.variant_name}.")
            return

    async def _get_transcripts_seqvar(self):
        """Get all transcripts for the given sequence variant."""
        if not self.seqvar:
            logger.error(
                "No sequence variant specified. Assure that the variant is resolved before fetching transcripts."
            )
            return
        try:
            mehari_client = MehariClient()
            response = await mehari_client.get_seqvar_transcripts(self.seqvar)
            if response["result"]:
                self.seqvar_transcripts = response["result"]
            else:
                self.seqvar_transcripts = None

            if self.seqvar_transcripts and len(self.seqvar_transcripts) > 0:
                self.gene_hgnc_id = self.seqvar_transcripts[0]["gene_id"]

                for transcript in self.seqvar_transcripts:
                    self.HGVSs.append(transcript["feature_id"])

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
