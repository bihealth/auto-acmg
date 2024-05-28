"""Utility functions for the AutoACMG and AutoPVS1."""

from typing import Dict, List, Optional, Tuple

from loguru import logger

from src.api.mehari import MehariClient
from src.core.config import Config
from src.defs.auto_pvs1 import SeqVarConsequence, SeqvarConsequenceMapping, TranscriptInfo
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException
from src.defs.mehari import TranscriptGene, TranscriptSeqvar
from src.defs.seqvar import SeqVar


class SeqVarTranscriptsHelper:
    """Transcript information for a sequence variant."""

    def __init__(self, seqvar: SeqVar, *, config: Optional[Config] = None):
        self.config: Config = config or Config()
        self.seqvar: SeqVar = seqvar

        # Attributes to be set
        self.HGVSs: List[str] = []
        self.HGNC_id: str = ""
        self.seqvar_ts_info: List[TranscriptSeqvar] = []
        self.seqvar_transcript: TranscriptSeqvar | None = None
        self.gene_ts_info: List[TranscriptGene] = []
        self.gene_transcript: TranscriptGene | None = None
        self.consequence: SeqVarConsequence = SeqVarConsequence.NotSet

    def get_ts_info(
        self,
    ) -> Tuple[
        TranscriptSeqvar | None,
        TranscriptGene | None,
        List[TranscriptSeqvar],
        List[TranscriptGene],
        SeqVarConsequence,
    ]:
        """Return the transcript information.

        Returns:
            Tuple[TranscriptSeqvar | None, TranscriptGene | None, List[TranscriptSeqvar],
            List[TranscriptGene], SeqVarConsequence]:
            The sequence variant transcript,
            gene transcript, and the consequence of the sequence variant.
        """
        return (
            self.seqvar_transcript,
            self.gene_transcript,
            self.seqvar_ts_info,
            self.gene_ts_info,
            self.consequence,
        )

    def initialize(self):
        """Get all transcripts for the given sequence variant from Mehari."""
        try:
            # Get transcripts from Mehari
            mehari_client = MehariClient(api_base_url=self.config.api_base_url_mehari)
            response_seqvar = mehari_client.get_seqvar_transcripts(self.seqvar)
            if not response_seqvar:
                self.seqvar_ts_info = []
            else:
                self.seqvar_ts_info = response_seqvar.result

            if not self.seqvar_ts_info or len(self.seqvar_ts_info) == 0:
                self.seqvar_transcript = None
                self.gene_transcript = None
                self.consequence = SeqVarConsequence.NotSet
                logger.warning("No transcripts found for the sequence variant.")
                return

            # Get HGNC ID and HGVSs
            self.HGNC_id = self.seqvar_ts_info[0].gene_id
            for transcript in self.seqvar_ts_info:
                self.HGVSs.append(transcript.feature_id)

            # Get gene transcripts from Mehari
            response_gene = mehari_client.get_gene_transcripts(
                self.HGNC_id, self.seqvar.genome_release
            )
            if not response_gene:
                self.gene_ts_info = []
            else:
                self.gene_ts_info = response_gene.transcripts

            # Choose the most suitable transcript for the PVS1 prediction
            if self.seqvar_ts_info and self.gene_ts_info:
                self.seqvar_transcript, self.gene_transcript = self._choose_transcript(
                    self.HGVSs, self.seqvar_ts_info, self.gene_ts_info
                )
                self.consequence = self._get_consequence(self.seqvar_transcript)
            else:
                self.seqvar_transcript = None
                self.gene_transcript = None
                self.consequence = SeqVarConsequence.NotSet

        except AutoAcmgBaseException as e:
            logger.error("Failed to get transcripts for the sequence variant. Error: {}", e)
            raise AlgorithmError("Failed to get transcripts for the sequence variant.") from e

    @staticmethod
    def _get_consequence(seqvar_transcript: TranscriptSeqvar | None) -> SeqVarConsequence:
        """Get the consequence of the sequence variant.

        Args:
            seqvar_transcript: The sequence variant transcript.

        Returns:
            SeqVarConsequence: The consequence of the sequence variant.
        """
        if not seqvar_transcript:
            return SeqVarConsequence.NotSet
        else:
            for consequence in seqvar_transcript.consequences:
                if consequence in SeqvarConsequenceMapping:
                    return SeqvarConsequenceMapping[consequence]
            return SeqVarConsequence.NotSet

    @staticmethod
    def _choose_transcript(
        hgvss: List[str],
        seqvar_transcripts: List[TranscriptSeqvar],
        gene_transcripts: List[TranscriptGene],
    ) -> Tuple[TranscriptSeqvar | None, TranscriptGene | None]:
        """Choose the most suitable transcript for the PVS1 prediction.

        Note:
            The first consideration is the MANE transcript, if available,
            and then the length of the exons.
        """
        logger.debug("Choosing the most suitable transcript for the PVS1 prediction.")
        transcripts_mapping: Dict[str, TranscriptInfo] = {}
        seqvar_transcript = None
        gene_transcript = None
        mane_transcripts: List[str] = []
        exon_lengths: Dict[str, int] = {}

        # Setup mapping from HGVS to pair of transcripts
        for hgvs in hgvss:
            transcripts_mapping[hgvs] = TranscriptInfo(seqvar=None, gene=None)
            for seqvar_ts in seqvar_transcripts:
                if seqvar_ts.feature_id == hgvs:
                    transcripts_mapping[hgvs].seqvar = seqvar_ts
                    break
            for gene_ts in gene_transcripts:
                if gene_ts.id == hgvs:
                    transcripts_mapping[hgvs].gene = gene_ts
                    break

        # Find MANE transcripts and calculate exon lengths
        for hgvs, transcript in transcripts_mapping.items():
            if transcript.seqvar and transcript.gene:
                if "ManeSelect" in transcript.seqvar.feature_tag:
                    mane_transcripts.append(hgvs)
                cds_sizes = [
                    exon.altEndI - exon.altStartI
                    for exon in transcript.gene.genomeAlignments[0].exons
                ]
                exon_lengths[hgvs] = sum(cds_sizes)

        # Choose the most suitable transcript
        if len(mane_transcripts) == 1:
            logger.debug("The MANE transcript found: {}", mane_transcripts[0])
            seqvar_transcript = transcripts_mapping[mane_transcripts[0]].seqvar
            gene_transcript = transcripts_mapping[mane_transcripts[0]].gene
        else:
            logger.debug("Choosing the longest transcript.")
            lookup_group = mane_transcripts if mane_transcripts else list(exon_lengths.keys())
            # If no overlapping transcripts, return None
            if not lookup_group or not exon_lengths:
                logger.warning("No overlapping transcripts found.")
                return None, None
            max_length_transcript = max(lookup_group, key=lambda x: exon_lengths[x])
            seqvar_transcript = transcripts_mapping[max_length_transcript].seqvar
            gene_transcript = transcripts_mapping[max_length_transcript].gene
        return seqvar_transcript, gene_transcript
