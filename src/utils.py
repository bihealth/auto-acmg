"""Utility functions for the AutoACMG and AutoPVS1."""

import itertools
import re
from typing import Dict, List, Optional, Tuple

from biocommons.seqrepo import SeqRepo  # type: ignore
from loguru import logger

from lib.maxentpy import maxent
from lib.maxentpy.maxent import load_matrix3, load_matrix5
from src.api.annonars import AnnonarsClient
from src.api.mehari import MehariClient
from src.core.config import Config
from src.defs.annonars_variant import VariantResult
from src.defs.auto_acmg import SpliceType
from src.defs.auto_pvs1 import SeqVarConsequence, SeqvarConsequenceMapping, TranscriptInfo
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException
from src.defs.genome_builds import CHROM_REFSEQ_37, CHROM_REFSEQ_38, GenomeRelease
from src.defs.mehari import TranscriptGene, TranscriptSeqvar
from src.defs.seqvar import SeqVar


class SplicingPrediction:
    def __init__(self, seqvar: SeqVar, *, config: Optional[Config] = None):
        self.donor_threshold = 3
        self.acceptor_threshold = 3
        self.percent_threshold = 0.7
        self.seqvar = seqvar
        self.config: Config = config or Config()
        self.annonars_client = AnnonarsClient(api_base_url=self.config.api_base_url_annonars)
        self.sr = SeqRepo(self.config.seqrepo_data_dir)

        self.maxentscore_ref = -1.00
        self.maxentscore_alt = -1.00
        self.maxent_foldchange = 1.00
        self.matrix5 = load_matrix5()
        self.matrix3 = load_matrix3()

    def get_sequence(self, start: int, end: int) -> str:
        """Retrieve the sequence for the specified range."""
        # Get the RefSeq chromosome name from the SeqVar chromosome
        if self.seqvar.genome_release == GenomeRelease.GRCh37:
            chrom = CHROM_REFSEQ_37[self.seqvar.chrom]
        elif self.seqvar.genome_release == GenomeRelease.GRCh38:
            chrom = CHROM_REFSEQ_38[self.seqvar.chrom]
        else:
            logger.error("Invalid genome release: {}", self.seqvar.genome_release)
            raise AlgorithmError("Invalid genome release.")
        try:
            seq = self.sr[chrom][start:end]
            return seq
        except Exception as e:
            logger.error("Failed to get sequence for {}:{}-{}. Error: {}", chrom, start, end, e)
            raise AlgorithmError("Failed to get sequence for the specified range.") from e

    def _calculate_maxentscore(
        self, refseq: str, altseq: str, splice_type: str
    ) -> Tuple[float, float, float]:
        """
        --- Calculate the maxentscan socre ---
        When a mutation occurs, if the WT score is above the threshold and
        the score variation (between WT and Mutant) is under -10% for HSF (-30% for MaxEnt)
        we consider that the mutation breaks the splice site.
        In the other case, if the WT score is under the threshold and
        the score variation is above +10% for HSF (+30% for MaxEnt) we consider that
        the mutation creates a new splice site.
        """
        maxentscore_ref = maxentscore_alt = -1.00
        if splice_type == "donor":
            if len(refseq) == 9:
                maxentscore_ref = maxent.score5(refseq, matrix=self.matrix5)
            if len(altseq) == 9:
                maxentscore_alt = maxent.score5(altseq, matrix=self.matrix5)
        elif splice_type == "acceptor":
            if len(refseq) == 23:
                maxentscore_ref = maxent.score3(refseq, matrix=self.matrix3)
            if len(altseq) == 23:
                maxentscore_alt = maxent.score3(altseq, matrix=self.matrix3)

        maxent_foldchange = maxentscore_alt / maxentscore_ref

        return round(maxentscore_ref, 2), round(maxentscore_alt, 2), round(maxent_foldchange, 2)

    def get_cryptic_ss(self, refseq: str, splice_type: SpliceType) -> List[Tuple[int, str, float]]:
        """
        Get cryptic splice sites around the variant position.

        Args:
            refseq: The reference sequence in the specified region.
            splice_type: The type of splice site, either "donor" or "acceptor".

        Returns:
            List[Tuple[int, str, float]]: List of cryptic splice sites with position, context, and score.
        """
        if splice_type == SpliceType.Unknown:
            logger.warning("Unknown splice type. Cannot predict cryptic splice sites.")
            return []

        refscore = self.maxentscore_ref
        search_flank = 20  # Flank size for checking positions around the variant
        cryptic_sites = []

        for offset in range(-search_flank, search_flank):
            pos = self.seqvar.pos + offset
            if splice_type == SpliceType.Donor:
                splice_context = self.get_sequence(pos, pos + 9)
                alt_index = self.seqvar.pos - pos
                if 0 < alt_index < 9:
                    splice_context = (
                        splice_context[:alt_index]
                        + self.seqvar.insert
                        + splice_context[
                            alt_index + len(self.seqvar.insert) : 9 - len(self.seqvar.insert)
                        ]
                    )
                if len(splice_context) == 9:
                    maxentscore = maxent.score5(splice_context, matrix=self.matrix5)
                else:
                    maxentscore = 0
                if (
                    splice_context[3:5] in ["GT", refseq[3:5]]
                    and maxentscore > 1
                    and (
                        maxentscore >= self.donor_threshold
                        or maxentscore / refscore >= self.percent_threshold
                    )
                ):
                    cryptic_sites.append((pos, splice_context, maxentscore))

            elif splice_type == SpliceType.Acceptor:
                splice_context = self.get_sequence(pos, pos + 23)
                alt_index = self.seqvar.pos - pos - 1
                if 0 < alt_index < 23:
                    splice_context = (
                        splice_context[:alt_index]
                        + self.seqvar.insert
                        + splice_context[
                            alt_index + len(self.seqvar.insert) : 23 - len(self.seqvar.insert)
                        ]
                    )
                if len(splice_context) == 23:
                    maxentscore = maxent.score3(splice_context, matrix=self.matrix3)
                else:
                    maxentscore = 0
                if (
                    splice_context[18:20] in ["AG", refseq[18:20]]
                    and maxentscore > 1
                    and (
                        maxentscore >= self.acceptor_threshold
                        or maxentscore / refscore >= self.percent_threshold
                    )
                ):
                    cryptic_sites.append((pos, splice_context, maxentscore))

        return cryptic_sites


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
