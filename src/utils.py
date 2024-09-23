"""Utility functions for the AutoACMG and AutoPVS1."""

from typing import Dict, List, Tuple, Union

from biocommons.seqrepo import SeqRepo
from loguru import logger

from lib.maxentpy import maxent
from lib.maxentpy.maxent import load_matrix3, load_matrix5
from src.api.reev.annonars import AnnonarsClient
from src.api.reev.mehari import MehariClient
from src.core.config import settings
from src.defs.auto_acmg import GenomicStrand, SpliceType, TranscriptInfo
from src.defs.auto_pvs1 import SeqvarConsequenceMapping, SeqVarPVS1Consequence
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException
from src.defs.genome_builds import CHROM_REFSEQ_37, CHROM_REFSEQ_38, GenomeRelease
from src.defs.mehari import Exon, TranscriptGene, TranscriptSeqvar, TranscriptStrucvar
from src.defs.seqvar import SeqVar
from src.defs.strucvar import StrucVar


class AutoACMGHelper:
    """Helper class for the AutoACMG algorithm."""

    def __init__(self):
        #: Annonars client for the API.
        self.annonars_client = AnnonarsClient(api_base_url=settings.AUTO_ACMG_API_ANNONARS_URL)


class SplicingPrediction:
    """Splicing prediction for a sequence variant."""

    def __init__(
        self,
        seqvar: SeqVar,
        *,
        strand: GenomicStrand,
        consequences: List[str] = [],
        exons: List[Exon],
    ):
        self.donor_threshold = 3
        self.acceptor_threshold = 3
        self.percent_threshold = 0.7
        self.seqvar = seqvar
        self.strand = strand
        self.exons = exons
        self.splice_type = self.determine_splice_type(consequences)
        self.annonars_client = AnnonarsClient(api_base_url=settings.AUTO_ACMG_API_ANNONARS_URL)
        self.sr = SeqRepo(settings.AUTO_ACMG_SEQREPO_DATA_DIR)

        self.maxentscore_ref = -1.00
        self.maxentscore_alt = -1.00
        self.maxent_foldchange = 1.00
        self.matrix5 = load_matrix5()
        self.matrix3 = load_matrix3()
        self._initialize_maxentscore()

    def _initialize_maxentscore(self):
        """
        Initialize the MaxEntScan scores for the sequence variant.

        Find the reference sequence, generate the alternative sequence, and calculate the
        MaxEntScan score for the reference and alternative sequences.
        """
        refseq_start, _, refseq = self._find_reference_sequence()
        altseq = self._generate_alt_sequence(refseq, refseq_start)
        self.maxentscore_ref, self.maxentscore_alt, self.maxent_foldchange = (
            self._calculate_maxentscore(refseq, altseq, self.splice_type)
        )

    def _find_reference_sequence(self) -> Tuple[int, int, str]:
        """Find the reference sequence based on the exon positions and splice type."""
        if self.strand == GenomicStrand.Plus:
            return self._find_refseq_plus_strand()
        elif self.strand == GenomicStrand.Minus:
            return self._find_refseq_minus_strand()
        else:
            logger.error("Invalid genomic strand: {}", self.strand.name)
            raise AlgorithmError("Invalid genomic strand.")

    def _find_refseq_plus_strand(self) -> Tuple[int, int, str]:
        refseq = ""
        refseq_start, refseq_end = 0, 0
        for i, exon in enumerate(self.exons):
            if self.splice_type == SpliceType.Donor:
                if (
                    exon.altEndI <= self.seqvar.pos
                    and self.exons[i + 1].altStartI >= self.seqvar.pos
                ):
                    refseq_start, refseq_end = exon.altEndI - 3, exon.altEndI + 6
                    refseq = self.get_sequence(refseq_start, refseq_end)
                    break
            elif self.splice_type == SpliceType.Acceptor:
                if exon.altStartI >= self.seqvar.pos:
                    refseq_start, refseq_end = exon.altStartI - 3, exon.altStartI + 20
                    refseq = self.get_sequence(refseq_start, refseq_end)
                    break
        return refseq_start, refseq_end, refseq

    def _find_refseq_minus_strand(self) -> Tuple[int, int, str]:
        refseq = ""
        refseq_start, refseq_end = 0, 0
        for i, exon in enumerate(self.exons):
            if self.splice_type == SpliceType.Donor:
                if exon.altStartI >= self.seqvar.pos:
                    refseq_start, refseq_end = exon.altStartI - 6, exon.altStartI + 3
                    refseq = self.reverse_complement(self.get_sequence(refseq_start, refseq_end))
                    break
            elif self.splice_type == SpliceType.Acceptor:
                if (
                    exon.altEndI <= self.seqvar.pos
                    and self.exons[i + 1].altStartI >= self.seqvar.pos
                ):
                    refseq_start, refseq_end = exon.altEndI - 20, exon.altEndI + 3
                    refseq = self.reverse_complement(self.get_sequence(refseq_start, refseq_end))
                    break
        return refseq_start, refseq_end, refseq

    def _generate_alt_sequence(self, refseq: str, refseq_start: int) -> str:
        """Generate the alternative sequence with the variant nucleotide(s) inserted."""
        altseq = ""
        altseq_index = self.seqvar.pos - refseq_start - 1
        # Splice length is 9 for Donor and 23 for Acceptor
        splice_length = 9 if self.splice_type == SpliceType.Donor else 23
        for i in range(len(refseq)):
            if len(altseq) == splice_length:
                break
            if i < altseq_index:
                altseq += refseq[i]
            elif i == altseq_index:
                # Insert the variant nucleotide(s) into the reference sequence
                altseq += self.seqvar.insert
                if len(altseq) == splice_length:
                    break
                elif len(altseq) > splice_length:
                    altseq = altseq[:splice_length]
                    break
                else:
                    # Add the remaining nucleotides from the reference sequence
                    for j in range(i + 1, len(refseq)):
                        altseq += refseq[j]
                        if len(altseq) == splice_length:
                            break
        return altseq

    def _calculate_maxentscore(
        self, refseq: str, altseq: str, splice_type: SpliceType
    ) -> Tuple[float, float, float]:
        """
        Calculate the MaxEntScan score for the reference and alternative sequences.

        When a mutation occurs, if the WT score is above the threshold and
        the score variation (between WT and Mutant) is under -10% for HSF (-30% for MaxEnt)
        we consider that the mutation breaks the splice site.
        In the other case, if the WT score is under the threshold and
        the score variation is above +10% for HSF (+30% for MaxEnt) we consider that
        the mutation creates a new splice site.

        Args:
            refseq: The reference sequence.
            altseq: The alternative sequence.
            splice_type: The type of splice site, either "donor" or "acceptor".

        Returns:
            Tuple[float, float, float]: The MaxEntScan score for the reference sequence,
            the alternative sequence, and the fold change.
        """
        if splice_type == SpliceType.Unknown:
            logger.warning("Unknown splice type. Cannot calculate MaxEntScan score.")
            return -1.00, -1.00, 1.00

        maxentscore_ref, maxentscore_alt = -1.00, -1.00
        if splice_type == SpliceType.Donor:
            if len(refseq) == 9:
                maxentscore_ref = maxent.score5(refseq, matrix=self.matrix5)
            if len(altseq) == 9:
                maxentscore_alt = maxent.score5(altseq, matrix=self.matrix5)
        elif splice_type == SpliceType.Acceptor:
            if len(refseq) == 23:
                maxentscore_ref = maxent.score3(refseq, matrix=self.matrix3)
            if len(altseq) == 23:
                maxentscore_alt = maxent.score3(altseq, matrix=self.matrix3)

        maxent_foldchange = maxentscore_alt / maxentscore_ref
        return (
            round(maxentscore_ref, 2),
            round(maxentscore_alt, 2),
            round(maxent_foldchange, 2),
        )

    @staticmethod
    def determine_splice_type(consequences: List[str]) -> SpliceType:
        """Determine the splice type based on the consequence."""
        splice_type = SpliceType.Unknown
        for consequence in consequences:
            match consequence:
                case "splice_acceptor_variant":
                    splice_type = SpliceType.Acceptor
                    break
                case "splice_donor_variant":
                    splice_type = SpliceType.Donor
                    break
                case _:
                    continue
        return splice_type

    @staticmethod
    def reverse_complement(seq: str) -> str:
        """Retrun a reverse complementary seq"""
        nt_complement = {
            "A": "T",
            "C": "G",
            "T": "A",
            "G": "C",
            "a": "t",
            "c": "g",
            "t": "a",
            "g": "c",
        }
        reverse_seq = list(reversed(seq))
        rev_comp_seq = [nt_complement[k] for k in reverse_seq]
        return "".join(rev_comp_seq)

    def get_sequence(self, start: int, end: int) -> str:
        """
        Retrieve the sequence for the specified range.
        The reference sequence is chosen based on the chromosome and genome release from the
        sequence variant, which was provided by the initialisation of the class.

        Args:
            start: The start position of the sequence.
            end: The end position of the sequence.

        Returns:
            str: The sequence for the specified range.

        Raises:
            AlgorithmError: If the sequence cannot be retrieved.
        """
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
            if self.strand == GenomicStrand.Minus:
                seq = self.reverse_complement(seq)
            return seq
        except Exception as e:
            logger.error("Failed to get sequence for {}:{}-{}. Error: {}", chrom, start, end, e)
            raise AlgorithmError("Failed to get sequence for the specified range.") from e

    def get_cryptic_ss(self, refseq: str, splice_type: SpliceType) -> List[Tuple[int, str, float]]:
        """
        Get cryptic splice sites around the variant position.
        """
        if splice_type == SpliceType.Unknown:
            logger.warning("Unknown splice type. Cannot predict cryptic splice sites.")
            return []

        refscore = self.maxentscore_ref
        search_flank = 20  # Flank size for checking positions around the variant
        cryptic_sites = []

        for offset in range(-search_flank, search_flank + 1):
            pos = self.seqvar.pos + offset
            if splice_type == SpliceType.Donor:
                cryptic_sites.extend(self._find_cryptic_donor_sites(pos, refseq, refscore))
            elif splice_type == SpliceType.Acceptor:
                cryptic_sites.extend(self._find_cryptic_acceptor_sites(pos, refseq, refscore))

        return cryptic_sites

    def _find_cryptic_donor_sites(
        self, pos: int, refseq: str, refscore: float
    ) -> List[Tuple[int, str, float]]:
        """Find cryptic donor sites."""
        splice_context = self.get_sequence(pos, pos + 9)
        if self.strand == GenomicStrand.Minus:
            splice_context = self.reverse_complement(splice_context)
        maxentscore = maxent.score5(splice_context, matrix=self.matrix5)

        if (
            splice_context[3:5] in ["GT", refseq[3:5]]
            and maxentscore > 1
            and (
                maxentscore >= self.donor_threshold
                or maxentscore / refscore >= self.percent_threshold
            )
        ):
            return [(pos, splice_context, maxentscore)]
        return []

    def _find_cryptic_acceptor_sites(
        self, pos: int, refseq: str, refscore: float
    ) -> List[Tuple[int, str, float]]:
        """Find cryptic acceptor sites."""
        splice_context = self.get_sequence(pos, pos + 23)
        if self.strand == GenomicStrand.Minus:
            splice_context = self.reverse_complement(splice_context)
        maxentscore = maxent.score3(splice_context, matrix=self.matrix3)

        if (
            splice_context[18:20] in ["AG", refseq[18:20]]
            and maxentscore > 1
            and (
                maxentscore >= self.acceptor_threshold
                or maxentscore / refscore >= self.percent_threshold
            )
        ):
            return [(pos, splice_context, maxentscore)]
        return []


class SeqVarTranscriptsHelper:
    """Transcript information for a sequence variant."""

    def __init__(self, seqvar: SeqVar):
        self.seqvar: SeqVar = seqvar

        # Attributes to be set
        self.HGVSs: List[str] = []
        self.HGNC_id: str = ""
        self.seqvar_ts_info: List[TranscriptSeqvar] = []
        self.seqvar_transcript: TranscriptSeqvar | None = None
        self.gene_ts_info: List[TranscriptGene] = []
        self.gene_transcript: TranscriptGene | None = None
        self.consequence: SeqVarPVS1Consequence = SeqVarPVS1Consequence.NotSet

    def get_ts_info(
        self,
    ) -> Tuple[
        TranscriptSeqvar | None,
        TranscriptGene | None,
        List[TranscriptSeqvar],
        List[TranscriptGene],
        SeqVarPVS1Consequence,
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
            mehari_client = MehariClient(api_base_url=settings.AUTO_ACMG_API_MEHARI_URL)
            response_seqvar = mehari_client.get_seqvar_transcripts(self.seqvar)
            if not response_seqvar:
                self.seqvar_ts_info = []
            else:
                self.seqvar_ts_info = response_seqvar.result

            if not self.seqvar_ts_info or len(self.seqvar_ts_info) == 0:
                self.seqvar_transcript = None
                self.gene_transcript = None
                self.consequence = SeqVarPVS1Consequence.NotSet
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
                self.consequence = SeqVarPVS1Consequence.NotSet

        except AutoAcmgBaseException as e:
            raise AlgorithmError("Failed to get transcripts for the sequence variant.") from e

    @staticmethod
    def _get_consequence(
        seqvar_transcript: TranscriptSeqvar | None,
    ) -> SeqVarPVS1Consequence:
        """Get the consequence of the sequence variant.

        Args:
            seqvar_transcript: The sequence variant transcript.

        Returns:
            SeqVarConsequence: The consequence of the sequence variant.
        """
        if not seqvar_transcript:
            return SeqVarPVS1Consequence.NotSet
        else:
            for consequence in seqvar_transcript.consequences:
                if consequence in SeqvarConsequenceMapping:
                    return SeqvarConsequenceMapping[consequence]
            return SeqVarPVS1Consequence.NotSet

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


class StrucVarTranscriptsHelper:
    """Transcript information for a structural variant."""

    def __init__(self, strucvar: StrucVar):
        self.strucvar: StrucVar = strucvar

        # Attributes to be set
        self.strucvar_ts_info: List[TranscriptStrucvar] = []
        self.strucvar_transcript: TranscriptStrucvar | None = None
        self.gene_ts_info: List[TranscriptGene] = []
        self.gene_transcript: TranscriptGene | None = None

    def get_ts_info(
        self,
    ) -> Tuple[
        TranscriptStrucvar | None,
        TranscriptGene | None,
        List[TranscriptStrucvar],
        List[TranscriptGene],
    ]:
        """Return the transcript information.

        Returns:
            Tuple[TranscriptStrucvar | None, TranscriptGene | None, List[TranscriptStrucvar],
            List[TranscriptGene]]:
            The structural variant transcript,
            gene transcript, and the consequence of the sequence variant.
        """
        return (
            self.strucvar_transcript,
            self.gene_transcript,
            self.strucvar_ts_info,
            self.gene_ts_info,
        )

    def initialize(self):
        """Get all transcripts for the given structural variant from Mehari."""
        try:
            # Get transcripts from Mehari
            mehari_client = MehariClient(api_base_url=settings.AUTO_ACMG_API_MEHARI_URL)
            response_strucvar = mehari_client.get_strucvar_transcripts(self.strucvar)
            if not response_strucvar:
                self.strucvar_ts_info = []
            else:
                self.strucvar_ts_info = response_strucvar.result

            if not self.strucvar_ts_info or len(self.strucvar_ts_info) == 0:
                self.strucvar_transcript = None
                self.gene_transcript = None
                logger.warning("No transcripts found for the structural variant.")
                return
            else:
                self.strucvar_transcript = self.strucvar_ts_info[0]

            # Get the first HGNC ID
            self.HGNC_id = self.strucvar_transcript.hgnc_id
            # Get gene transcripts from Mehari
            response_gene = mehari_client.get_gene_transcripts(
                self.HGNC_id, self.strucvar.genome_release
            )
            if not response_gene:
                self.gene_ts_info = []
            else:
                self.gene_ts_info = response_gene.transcripts

            # Choose the most suitable gene transcript
            self.gene_transcript = self._choose_transcript(self.gene_ts_info)

        except AutoAcmgBaseException as e:
            raise AlgorithmError("Failed to get transcripts for the structural variant.") from e

    @staticmethod
    def _choose_transcript(
        gene_transcripts: List[TranscriptGene],
    ) -> Union[TranscriptGene, None]:
        """Choose the most suitable transcript for the PVS1 prediction.

        Note:
            The first consideration is the MANE transcript, if available,
            and then the length of the exons.

        Args:
            gene_transcripts: The gene transcripts.

        Returns:
            TranscriptGene | None: The most suitable transcript (with ManeSelect tag) for the PVS1
            prediction.
        """
        if not gene_transcripts:
            return None
        mane_transcripts: List[str] = []
        exon_lengths: Dict[str, int] = {}
        transcripts_mapping: Dict[str, TranscriptGene] = {}
        for transcript in gene_transcripts:
            transcripts_mapping[transcript.id] = transcript
            if transcript.tags and "TRANSCRIPT_TAG_MANE_SELECT" in transcript.tags:
                mane_transcripts.append(transcript.id)
            cds_sizes = [
                exon.altEndI - exon.altStartI for exon in transcript.genomeAlignments[0].exons
            ]
            exon_lengths[transcript.id] = sum(cds_sizes)

        if len(mane_transcripts) == 1:
            logger.debug("The MANE transcript found: {}", mane_transcripts[0])
            return transcripts_mapping[mane_transcripts[0]]
        else:
            logger.debug("Choosing the longest transcript.")
            lookup_group = mane_transcripts if mane_transcripts else list(exon_lengths.keys())
            if not lookup_group or not exon_lengths:
                logger.warning("No overlapping transcripts found.")
                return None
            max_length_transcript = max(lookup_group, key=lambda x: exon_lengths[x])
            return transcripts_mapping[max_length_transcript]
