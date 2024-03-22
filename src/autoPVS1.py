"""Implementations of the PVS1 algorithm."""

import logging
from typing import Dict, List, Tuple, Union

from src.api.mehari import MehariClient
from src.core.config import settings
from src.enums import PVS1Prediction, SeqVarConsequence
from src.genome_builds import GenomeRelease
from src.models.autopvs1 import TranscriptInfo
from src.models.mehari_gene import GeneTranscripts, TranscriptGene
from src.models.mehari_seqvar import TranscriptSeqvar, TranscriptsSeqVar
from src.seqvar import SeqVar, SeqVarResolver
from src.seqvar_pvs1 import SeqVarPVS1

# Setup logging
logging_level = logging.DEBUG if settings.DEBUG else logging.INFO
logging.basicConfig(level=logging_level)
logger = logging.getLogger(__name__)


#: Mapping of consequence from transcript info to SeqVarConsequence
SeqvarConsequenceMapping = {
    "intergenic_variant": SeqVarConsequence.NotSet,
    "intron_variant": SeqVarConsequence.NotSet,
    "upstream_gene_variant": SeqVarConsequence.NotSet,
    "downstream_gene_variant": SeqVarConsequence.NotSet,
    "5_prime_utr_variant": SeqVarConsequence.NotSet,
    "3_prime_utr_variant": SeqVarConsequence.NotSet,
    "splice_region_variant": SeqVarConsequence.SpliceSites,  # Can affect splicing
    "splice_donor_variant": SeqVarConsequence.SpliceSites,  # Canonical splice site
    "splice_acceptor_variant": SeqVarConsequence.SpliceSites,  # Canonical splice site
    "frameshift_variant": SeqVarConsequence.NonsenseFrameshift,  # Loss of function
    "transcript_ablation": SeqVarConsequence.NotSet,  # Severe effect, but not specifically classified here
    "transcript_amplification": SeqVarConsequence.NotSet,
    "inframe_insertion": SeqVarConsequence.NotSet,  # Not necessarily loss of function
    "inframe_deletion": SeqVarConsequence.NotSet,  # Not necessarily loss of function
    "synonymous_variant": SeqVarConsequence.NotSet,  # Usually benign, but exceptions exist
    "stop_retained_variant": SeqVarConsequence.NotSet,
    "missense_variant": SeqVarConsequence.NotSet,  # Not loss of function in a direct way
    "initiator_codon_variant": SeqVarConsequence.InitiationCodon,  # Affects start codon
    "stop_gained": SeqVarConsequence.NonsenseFrameshift,  # Nonsense variant
    "stop_lost": SeqVarConsequence.NotSet,  # Could be significant, but not classified here as Nonsense/Frameshift
    "mature_mirna_variant": SeqVarConsequence.NotSet,
    "non_coding_exon_variant": SeqVarConsequence.NotSet,  # Impact unclear
    "nc_transcript_variant": SeqVarConsequence.NotSet,
    "incomplete_terminal_codon_variant": SeqVarConsequence.NotSet,  # Rarely significant
    "nmd_transcript_variant": SeqVarConsequence.NotSet,  # Effect on NMD, not directly LOF
    "coding_sequence_variant": SeqVarConsequence.NotSet,  # Ambiguous
    "tfbs_ablation": SeqVarConsequence.NotSet,  # Regulatory, not LOF
    "tfbs_amplification": SeqVarConsequence.NotSet,
    "tf_binding_site_variant": SeqVarConsequence.NotSet,
    "regulatory_region_ablation": SeqVarConsequence.NotSet,  # Regulatory, not LOF
    "regulatory_region_variant": SeqVarConsequence.NotSet,
    "regulatory_region_amplification": SeqVarConsequence.NotSet,
}


class AutoPVS1:
    """AutoPVS1 algorithm for PVS1 criteria prediction."""

    def __init__(self, variant_name: str, genome_release: GenomeRelease = GenomeRelease.GRCh38):
        self.variant_name = variant_name
        self.genome_release = genome_release

    async def resolve_variant(self) -> SeqVar | None:
        """Resolve the variant."""
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
            self.seqvar: SeqVar = variant
            self.HGVSs: List[str] = []
            self.seqvar_ts_info: List[TranscriptSeqvar] | None = None
            self.gene_ts_info: List[TranscriptGene] | None = None
            self.seqvar_transcript: TranscriptSeqvar | None = None
            self.gene_transcript: TranscriptGene | None = None
            self.gene_hgnc_id: str = ""
            self.pvs1: SeqVarPVS1 | None = None
            self.consequence: SeqVarConsequence = SeqVarConsequence.NotSet
            self.prediction: PVS1Prediction = PVS1Prediction.NotPVS1

            logger.debug(f"Retrieving transcripts.")
            await self._get_transcripts_info()
            if not self.seqvar_transcript or not self.gene_transcript:
                logger.error("No transcripts found for the variant.")
                return
            else:
                try:
                    self.consequence = self._get_consequence(self.seqvar_transcript)
                    self.pvs1 = SeqVarPVS1(
                        self.seqvar, self.seqvar_transcript, self.gene_transcript, self.consequence
                    )
                    await self.pvs1.verify_PVS1()
                    self.prediction = self.pvs1.prediction
                    logger.info(
                        f"PVS1 prediction for {self.pvs1.seqvar.user_representation}: {self.pvs1.prediction}"
                    )
                except Exception as e:
                    logger.error(
                        f"Failed to predict PVS1 for variant {self.seqvar.user_representation}."
                    )
                    logger.error(e)
        elif isinstance(variant, str):
            # TODO: Add Structure variants PVS1 prediction
            pass
        else:
            logger.error(f"Failed to resolve variant {self.variant_name}.")
            return

    @staticmethod
    def _get_consequence(seqvar_transcript: TranscriptSeqvar | None) -> SeqVarConsequence:
        """Get the consequence of the sequence variant."""
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
        """
        Choose the most suitable transcript for the PVS1 prediction.
        The first consideration is the MANE transcript, if available,
        and then the length of the exons.
        """
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
            seqvar_transcript = transcripts_mapping[mane_transcripts[0]].seqvar
            gene_transcript = transcripts_mapping[mane_transcripts[0]].gene
        else:
            lookup_group = mane_transcripts if mane_transcripts else list(exon_lengths.keys())
            max_length_transcript = max(lookup_group, key=lambda x: exon_lengths[x])
            seqvar_transcript = transcripts_mapping[max_length_transcript].seqvar
            gene_transcript = transcripts_mapping[max_length_transcript].gene
        return seqvar_transcript, gene_transcript

    async def _get_transcripts_info(self):
        """Get all transcripts for the given sequence variant from Mehari."""
        if not self.seqvar:
            logger.error(
                "No sequence variant specified. Assure that the variant is resolved before fetching transcripts."
            )
            return
        try:
            mehari_client = MehariClient()
            response_seqvar = await mehari_client.get_seqvar_transcripts(self.seqvar)
            if not response_seqvar:
                self.seqvar_ts_info = None
            else:
                self.seqvar_ts_info = response_seqvar.result

            if self.seqvar_ts_info and len(self.seqvar_ts_info) > 0:
                self.gene_hgnc_id = self.seqvar_ts_info[0].gene_id

                for transcript in self.seqvar_ts_info:
                    self.HGVSs.append(transcript.feature_id)

                response_gene = await mehari_client.get_gene_transcripts(
                    self.gene_hgnc_id, self.seqvar.genome_release
                )

                if not response_gene:
                    self.gene_ts_info = None
                else:
                    self.gene_ts_info = response_gene.transcripts

            if self.seqvar_ts_info and self.gene_ts_info:
                self.seqvar_transcript, self.gene_transcript = self._choose_transcript(
                    self.HGVSs, self.seqvar_ts_info, self.gene_ts_info
                )
            else:
                self.seqvar_transcript = None
                self.gene_transcript = None

        except Exception as e:
            logger.error(
                f"Failed to get transcripts for variant {self.seqvar.user_representation}."
            )
            logger.error(e)
