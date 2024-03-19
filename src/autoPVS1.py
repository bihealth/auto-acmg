"""Implementations of the PVS1 algorithm."""

import logging
from typing import Any, Dict, List, Tuple

from src.api.mehari import MehariClient
from src.core.config import settings
from src.genome_builds import GenomeRelease
from src.seqvar import SeqVar, SeqVarResolver
from src.seqvar_pvs1 import SeqVarPVS1
from src.types import PVS1Prediction, SeqVarConsequence

# Setup logging
logging_level = logging.DEBUG if settings.DEBUG else logging.INFO
logging.basicConfig(level=logging_level)
logger = logging.getLogger(__name__)


#: TODO: Review the mapping
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
            self.seqvar = variant
            self.HGVSs: List[str] = []
            self.seqvar_ts_info = None
            self.gene_ts_info = None
            self.seqvar_transcript = None
            self.gene_transcript = None
            self.gene_hgnc_id = None
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
                    pvs1 = SeqVarPVS1(
                        self.seqvar, self.seqvar_transcript, self.gene_transcript, self.consequence
                    )
                    await pvs1.verify_PVS1()
                    self.prediction = pvs1.prediction
                    logger.info(
                        f"PVS1 prediction for {pvs1.seqvar.user_representation}: {pvs1.prediction}"
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
    def _get_consequence(seqvar_transcript: Dict[Any, Any]) -> SeqVarConsequence:
        """Get the consequence of the sequence variant."""
        if not seqvar_transcript:
            return SeqVarConsequence.NotSet
        else:
            for consequence in seqvar_transcript["consequences"]:
                if consequence in SeqvarConsequenceMapping:
                    return SeqvarConsequenceMapping[consequence]
            return SeqVarConsequence.NotSet

    @staticmethod
    def _choose_transcript(
        hgvss: List[str], seqvar_transcripts, gene_transcripts
    ) -> Tuple[Dict[Any, Any], Dict[Any, Any]] | Tuple[None, None]:
        """
        Choose the most suitable transcript for the PVS1 prediction.
        The first consideration is the MANE transcript, if available,
        and then the length of the exons.
        """
        transcripts = {}
        seqvar_transcript = None
        gene_transcript = None
        mane_transcripts = []
        exon_lengths = {}

        # Setup mapping from HGVS to pair of transcripts
        for hgvs in hgvss:
            transcripts[hgvs] = {"seqvar": None, "gene": None}
            for transcript in seqvar_transcripts:
                if transcript["feature_id"] == hgvs:
                    transcripts[hgvs]["seqvar"] = transcript
                    break
            for transcript in gene_transcripts:
                if transcript["id"] == hgvs:
                    transcripts[hgvs]["gene"] = transcript
                    break

        # Find MANE transcripts and calculate exon lengths
        for hgvs, transcript in transcripts.items():
            if transcript["seqvar"] and transcript["gene"]:
                if "ManeSelect" in transcript["seqvar"]["feature_tag"]:
                    mane_transcripts.append(hgvs)
                cds_sizes = [
                    exon["altEndI"] - exon["altStartI"]
                    for exon in transcript["gene"]["genomeAlignments"][0][
                        "exons"
                    ]  # TODO: Change 0 to proper Genome Build
                ]
                exon_lengths[hgvs] = sum(cds_sizes)

        # Choose the most suitable transcript
        if len(mane_transcripts) == 1:
            seqvar_transcript = transcripts[mane_transcripts[0]]["seqvar"]
            gene_transcript = transcripts[mane_transcripts[0]]["gene"]
        else:
            lookup_group = mane_transcripts if mane_transcripts else exon_lengths.items()
            max_length = max(lookup_group, key=lambda x: x[1])
            seqvar_transcript = transcripts[max_length[0]]["seqvar"]
            gene_transcript = transcripts[max_length[0]]["gene"]
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
            response = await mehari_client.get_seqvar_transcripts(self.seqvar)
            if response["result"]:
                self.seqvar_ts_info = response["result"]
            else:
                self.seqvar_ts_info = None

            if self.seqvar_ts_info and len(self.seqvar_ts_info) > 0:
                self.gene_hgnc_id = self.seqvar_ts_info[0]["gene_id"]

                for transcript in self.seqvar_ts_info:
                    self.HGVSs.append(transcript["feature_id"])

                result = await mehari_client.get_gene_transcripts(
                    self.gene_hgnc_id, self.seqvar.genome_release
                )

                if result["transcripts"]:
                    self.gene_ts_info = result["transcripts"]
                else:
                    self.gene_ts_info = None

            self.seqvar_transcript, self.gene_transcript = self._choose_transcript(
                self.HGVSs, self.seqvar_ts_info, self.gene_ts_info
            )

        except Exception as e:
            logger.error(
                f"Failed to get transcripts for variant {self.seqvar.user_representation}."
            )
            logger.error(e)
