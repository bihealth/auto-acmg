"""PVS1 criteria for Sequence Variants (SeqVar)."""

import re
from typing import Dict, List, Tuple

import typer

from src.api.annonars import AnnonarsClient
from src.api.mehari import MehariClient
from src.core.exceptions import InvalidAPIResposeError
from src.seqvar import SeqVar
from src.types.autopvs1 import TranscriptInfo
from src.types.enums import PVS1Prediction, PVS1PredictionSeqVarPath, SeqVarConsequence
from src.types.mehari import CdsPos, Exon, TranscriptGene, TranscriptSeqvar

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


class SeqVarPVS1Helpers:
    """Helper methods for PVS1 criteria for transcript."""

    @staticmethod
    def _get_pHGVS_termination(pHGVS: str) -> int:
        """
        Get termination position from pHGVS.
        **Note:** If the position is not found, return -1.
        Examples:
        - NM_031475.2:p.Gln98*
        - NM_031475.2:p.Ala586Glyfs*73
        - NP_000305.3:p.Arg378SerfsTer5
        - p.Arg97Glyfs*26 (alternatively p.Arg97GlyfsTer26, or short p.Arg97fs)
        """
        if "fs" in pHGVS:  # If frameshift
            pattern1 = re.compile(r"p\.\D+(\d+)\D+fs(\*|X|Ter)(\d+)")
            match1 = pattern1.search(pHGVS)
            pattern2 = re.compile(r"p\.\D+(\d+)fs")
            match2 = pattern2.search(pHGVS)

            if match1:
                # if int(match1.group(1)) / (self.transcript.cds_length/3) > 0.5:
                termination = int(match1.group(1)) + int(match1.group(3))
                # else:
                #    termination = int((self.transcript.cds_length/3)/2)
            elif match2:
                termination = int(match2.group(1))
            else:
                termination = -1

        elif [char in pHGVS for char in ["*", "X", "Ter"]]:  # If premature termination codon
            pattern = re.compile(r"p\.\D+(\d+)(\*|X|Ter)")
            match = pattern.search(pHGVS)
            termination = int(match.group(1)) if match else -1
        else:
            termination = -1
        return termination

    def _undergo_nmd(self, exons: List[Exon], pHGVS: str, hgnc_id: str) -> bool:
        """
        Nonsense-mediated decay (NMD) classification. Return if the variant
        undergoes NMD prediction.
        **Rule:** If the variant is located in the last exon or in the last 50 nucleotides
        of the penultimate exon, it is NOT predicted to undergo NMD.
        """
        new_stop_codon = self._get_pHGVS_termination(pHGVS)
        cds_sizes = [exon.altEndI - exon.altStartI for exon in exons]
        if hgnc_id == "HGNC:4284":  # Hearing Loss Guidelines GJB2
            return True
        elif len(cds_sizes) <= 1:
            return False
        else:
            nmd_cutoff = sum(cds_sizes[:-1]) - min(50, cds_sizes[-2])
            return new_stop_codon * 3 <= nmd_cutoff

    @staticmethod
    def _in_biologically_relevant_transcript(transcript_tags: List[str]) -> bool:
        """
        Check if the exon with SeqVar is in a biologically relevant transcript.
        **Rule:** If the variant is located in a transcript with a MANE Select tag, it is
        predicted to be in a biologically relevant transcript.
        """
        # Ensure that necessary data is available
        return "ManeSelect" in transcript_tags

    @staticmethod
    def _critical4protein_function(seqvar: SeqVar, cds_pos: int | None, exons: List[Exon]) -> bool:
        """
        Check if the truncated/altered region is critical for the protein function.
        **Rule:** Affect of the variant on the protein function is indicated
        by experimental or clinical evidence:
        - Presence of Pathogenic variants downstream of the new stop codon
        """
        if not cds_pos:
            return False
        start_pos = 0
        for exon in exons:
            if exon.altCdsStartI < cds_pos and exon.altCdsEndI > cds_pos:
                start_pos = exon.altStartI + (cds_pos - exon.altCdsStartI)
                break
        end_pos = exons[-1].altEndI

        try:
            annonars_client = AnnonarsClient()
            response = annonars_client.get_variant_from_range(seqvar, start_pos, end_pos)
            if (
                response
                and response.result.gnomad_genomes
                and response.result.gnomad_genomes[0].vep
            ):
                pathogenic_variants = 0
                for variant in response.result.gnomad_genomes[0].vep:
                    # TODO: count only pathogenic variants
                    if variant.consequence in ["pathogenic"]:
                        pathogenic_variants += 1
                # TODO: Proove that this is the correct threshold
                if (
                    pathogenic_variants > 2
                    and pathogenic_variants / len(response.result.gnomad_genomes[0].vep) > 0.005
                ):
                    return True
                return False
            else:
                raise InvalidAPIResposeError("Failed to get variant from range.")
        except Exception as e:
            typer.secho(
                f"Failed to get variant from range for variant {seqvar.user_representation}.",
                err=True,
                fg=typer.colors.RED,
            )
            typer.secho(e, err=True, fg=typer.colors.RED)
            return False

    @staticmethod
    def _lof_is_frequent_in_population(seqvar: SeqVar, start: int, end: int) -> bool:
        """
        Check if the LoF variants in the exon are frequent in the general population.
        **Rule:** If the LoF variants in the exon > 0.1% in the general population, it is
        predicted to be frequent in the general population.
        """
        try:
            annonars_client = AnnonarsClient()
            response = annonars_client.get_variant_from_range(seqvar, start, end)
            if (
                response
                and response.result.gnomad_genomes
                and response.result.gnomad_genomes[0].vep
            ):
                lof_variants = 0
                all_variants = len(response.result.gnomad_genomes[0].vep)
                for variant in response.result.gnomad_genomes[0].vep:
                    if variant.consequence in ["frameshift_variant", "stop_gained"]:
                        lof_variants += 1
                # TODO: Proove that this is the correct threshold
                # The guideline does not specify the exact threshold. 0.1% was taken from AutoPVS1
                if lof_variants / all_variants > 0.001:
                    return True
                else:
                    return False
            else:
                raise InvalidAPIResposeError("Failed to get variant from range.")
        except Exception as e:
            typer.secho(
                f"Failed to get variant from range for variant {seqvar.user_representation}.",
                err=True,
                fg=typer.colors.RED,
            )
            typer.secho(e, err=True, fg=typer.colors.RED)
            return False

    @staticmethod
    def _lof_removes_more_then_10_percent_of_protein(pHGVS: str, exons: List[Exon]) -> bool:
        """Check if the LoF variant removes more than 10% of the protein."""
        cds_length = sum([exon.altEndI - exon.altStartI for exon in exons])
        pattern = re.compile(r"p\.\D+(\d+)(\D+fs)?(\*|X|Ter)(\d+)?")
        match = pattern.search(pHGVS)
        codon_offset = int(match.group(1)) if match else -1
        codon_length = cds_length / 3
        if codon_offset > 0 and (codon_length - codon_offset) / codon_length > 0.1:
            return True
        else:
            return False

    @staticmethod
    def _exon_skipping_or_cryptic_ss_disruption() -> bool:
        """Check if the variant causes exon skipping or cryptic splice site disruption."""
        # TODO: Implement this method
        return False

    @staticmethod
    def _alternative_start_codon() -> bool:
        """Check if the variant introduces an alternative start codon."""
        # TODO: Implement this method
        return False

    @staticmethod
    def _upstream_pathogenic_variant() -> bool:
        """Check if the variant is an upstream pathogenic variant."""
        # TODO: Implement this method
        return False


class SeqVarTranscriptsHelper:
    """Transcript information for a sequence variant."""

    def __init__(self, seqvar: SeqVar):
        self.seqvar: SeqVar = seqvar

        # Attributes to be set
        self.HGVSs: List[str] = []
        self.HGNC_id: str = ""
        self.seqvar_ts_info: List[TranscriptSeqvar] | None = None
        self.seqvar_transcript: TranscriptSeqvar | None = None
        self.gene_ts_info: List[TranscriptGene] | None = None
        self.gene_transcript: TranscriptGene | None = None
        self.consequence: SeqVarConsequence = SeqVarConsequence.NotSet

    def get_ts_info(
        self,
    ) -> Tuple[TranscriptSeqvar | None, TranscriptGene | None, SeqVarConsequence]:
        """Return the transcript information."""
        return self.seqvar_transcript, self.gene_transcript, self.consequence

    def initialize(self):
        """Get all transcripts for the given sequence variant from Mehari."""
        if not self.seqvar:
            typer.secho(
                "No sequence variant specified. Assure that the variant is resolved before fetching transcripts.",
                err=True,
                fg=typer.colors.RED,
            )
            return
        try:
            # Get transcripts from Mehari
            mehari_client = MehariClient()
            response_seqvar = mehari_client.get_seqvar_transcripts(self.seqvar)
            if not response_seqvar:
                self.seqvar_ts_info = None
            else:
                self.seqvar_ts_info = response_seqvar.result

            if not self.seqvar_ts_info or len(self.seqvar_ts_info) == 0:
                self.seqvar_transcript = None
                self.gene_transcript = None
                self.consequence = SeqVarConsequence.NotSet
                typer.secho(
                    f"No transcripts found for variant {self.seqvar.user_representation}.",
                    err=True,
                    fg=typer.colors.RED,
                )
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
                self.gene_ts_info = None
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

        except Exception as e:
            typer.secho(
                f"Failed to get transcripts for variant {self.seqvar.user_representation}.",
                err=True,
                fg=typer.colors.RED,
            )
            typer.secho(e, err=True)

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


class SeqVarPVS1(SeqVarPVS1Helpers):
    """PVS1 criteria for transcript."""

    def __init__(self, seqvar: SeqVar):
        # Attributes to be set
        self.seqvar = seqvar

        # Attributes to be computed
        self._seqvar_transcript: TranscriptSeqvar | None = None
        self._gene_transcript: TranscriptGene | None = None
        self._consequence: SeqVarConsequence = SeqVarConsequence.NotSet
        self.HGVS: str = ""
        self.pHGVS: str = ""
        self.tHGVS: str = ""
        self.HGNC_id: str = ""
        self.transcript_tags: List[str] = []
        self.exons: List[Exon] = []
        self.cds_pos: int | None = None

        # Prediction attributes
        self.prediction: PVS1Prediction = PVS1Prediction.NotPVS1
        self.prediction_path: PVS1PredictionSeqVarPath = PVS1PredictionSeqVarPath.NotSet

    def initialize(self):
        """Setup the PVS1 class."""
        # Fetch transcript data
        seqvar_transcript = SeqVarTranscriptsHelper(self.seqvar)
        seqvar_transcript.initialize()
        self._seqvar_transcript, self._gene_transcript, self._consequence = (
            seqvar_transcript.get_ts_info()
        )

        if (
            not self._seqvar_transcript
            or not self._gene_transcript
            or self._consequence == SeqVarConsequence.NotSet
        ):
            typer.secho(
                f"Failed to setup transcripts data.",
                err=True,
                fg=typer.colors.RED,
            )
            return

        # Set attributes
        self.HGVS = self._gene_transcript.id
        self.pHGVS = self.HGVS + ":" + (self._seqvar_transcript.hgvs_p or "")
        self.tHGVS = self.HGVS + ":" + (self._seqvar_transcript.hgvs_t or "")
        self.HGNC_id = self._seqvar_transcript.gene_id
        self.transcript_tags = self._seqvar_transcript.feature_tag
        self.exons = self._gene_transcript.genomeAlignments[0].exons
        self.cds_pos = (
            self._seqvar_transcript.cds_pos.ord
            if isinstance(self._seqvar_transcript.cds_pos, CdsPos)
            else None
        )

    def verify_PVS1(self):
        """Make the PVS1 prediction."""
        if (
            not self._seqvar_transcript
            or not self._gene_transcript
            or self._consequence == SeqVarConsequence.NotSet
        ):
            typer.secho(
                f"Transcript data is not set. Assure to call initialize() before verify_PVS1().",
                err=True,
                fg=typer.colors.RED,
            )
            return

        if self._consequence == SeqVarConsequence.NonsenseFrameshift:
            if self.HGNC_id == "HGNC:9588":  # Follow guidelines for PTEN
                if self._get_pHGVS_termination(self.pHGVS) < 374:
                    self.prediction = PVS1Prediction.PVS1
                    return

            if self._undergo_nmd(self.exons, self.pHGVS, self.HGNC_id):
                if self._in_biologically_relevant_transcript(self.transcript_tags):
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.NF1
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.NF2
            else:
                if self._critical4protein_function(self.seqvar, self.cds_pos, self.exons):
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionSeqVarPath.NF3
                else:
                    if self._lof_is_frequent_in_population(
                        self.seqvar,
                        self._gene_transcript.genomeAlignments[0].cdsStart,
                        self._gene_transcript.genomeAlignments[0].cdsEnd,
                    ) or not self._in_biologically_relevant_transcript(self.transcript_tags):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionSeqVarPath.NF4
                    else:
                        if self._lof_removes_more_then_10_percent_of_protein(
                            self.pHGVS, self.exons
                        ):
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionSeqVarPath.NF5
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionSeqVarPath.NF6

        elif self._consequence == SeqVarConsequence.SpliceSites:
            if self._exon_skipping_or_cryptic_ss_disruption() and self._undergo_nmd(
                self.exons, self.pHGVS, self.HGNC_id
            ):
                if self._in_biologically_relevant_transcript(self.transcript_tags):
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.SS1
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.SS2
            elif self._exon_skipping_or_cryptic_ss_disruption() and not self._undergo_nmd(
                self.exons, self.pHGVS, self.HGNC_id
            ):
                if self._critical4protein_function(self.seqvar, self.cds_pos, self.exons):
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionSeqVarPath.SS3
                else:
                    if self._lof_is_frequent_in_population(
                        self.seqvar,
                        self._gene_transcript.genomeAlignments[0].cdsStart,
                        self._gene_transcript.genomeAlignments[0].cdsEnd,
                    ) or not self._in_biologically_relevant_transcript(self.transcript_tags):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionSeqVarPath.SS4
                    else:
                        if self._lof_removes_more_then_10_percent_of_protein(
                            self.pHGVS, self.exons
                        ):
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionSeqVarPath.SS5
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionSeqVarPath.SS6
            else:
                if self._critical4protein_function(self.seqvar, self.cds_pos, self.exons):
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionSeqVarPath.SS10
                else:
                    if self._lof_is_frequent_in_population(
                        self.seqvar,
                        self._gene_transcript.genomeAlignments[0].cdsStart,
                        self._gene_transcript.genomeAlignments[0].cdsEnd,
                    ) or not self._in_biologically_relevant_transcript(self.transcript_tags):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionSeqVarPath.SS7
                    else:
                        if self._lof_removes_more_then_10_percent_of_protein(
                            self.pHGVS, self.exons
                        ):
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionSeqVarPath.SS8
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionSeqVarPath.SS9

        elif self._consequence == SeqVarConsequence.InitiationCodon:
            if self._alternative_start_codon():
                self.prediction = PVS1Prediction.NotPVS1
                self.prediction_path = PVS1PredictionSeqVarPath.IC3
            else:
                if self._upstream_pathogenic_variant():
                    self.prediction = PVS1Prediction.PVS1_Moderate
                    self.prediction_path = PVS1PredictionSeqVarPath.IC1
                else:
                    self.prediction = PVS1Prediction.PVS1_Supporting
                    self.prediction_path = PVS1PredictionSeqVarPath.IC2
        else:
            self.prediction = PVS1Prediction.NotPVS1
            typer.secho(
                f"Consequence {self._consequence} is not supported for PVS1 prediction.",
                err=True,
                fg=typer.colors.RED,
            )

    def get_prediction(self) -> Tuple[PVS1Prediction, PVS1PredictionSeqVarPath]:
        """Return the PVS1 prediction."""
        return self.prediction, self.prediction_path
