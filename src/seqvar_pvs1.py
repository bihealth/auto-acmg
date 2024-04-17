"""PVS1 criteria for Sequence Variants (SeqVar)."""

import json
import re
from typing import Dict, List, Tuple

import typer

from src.api.annonars import AnnonarsClient
from src.api.mehari import MehariClient
from src.defs.autopvs1 import (
    PVS1Prediction,
    PVS1PredictionSeqVarPath,
    SeqVarConsequence,
    SeqvarConsequenceMapping,
    TranscriptInfo,
)
from src.defs.exceptions import InvalidAPIResposeError
from src.defs.mehari import CdsPos, Exon, TranscriptGene, TranscriptSeqvar
from src.defs.seqvar import SeqVar


class SeqVarPVS1Helper:
    """Helper methods for PVS1 criteria for sequence variants."""

    @staticmethod
    def _get_pHGVS_termination(pHGVS: str) -> int:
        """Gets the termination position from a protein HGVS (p.HGVS) notation.

        Note:
            If the position is not found, returns -1.

        Args:
            pHGVS: A string containing the protein HGVS notation.

        Returns:
            int: The termination position extracted from the pHGVS string, or -1 if not found.

        Examples:
            >>> get_termination_position("NM_031475.2:p.Gln98*")
            98
            >>> get_termination_position("NM_031475.2:p.Ala586Glyfs*73")
            586
            >>> get_termination_position("NP_000305.3:p.Arg378SerfsTer5")
            378
            >>> get_termination_position("p.Arg97Glyfs*26")
            97 + 26 = 123
            >>> get_termination_position("p.Arg97GlyfsTer26")
            97 + 26 = 123
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
        """Classifies if the variant undergoes Nonsense-mediated decay (NMD).

        Note:
            Rule:
            If the variant is located in the last exon or in the last 50 nucleotides of the
            penultimate exon, it is NOT predicted to undergo NMD.

            Important:
            For the GJB2 gene (HGNC:4284), the variant is always predicted to undergo NMD.

        Args:
            exons: A list of exons of the gene.
            pHGVS: A string containing the protein HGVS notation.
            hgnc_id: The HGNC ID of the gene.

        Returns:
            bool: True if the variant is predicted to undergo NMD, False otherwise.
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
        """Checks if the exon with SeqVar is in a biologically relevant transcript.

        Note:
            Rule:
            If the variant is located in a transcript with a MANE Select tag, it is
            considered to be in a biologically relevant transcript.

        Args:
            transcript_tags: A list of tags for the transcript.

        Returns:
            bool: True if the variant is in a biologically relevant transcript, False otherwise.
        """
        return "ManeSelect" in transcript_tags

    @staticmethod
    def _critical4protein_function(seqvar: SeqVar, cds_pos: int | None, exons: List[Exon]) -> bool:
        """Checks if the truncated or altered region is critical for the protein function.

        This method assesses the impact of a sequence variant based on the presence of pathogenic
        variants downstream of the new stop codon, utilizing both experimental and clinical evidence.

        Note:
            The significance of a truncated or altered region is determined by the presence and
            frequency of pathogenic variants downstream from the new stop codon.

        Implementation:
            The method implements the rule by:
            - Fetching variants from the range of the altered region downstream from the new position.
            - Counting the number of pathogenic variants in that region.
            - Considering the region critical if there are more than two pathogenic variants and
            their frequency exceeds 0.5%.

        Args:
            seqvar (SeqVar): The sequence variant being analyzed.
            cds_pos (int): The position of the variant in the coding sequence.
            exons (list of Exon): A list of exons of the gene where the variant occurs.

        Returns:
            bool: True if the altered region is critical for the protein function, otherwise False.

        Raises:
            InvalidAPIResponseError: If the API response is invalid or cannot be processed.
        """
        if not cds_pos:
            return False
        # Get the range of the altered region
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
                f"Failed to get variant from range for variant {seqvar.user_repr}.",
                err=True,
                fg=typer.colors.RED,
            )
            typer.secho(e, err=True, fg=typer.colors.RED)
            return False

    @staticmethod
    def _lof_is_frequent_in_population(seqvar: SeqVar, start: int, end: int) -> bool:
        """Checks if the Loss-of-Function (LoF) variants in the exon are frequent in the general population.

        This function determines the frequency of LoF variants within a specified genomic region and evaluates
        whether this frequency exceeds a defined threshold indicative of common occurrence in the general population.

        Note:
            A LoF variant is considered frequent if its occurrence in the general population exceeds 0.1%.
            This threshold is set based on guidelines from the AutoPVS1 software.

        Implementation:
            The function implements the rule by:
            - Fetching variants from the specified range of the altered region.
            - Counting the number of LoF variants, specifically those classified as 'frameshift_variant' or 'stop_gained'.
            - Determining the region as having frequent LoF variants if their frequency exceeds 0.1%.

        Args:
            seqvar (SeqVar): The sequence variant being analyzed.
            start (int): The start position of the altered region in the genomic sequence.
            end (int): The end position of the altered region in the genomic sequence.

        Returns:
            bool: True if the LoF variant frequency is greater than 0.1%, False otherwise.

        Raises:
            InvalidAPIResponseError: If the API response is invalid or cannot be processed.
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
                f"Failed to get variant from range for variant {seqvar.user_repr}.",
                err=True,
                fg=typer.colors.RED,
            )
            typer.secho(e, err=True, fg=typer.colors.RED)
            return False

    @staticmethod
    def _lof_removes_more_then_10_percent_of_protein(pHGVS: str, exons: List[Exon]) -> bool:
        """Check if the LoF variant removes more than 10% of the protein.

        Note:
            Rule:
            A LoF variant is considered to remove more than 10% of the protein if the variant
            removes more than 10% of the protein:)

            Implementation:
            The rule is implemented by:
            - Calculating the length of the coding sequence (based on pHGVS).
            - Calculating the length of the protein based on exons information.
            - If the variant removes more than 10% of the protein, the rule is met.

        Args:
            pHGVS: A string containing the protein HGVS notation.
            exons: A list of exons of the gene.

        Returns:
            bool: True if the LoF variant removes more than 10% of the protein, False otherwise.
        """
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
        """Return the transcript information.

        Returns:
            Tuple[TranscriptSeqvar | None, TranscriptGene | None, SeqVarConsequence]: The sequence variant transcript,
            gene transcript, and the consequence of the sequence variant.
        """
        return self.seqvar_transcript, self.gene_transcript, self.consequence

    def initialize(self):
        """Get all transcripts for the given sequence variant from Mehari."""
        # Should never happen, since __init__ is called with seqvar
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
                    f"No transcripts found for variant {self.seqvar.user_repr}.",
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
                f"Failed to get transcripts for variant {self.seqvar.user_repr}.",
                err=True,
                fg=typer.colors.RED,
            )
            typer.secho(e, err=True)

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
            # If no overlapping transcripts, return None
            if not lookup_group or not exon_lengths:
                return None, None
            max_length_transcript = max(lookup_group, key=lambda x: exon_lengths[x])
            seqvar_transcript = transcripts_mapping[max_length_transcript].seqvar
            gene_transcript = transcripts_mapping[max_length_transcript].gene
        return seqvar_transcript, gene_transcript


class SeqVarPVS1(SeqVarPVS1Helper):
    """Handles the PVS1 criteria assessment for sequence variants."""

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
        """Setup the PVS1 class.

        Fetches the transcript data and sets the attributes. Use this method before making
        predictions.
        """
        # Fetch transcript data
        seqvar_transcript_helper = SeqVarTranscriptsHelper(self.seqvar)
        seqvar_transcript_helper.initialize()
        self._seqvar_transcript, self._gene_transcript, self._consequence = (
            seqvar_transcript_helper.get_ts_info()
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
        """Make the PVS1 prediction.

        The prediction is based on the PVS1 criteria for sequence variants. The prediction
        and prediction path is stored in the prediction and prediction_path attributes.
        """
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
        """Return the PVS1 prediction.

        Returns:
            Tuple[PVS1Prediction, PVS1PredictionSeqVarPath]: The PVS1 prediction and the path leading to the prediction.
        """
        return self.prediction, self.prediction_path
