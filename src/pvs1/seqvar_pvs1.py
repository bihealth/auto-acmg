"""PVS1 criteria for Sequence Variants (SeqVar)."""

import re
from typing import Dict, List, Tuple

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.api.mehari import MehariClient
from src.defs.autopvs1 import (
    AlteredRegionMode,
    CdsInfo,
    PVS1Prediction,
    PVS1PredictionSeqVarPath,
    SeqVarConsequence,
    SeqvarConsequenceMapping,
    TranscriptInfo,
)
from src.defs.exceptions import AlgorithmError, InvalidAPIResposeError
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
        logger.debug("Getting termination position from pHGVS: {}", pHGVS)
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
        logger.debug("Termination position: {}", termination)
        return termination

    @staticmethod
    def _calculate_altered_region(
        cds_pos: int, exons: List[Exon], mode: AlteredRegionMode
    ) -> Tuple[int, int]:
        """Calculates the altered region's start and end positions.

        Args:
            cds_pos: The position of the variant in the coding sequence.
            exons: A list of exons of the gene where the variant occurs.
            mode: The mode to calculate the altered region.

        Returns:
            Tuple[int, int]: The start and end positions of the altered region.
        """
        logger.debug(
            "Calculating altered region for CDS position: {} for the mode: {}.", cds_pos, mode
        )
        if mode == AlteredRegionMode.Downstream:
            start_pos = exons[0].altStartI
            for exon in exons:
                if exon.altCdsStartI <= cds_pos <= exon.altCdsEndI:
                    start_pos = exon.altStartI + (cds_pos - exon.altCdsStartI)
                    break
                if exon.altCdsEndI < cds_pos:
                    start_pos = exon.altEndI + (cds_pos - exon.altCdsEndI)
            end_pos = exons[-1].altEndI
            logger.debug("Altered region: {} - {}", start_pos, end_pos)
            return start_pos, end_pos
        elif mode == AlteredRegionMode.Exon:
            # Get the range of the altered exon
            start_pos = exons[0].altStartI
            end_pos = exons[-1].altEndI
            for exon in exons:
                # Variant is in the exon
                if exon.altCdsStartI <= cds_pos and exon.altCdsEndI >= cds_pos:
                    start_pos = exon.altStartI + (cds_pos - exon.altCdsStartI)
                    end_pos = exon.altEndI
                    break
                # Variant is in the intron
                if exon.altCdsEndI < cds_pos:
                    start_pos = exon.altEndI + (cds_pos - exon.altCdsEndI)
                if exon.altCdsStartI > cds_pos:
                    end_pos = exon.altStartI
                    break
            logger.debug("Altered region: {} - {}", start_pos, end_pos)
            return start_pos, end_pos

    @staticmethod
    def _count_pathogenic_variants(seqvar: SeqVar, start_pos: int, end_pos: int) -> Tuple[int, int]:
        """Counts pathogenic variants in the specified range.

        Args:
            seqvar: The sequence variant being analyzed.
            start_pos: The start position of the range.
            end_pos: The end position of the range.

        Returns:
            Tuple[int, int]: The number of pathogenic variants and the total number of variants.

        Raises:
            InvalidAPIResposeError: If the API response is invalid or cannot be processed.
        """
        logger.debug("Counting pathogenic variants in the range {} - {}.", start_pos, end_pos)
        annonars_client = AnnonarsClient()
        response = annonars_client.get_variant_from_range(seqvar, start_pos, end_pos)
        if response and response.result.clinvar:
            pathogenic_variants = [
                v
                for v in response.result.clinvar
                if v.referenceAssertions
                and v.referenceAssertions[0].clinicalSignificance
                in [
                    "CLINICAL_SIGNIFICANCE_LIKELY_PATHOGENIC",
                    "CLINICAL_SIGNIFICANCE_PATHOGENIC",
                ]
            ]
            logger.debug(
                "Pathogenic variants: {}, Total variants: {}",
                len(pathogenic_variants),
                len(response.result.clinvar),
            )
            return len(pathogenic_variants), len(response.result.clinvar)
        else:
            logger.error("Failed to get variant from range. No ClinVar data.")
            raise InvalidAPIResposeError("Failed to get variant from range. No ClinVar data.")

    @staticmethod
    def _get_consequence(val: SeqVarConsequence) -> List[str]:
        """Get the VEP consequence of the sequence variant by value.

        Args:
            val: The value of the consequence.

        Returns:
            List[str]: The VEP consequences of the sequence variant.
        """
        return [key for key, value in SeqvarConsequenceMapping.items() if value == val]

    def _count_lof_variants(self, seqvar: SeqVar, start_pos: int, end_pos: int) -> Tuple[int, int]:
        """Counts Loss-of-Function (LoF) variants in the specified range.

        Args:
            seqvar: The sequence variant being analyzed.
            start_pos: The start position of the range.
            end_pos: The end position of the range.

        Returns:
            Tuple[int, int]: The number of frequent LoF variants and the total number of LoF variants.
        """
        logger.debug("Counting LoF variants in the range {} - {}.", start_pos, end_pos)
        annonars_client = AnnonarsClient()
        response = annonars_client.get_variant_from_range(seqvar, start_pos, end_pos)
        if response and response.result.gnomad_genomes:
            frequent_lof_variants = 0
            lof_variants = 0
            for variant in response.result.gnomad_genomes:
                if not variant.vep:
                    continue
                for vep in variant.vep:
                    if vep.consequence in self._get_consequence(
                        SeqVarConsequence.NonsenseFrameshift
                    ):
                        lof_variants += 1
                        if not variant.alleleCounts:
                            continue
                        for allele in variant.alleleCounts:
                            if allele.afPopmax and allele.afPopmax > 0.001:
                                frequent_lof_variants += 1
                                break
            logger.debug(
                "Frequent LoF variants: {}, Total LoF variants: {}",
                frequent_lof_variants,
                lof_variants,
            )
            return frequent_lof_variants, lof_variants
        else:
            logger.error("Failed to get variant from range. No gnomAD genomes data.")
            raise InvalidAPIResposeError(
                "Failed to get variant from range. No gnomAD genomes data."
            )

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
        logger.debug("Checking if the variant undergoes NMD.")
        new_stop_codon = self._get_pHGVS_termination(pHGVS)
        cds_sizes = [exon.altEndI - exon.altStartI for exon in exons]
        if hgnc_id == "HGNC:4284":  # Hearing Loss Guidelines GJB2
            logger.debug("Variant is in the GJB2 gene. Predicted to undergo NMD.")
            return True
        elif len(cds_sizes) <= 1:
            logger.debug("Only one (or 0) exon. Predicted to undergo NMD.")
            return False
        else:
            nmd_cutoff = sum(cds_sizes[:-1]) - min(50, cds_sizes[-2])
            logger.debug(
                "New stop codon (*3): {}, NMD cutoff: {}.",
                new_stop_codon * 3,
                nmd_cutoff,
            )
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
        logger.debug("Checking if the variant is in a biologically relevant transcript.")
        return "ManeSelect" in transcript_tags

    def _critical4protein_function(
        self, seqvar: SeqVar, cds_pos: int | None, exons: List[Exon]
    ) -> bool:
        """Checks if the truncated or altered region is critical for the protein function.

        This method assesses the impact of a sequence variant based on the presence of pathogenic
        variants downstream of the new stop codon, utilizing both experimental and clinical
        evidence.

        Note:
            The significance of a truncated or altered region is determined by the presence and
            frequency of pathogenic variants downstream from the new stop codon.
            Implementation:
                The method implements the rule by:
                - Calculating the range of the altered region.
                - Fetching variants from the specified range of the altered region.
                - Counting the number of pathogenic variants in that region, by iterating through
                  the clinvar data.
                - Considering the region critical if there are more than two pathogenic variants and
                  their frequency exceeds 5%.

        Args:
            seqvar: The sequence variant being analyzed.
            cds_pos: The position of the variant in the coding sequence.
            exons: A list of exons of the gene where the variant occurs.

        Returns:
            bool: True if the altered region is critical for the protein function, otherwise False.

        Raises:
            InvalidAPIResponseError: If the API response is invalid or cannot be processed.
        """
        logger.debug("Checking if the altered region is critical for the protein function.")
        if not cds_pos:
            logger.error("CDS position is not available. Cannot determine criticality.")
            raise AlgorithmError("CDS position is not available. Cannot determine criticality.")

        start_pos, end_pos = self._calculate_altered_region(
            cds_pos, exons, AlteredRegionMode.Downstream
        )
        try:
            pathogenic_variants, total_variants = self._count_pathogenic_variants(
                seqvar, start_pos, end_pos
            )
            if pathogenic_variants > 5 and pathogenic_variants / total_variants > 0.05:
                return True
            else:
                return False
        except Exception as e:
            logger.error("Failed to predict criticality for variant. Error: {}", e)
            raise AlgorithmError("Failed to predict criticality for variant. Error: {}", e)

    def _lof_is_frequent_in_population(
        self, seqvar: SeqVar, cds_pos: int | None, exons: List[Exon]
    ) -> bool:
        """Checks if the Loss-of-Function (LoF) variants in the exon are frequent in the general
        population.

        This function determines the frequency of LoF variants within a specified genomic region and
        evaluates whether this frequency exceeds a defined threshold indicative of common occurrence
        in the general population.

        Note:
            A LoF variant is considered frequent if its occurrence in the general population
            exceeds 0.1%. This threshold is set based on guidelines from the AutoPVS1 software.
            Implementation:
                The function implements the rule by:
                - Calculating the range of the altered region (exon with the variant).
                - Fetching variants from the specified range of the altered region.
                - Counting the number of LoF variants and frequent LoF variants in that region,
                  by iterating through the gnomAD genomes data.
                - Considering the LoF variants frequent in the general population if the frequency
                  of "frequent" LoF variants exceeds 10%.
        Args:
            seqvar: The sequence variant being analyzed.
            cds_pos: The position of the variant in the coding sequence.
            exons: A list of exons of the gene where the variant occurs.

        Returns:
            bool: True if the LoF variant frequency is greater than 0.1%, False otherwise.

        Raises:
            InvalidAPIResponseError: If the API response is invalid or cannot be processed.
        """
        logger.debug("Checking if LoF variants are frequent in the general population.")
        if not cds_pos:
            logger.error("CDS position is not available. Cannot determine LoF frequency.")
            raise AlgorithmError("CDS position is not available. Cannot determine LoF frequency.")

        start_pos, end_pos = self._calculate_altered_region(cds_pos, exons, AlteredRegionMode.Exon)
        try:
            frequent_lof_variants, lof_variants = self._count_lof_variants(
                seqvar, start_pos, end_pos
            )
            if frequent_lof_variants > 0 and frequent_lof_variants / lof_variants > 0.1:
                return True
            else:
                return False
        except Exception as e:
            logger.error("Failed to predict LoF frequency for variant. Error: {}", e)
            raise AlgorithmError("Failed to predict LoF frequency for variant. Error: {}", e)

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
        logger.debug("Checking if the LoF variant removes more than 10% of the protein.")
        cds_length = sum([exon.altEndI - exon.altStartI for exon in exons])
        pattern = re.compile(r"p\.\D+(\d+)(\D+fs)?(\*|X|Ter)(\d+)?")
        match = pattern.search(pHGVS)
        codon_offset = int(match.group(1)) if match else -1
        codon_length = cds_length / 3
        if codon_offset > 0 and (codon_length - codon_offset) / codon_length > 0.1:
            logger.debug(
                (
                    "LoF variant removes more than 10% of the protein lenght with "
                    "codon offset: {} and total length {}."
                ),
                codon_offset,
                codon_length,
            )
            return True
        else:
            logger.debug(
                (
                    "LoF variant does not remove more than 10% of the protein lenght with "
                    "codon offset: {} and total length {}."
                ),
                codon_offset,
                codon_length,
            )
            return False

    @staticmethod
    def _exon_skipping_or_cryptic_ss_disruption() -> bool:
        """Check if the variant causes exon skipping or cryptic splice site disruption."""
        # TODO: Implement this method
        return False

    @staticmethod
    def _alternative_start_codon(hgvs: str, cds_info: Dict[str, CdsInfo]) -> bool:
        """Check if the variant introduces an alternative start codon in other transcripts.

        Note:
            Rule:
                If the variant introduces an alternative start codon in other transcripts, it is
                considered to be non-pathogenic.
            Implementation:
                The rule is implemented by:
                - Iterating through all transcripts and checking if the coding sequence start
                  differs from the main transcript.
                - If the start codon differs, the rule is met.

        Args:
            hgvs: The main transcript ID.
            cds_info: A dictionary containing the CDS information for all transcripts.

        Returns:
            bool: True if the variant introduces an alternative start codon in other transcripts,
                False otherwise.
        """
        logger.debug("Checking if the variant introduces an alternative start codon.")
        if hgvs not in cds_info:
            logger.error("Main transcript ID {} not found in the dataset.", hgvs)
            raise ValueError(f"Main transcript ID {hgvs} not found in the dataset.")

        main_start_codon, main_cds_start = cds_info[hgvs].start_codon, cds_info[hgvs].cds_start
        # Check if other transcripts have alternative start codons
        alternative_starts = False
        for transcript_id, info in cds_info.items():
            if transcript_id == hgvs:
                continue
            if info.start_codon != main_start_codon and info.cds_start != main_cds_start:
                logger.debug("Alternative start codon found in transcript {}.", transcript_id)
                alternative_starts = True
                break
        return alternative_starts

    @staticmethod
    def _upstream_pathogenic_variant() -> bool:
        """Check if the transcript has an upstream pathogenic variant(s)."""
        return False


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
            mehari_client = MehariClient()
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

        except Exception as e:
            logger.error("Failed to get transcripts for the sequence variant. Error: {}", e)
            raise AlgorithmError("Failed to get transcripts for the sequence variant. Error: {}", e)

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


class SeqVarPVS1(SeqVarPVS1Helper):
    """Handles the PVS1 criteria assessment for sequence variants.

    Attributes:
        seqvar (SeqVar): The sequence variant being analyzed.
        _seqvar_transcript (TranscriptSeqvar | None): Associated transcript of the sequence variant.
        _gene_transcript (TranscriptGene | None): Associated gene transcript.
        _all_seqvar_ts (List[TranscriptSeqvar]): All sequence variant transcripts.
        _all_gene_ts (List[TranscriptGene]): All gene transcripts.
        _consequence (SeqVarConsequence): Consequence of the sequence variant.
        HGVS (str): HGVS notation of the gene.
        pHGVS (str): Protein HGVS notation.
        tHGVS (str): Transcript HGVS notation.
        HGNC_id (str): HGNC identifier for the gene.
        transcript_tags (List[str]): Tags associated with the transcript.
        exons (List[Exon]): List of exons from the gene transcript.
        cds_pos (int | None): Position of the coding sequence.
        cds_info (Dict[str, CdsInfo]): CDS information for all transcripts.
        prediction (PVS1Prediction): Prediction result based on PVS1 criteria.
        prediction_path (PVS1PredictionSeqVarPath): Pathway leading to the prediction decision.
    """

    def __init__(self, seqvar: SeqVar):
        # Attributes to be set
        self.seqvar = seqvar

        # Attributes to be computed
        self._seqvar_transcript: TranscriptSeqvar | None = None
        self._gene_transcript: TranscriptGene | None = None
        self._all_seqvar_ts: List[TranscriptSeqvar] = []
        self._all_gene_ts: List[TranscriptGene] = []
        self._consequence: SeqVarConsequence = SeqVarConsequence.NotSet
        self.HGVS: str = ""
        self.pHGVS: str = ""
        self.tHGVS: str = ""
        self.HGNC_id: str = ""
        self.transcript_tags: List[str] = []
        self.exons: List[Exon] = []
        self.cds_pos: int | None = None
        self.cds_info: Dict[str, CdsInfo] = {}

        # Prediction attributes
        self.prediction: PVS1Prediction = PVS1Prediction.NotPVS1
        self.prediction_path: PVS1PredictionSeqVarPath = PVS1PredictionSeqVarPath.NotSet

    @staticmethod
    def choose_hgvs_p(
        hgvs: str, seqvar_ts: TranscriptSeqvar, seqvar_transcripts: List[TranscriptSeqvar]
    ) -> str:
        """Choose the most suitable protein HGVS notation.

        Args:
            hgvs: The transcript HGVS notation.
            seqvar_ts: The sequence variant transcript.
            seqvar_transcripts: A list of all sequence variant transcripts.

        Returns:
            str: The most suitable protein HGVS notation.
        """
        logger.debug("Choosing the most suitable protein HGVS notation.")
        # Return pHGVS from the main transcript
        if seqvar_ts.hgvs_p and seqvar_ts.hgvs_p not in ["", "p.?"]:
            logger.debug("Protein HGVS found in the main transcript {}.", hgvs)
            return hgvs + ":" + seqvar_ts.hgvs_p
        # Choose the first transcript with a protein HGVS
        for transcript in seqvar_transcripts:
            if transcript.hgvs_p and transcript.hgvs_p not in ["", "p.?"]:
                logger.debug("Protein HGVS found in the transcript {}.", transcript.feature_id)
                return hgvs + ":" + transcript.hgvs_p
        return hgvs + ":p.?"

    def initialize(self):
        """Setup the PVS1 class.

        Fetches the transcript data and sets the attributes. Use this method before making
        predictions.
        """
        logger.debug("Setting up the SeqVarPVS1 class.")
        # Fetch transcript data
        seqvar_transcript_helper = SeqVarTranscriptsHelper(self.seqvar)
        seqvar_transcript_helper.initialize()
        (
            self._seqvar_transcript,
            self._gene_transcript,
            self._all_seqvar_ts,
            self._all_gene_ts,
            self._consequence,
        ) = seqvar_transcript_helper.get_ts_info()

        if (
            not self._seqvar_transcript
            or not self._gene_transcript
            or self._consequence == SeqVarConsequence.NotSet
        ):
            logger.error("Transcript data is not set. Cannot initialize the PVS1 class.")
            raise AlgorithmError("Transcript data is not set. Cannot initialize the PVS1 class.")

        # Set attributes
        logger.debug("Setting up the attributes for the PVS1 class.")
        self.HGVS = self._gene_transcript.id
        self.pHGVS = self.choose_hgvs_p(self.HGVS, self._seqvar_transcript, self._all_seqvar_ts)
        self.tHGVS = self.HGVS + ":" + (self._seqvar_transcript.hgvs_t or "")
        self.HGNC_id = self._seqvar_transcript.gene_id
        self.transcript_tags = self._seqvar_transcript.feature_tag
        self.exons = self._gene_transcript.genomeAlignments[0].exons
        self.cds_pos = (
            self._seqvar_transcript.cds_pos.ord
            if isinstance(self._seqvar_transcript.cds_pos, CdsPos)
            else None
        )
        self.cds_info = {
            ts.id: CdsInfo(
                start_codon=ts.startCodon,
                stop_codon=ts.stopCodon,
                cds_start=ts.genomeAlignments[0].cdsStart,
                cds_end=ts.genomeAlignments[0].cdsEnd,
                exons=ts.genomeAlignments[0].exons,
            )
            for ts in self._all_gene_ts
        }
        logger.debug("SeqVarPVS1 initialized successfully.")

    def verify_PVS1(self):
        """Make the PVS1 prediction.

        The prediction is based on the PVS1 criteria for sequence variants. The prediction
        and prediction path is stored in the prediction and prediction_path attributes.
        """
        logger.debug("Verifying the PVS1 criteria.")
        if (
            not self._seqvar_transcript
            or not self._gene_transcript
            or self._consequence == SeqVarConsequence.NotSet
        ):
            logger.error("Transcript data is not set. Did you forget to initialize the class?")
            raise AlgorithmError(
                "Transcript data is not set. Did you forget to initialize the class?"
            )

        if self._consequence == SeqVarConsequence.NonsenseFrameshift:
            if self.HGNC_id == "HGNC:9588":  # Follow guidelines for PTEN
                if self._get_pHGVS_termination(self.pHGVS) < 374:
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.PTEN
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
                        self.seqvar, self.cds_pos, self.exons
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
                        self.seqvar, self.cds_pos, self.exons
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
                        self.seqvar, self.cds_pos, self.exons
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
            if self._alternative_start_codon(self.HGVS, self.cds_info):
                self.prediction = PVS1Prediction.NotPVS1
                self.prediction_path = PVS1PredictionSeqVarPath.IC3
            else:
                if self._upstream_pathogenic_variant():
                    self.prediction = PVS1Prediction.PVS1_Moderate
                    self.prediction_path = PVS1PredictionSeqVarPath.IC1
                else:
                    self.prediction = PVS1Prediction.PVS1_Supporting
                    self.prediction_path = PVS1PredictionSeqVarPath.IC2

        elif self._consequence == SeqVarConsequence.NotSet:
            self.prediction = PVS1Prediction.NotPVS1
            self.prediction_path = PVS1PredictionSeqVarPath.NotSet
            logger.error("Consequence of the sequence variant is not set.")
            raise AlgorithmError("Consequence of the sequence variant is not set.")

    def get_prediction(self) -> Tuple[PVS1Prediction, PVS1PredictionSeqVarPath]:
        """Return the PVS1 prediction.

        Returns:
            Tuple[PVS1Prediction, PVS1PredictionSeqVarPath]: The PVS1 prediction and the path
            leading to the prediction.
        """
        return self.prediction, self.prediction_path
