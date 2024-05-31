"""PVS1 criteria for Sequence Variants (SeqVar)."""

import re
from typing import Dict, List, Optional, Tuple

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.auto_acmg import SpliceType
from src.defs.auto_pvs1 import (
    CdsInfo,
    GenomicStrand,
    PVS1Prediction,
    PVS1PredictionSeqVarPath,
    SeqVarConsequence,
    SeqvarConsequenceMapping,
)
from src.defs.exceptions import (
    AlgorithmError,
    AutoAcmgBaseException,
    InvalidAPIResposeError,
    MissingDataError,
)
from src.defs.mehari import Exon, ProteinPos, TranscriptGene, TranscriptSeqvar, TxPos
from src.defs.seqvar import SeqVar
from src.utils import SeqVarTranscriptsHelper, SplicingPrediction


class SeqVarPVS1Helper:
    """Helper methods for PVS1 criteria for sequence variants."""

    def __init__(self, *, config: Optional[Config] = None):
        self.config: Config = config or Config()
        self.annonars_client = AnnonarsClient(api_base_url=self.config.api_base_url_annonars)

    # The method is outdated
    @staticmethod
    def _choose_hgvs_p(
        hgvs: str, seqvar_ts: TranscriptSeqvar, seqvar_transcripts: List[TranscriptSeqvar]
    ) -> str:
        """
        Choose the most suitable protein HGVS notation.

        This method chooses the most suitable protein HGVS notation from available transcripts.
        The method is implemented as follows:
        - If the pHGVS is available in the main transcript, return it.
        - Else, choose the first transcript with a protein HGVS which is not empty.

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

    # The method is outdated
    @staticmethod
    def _get_pHGVS_termination(pHGVS: str) -> int:
        """
        Gets the termination position from a protein HGVS (p.HGVS) notation.

        This method extracts the termination position from the protein HGVS notation if the variant
        is a frameshift or a premature termination codon.

        Implementation of the rule:
        - If the variant is a frameshift ("fs" in pHGVS), the termination position is calculated
        by adding the frameshift position and the termination position.
        - If the variant is a premature termination codon ("*", "X", or "Ter" in pHGVS), the
        termination position is extracted from the pHGVS string.

        Note:
        If the position is not found, returns -1.

        Args:
            pHGVS: A string containing the protein HGVS notation.

        Returns:
            int: The termination position extracted from the pHGVS string, or -1 if not found.
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

    # The method is outdated
    @staticmethod
    def _get_cds_position(tHGVS: str) -> int:
        """
        Get the variant position from a transcript HGVS notation.

        This method extracts the variant position from the tHGVS notation.
        The method parses the tHGVS string using regular expressions and returns the position.

        Note:
            If the position is not found, returns -1.

        Args:
            tHGVS: A string containing the coding HGVS notation.

        Returns:
            int: The variant position extracted from the tHGVS string, or -1 if not found.
        """
        logger.debug("Getting variant position from tHGVS: {}", tHGVS)
        patterns = [
            re.compile(r"c\.(\d+)([\+\-]\d+)?(\D+)?"),  # General pattern
            re.compile(r"c\.(\*\d+)([\+\-]\d+)?(\D+)?"),  # Patterns with *
            re.compile(r"c\.(\d+_\d+)(\D+)?"),  # Patterns with range (e.g., c.123_456del)
            re.compile(r"c\.\*(\d+)([\+\-]?\d+)?(\D+)?"),  # Patterns with *
        ]

        position = -1  # default position if no match found

        for pattern in patterns:
            match = pattern.search(tHGVS)
            if match:
                # Extract numeric position ignoring other characters
                numeric_part = re.findall(r"\d+", match.group(1))
                if numeric_part:
                    position = int(numeric_part[0])
                break
        return position

    # The method is outdated
    @staticmethod
    def _calculate_5_prime_UTR_length(
        exons: List[Exon], cds_start: int, cds_end: int, strand: GenomicStrand
    ) -> int:
        """
        Calculates the length of the 5' UTR region.

        The method iterates through the exons of the gene and adds the length of the exons to the
        offset position until the coding sequence start position is reached.

        Args:
            exons: A list of exons of the gene.
            cds_start: The start position of the coding sequence.
            cds_end: The end position of the coding sequence.
            strand: The genomic strand of the gene.

        Returns:
            int: The length of the 5' UTR region.
        """
        logger.debug("Calculating the length of the 5' UTR region. Strand: {}", strand)
        utr_length = 0
        if strand == GenomicStrand.Plus:
            for exon in exons:
                if exon.altStartI < cds_start and exon.altEndI < cds_start:
                    utr_length += exon.altEndI - exon.altStartI
                elif exon.altStartI < cds_start < exon.altEndI:
                    utr_length += cds_start - exon.altStartI
                    break
        elif strand == GenomicStrand.Minus:
            for exon in exons[::-1]:
                if exon.altStartI > cds_end and exon.altEndI > cds_end:
                    utr_length += exon.altEndI - exon.altStartI
                elif exon.altStartI < cds_end < exon.altEndI:
                    utr_length += exon.altEndI - cds_end
                    break
        return utr_length

    @staticmethod
    def _calculate_altered_region(
        var_pos: int, exons: List[Exon], strand: GenomicStrand
    ) -> Tuple[int, int]:
        """
        Calculates the altered region's start and end positions.

        This method calculates the start and end positions of the altered region based on the
        position of the variant in the coding sequence and the exons of the gene.
        The method is implemented as follows:
        - If the genomic strand is plus, the start position is the variant position, and the end
        position is the last exon's end position.
        - If the genomic strand is minus, the start position is the first exon's start position,
        and the end position is the variant position.

        Args:
            var_pos: The position of the variant in the coding sequence.
            exons: A list of exons of the gene where the variant occurs.
            strand: The genomic strand of the gene.

        Returns:
            Tuple[int, int]: The start and end positions of the altered region.
        """
        logger.debug("Calculating altered region for variant position: {}.", var_pos)
        start_pos = var_pos if strand == GenomicStrand.Plus else exons[0].altStartI
        end_pos = exons[-1].altEndI if strand == GenomicStrand.Plus else var_pos
        logger.debug("Altered region: {} - {}", start_pos, end_pos)
        return start_pos, end_pos

    def _count_pathogenic_variants(
        self, seqvar: SeqVar, start_pos: int, end_pos: int
    ) -> Tuple[int, int]:
        """
        Counts pathogenic variants in the specified range.

        The method retrieves variants from the specified range and iterates through the ClinVar data
        of each variant to count the number of pathogenic variants and the total number of variants.

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
        if end_pos < start_pos:
            logger.error("End position is less than the start position.")
            logger.debug("Positions given: {} - {}", start_pos, end_pos)
            raise AlgorithmError("End position is less than the start position.")

        response = self.annonars_client.get_variant_from_range(seqvar, start_pos, end_pos)
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
        """
        Get the VEP consequence of the sequence variant by value.

        Args:
            val: The value of the consequence.

        Returns:
            List[str]: The VEP consequences of the sequence variant.
        """
        return [key for key, value in SeqvarConsequenceMapping.items() if value == val]

    def _count_lof_variants(self, seqvar: SeqVar, start_pos: int, end_pos: int) -> Tuple[int, int]:
        """
        Counts Loss-of-Function (LoF) variants in the specified range.

        The method retrieves variants from the specified range and iterates through the available
        data of each variant. The method counts the number of LoF variants and the number of
        frequent LoF variants in the specified range, based on the gnomAD genomes data (for the
        consequence of Nonsense and Frameshift variants) and the allele frequency (for the frequency
        of the LoF variants in the general population).

        Note:
            A LoF variant is considered frequent if its occurrence in the general population exceeds
            some threshold. We use a threshold of 0.1% to determine if the LoF variant is frequent.

        Args:
            seqvar: The sequence variant being analyzed.
            start_pos: The start position of the range.
            end_pos: The end position of the range.

        Returns:
            Tuple[int, int]: The number of frequent LoF variants and the total number of LoF variants.
        """
        logger.debug("Counting LoF variants in the range {} - {}.", start_pos, end_pos)
        if end_pos < start_pos:
            logger.error("End position is less than the start position.")
            logger.debug("Positions given: {} - {}", start_pos, end_pos)
            raise AlgorithmError("End position is less than the start position.")

        response = self.annonars_client.get_variant_from_range(seqvar, start_pos, end_pos)
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

    @staticmethod
    def _undergo_nmd(
        var_pos: int,
        hgnc_id: str,
        strand: Optional[GenomicStrand],
        exons: List[Exon],
    ) -> bool:
        """
        Classifies if the variant undergoes Nonsense-mediated decay (NMD).

        Implementation of the rule:
        The method checks if the variant is in the GJB2 gene and always predicts it to undergo NMD.
        If the variant is not in the GJB2 gene, the method checks if the new stop codon position
        including the 5' UTR length is less or equal to the NMD cutoff, and if it is indeed less or
        equal, the variant is predicted to undergo NMD. Otherwise, the variant is predicted to
        escape NMD.

        Note:
            Rule:
                If the variant is located in the last exon or in the last 50 nucleotides of the
                penultimate exon, it is NOT predicted to undergo NMD.
            Important:
                For the GJB2 gene (HGNC:4284), the variant is always predicted to undergo NMD.

        Args:
            var_pos: The position of the new stop codon !including the 5' UTR length!.
            hgnc_id: The HGNC ID of the gene.
            strand: The genomic strand of the gene.
            exons: A list of exons of the gene.

        Returns:
            bool: True if the variant undergoes NMD, False if variant escapes NMD.
        """
        logger.debug("Checking if the variant undergoes NMD.")
        if hgnc_id == "HGNC:4284":  # Hearing Loss Guidelines GJB2
            logger.debug("Variant is in the GJB2 gene. Predicted to undergo NMD.")
            return True
        if not strand or not exons:
            logger.error("Strand information or exons are not available. Cannot determine NMD.")
            raise MissingDataError(
                "Strand information or exons are not available. Cannot determine NMD."
            )

        tx_sizes = [exon.altCdsEndI - exon.altCdsStartI + 1 for exon in exons]
        if strand == GenomicStrand.Minus:
            tx_sizes = tx_sizes[::-1]  # Reverse the exons

        if len(tx_sizes) == 1:
            logger.debug("The only exon. Predicted to undergo NMD.")
            return False

        nmd_cutoff = sum(tx_sizes[:-1]) - min(50, tx_sizes[-2])
        logger.debug(
            "New stop codon: {}, NMD cutoff: {}.",
            var_pos,
            nmd_cutoff,
        )
        return var_pos <= nmd_cutoff

    @staticmethod
    def _in_biologically_relevant_transcript(transcript_tags: List[str]) -> bool:
        """
        Checks if the exon with SeqVar is in a biologically relevant transcript.

        Implementation of the rule:
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
        self,
        seqvar: SeqVar,
        exons: List[Exon],
        strand: Optional[GenomicStrand],
    ) -> bool:
        """
        Checks if the truncated or altered region is critical for the protein function.

        This method assesses the impact of a sequence variant based on the presence of pathogenic
        variants downstream of the new stop codon, utilizing both experimental and clinical
        evidence.

        Implementation of the rule:
        - Calculating the range of the altered region, based on the position of the variant in
        the coding sequence and the exons of the gene.
        - Fetching variants from the specified range of the altered region.
        - Counting the number of pathogenic variants in that region, by iterating through
        the clinvar data of each variant.
        - Considering the region critical if the frequency of pathogenic variants exceeds 5%.

        Note:
        The significance of a truncated or altered region is determined by the presence and
        frequency of pathogenic variants downstream from the new stop codon. We use a threshold
        of 5% to determine if the region is critical for the protein function.

        Args:
            seqvar: The sequence variant being analyzed.
            cds_pos: The position of the variant in the coding sequence.
            exons: A list of exons of the gene where the variant occurs.
            strand: The genomic strand of the gene.

        Returns:
            bool: True if the altered region is critical for the protein function, otherwise False.

        Raises:
            InvalidAPIResponseError: If the API response is invalid or cannot be processed.
        """
        logger.debug("Checking if the altered region is critical for the protein function.")
        if not strand or not exons:
            logger.error(
                "Genomic strand or exons are not available. " "Cannot determine criticality."
            )
            raise MissingDataError(
                "Genomic strand or exons are not available. " "Cannot determine criticality."
            )

        try:
            start_pos, end_pos = self._calculate_altered_region(seqvar.pos, exons, strand)
            pathogenic_variants, total_variants = self._count_pathogenic_variants(
                seqvar, start_pos, end_pos
            )
            if total_variants == 0:  # Avoid division by zero
                return False
            if pathogenic_variants / total_variants > 0.05:
                return True
            else:
                return False
        except AutoAcmgBaseException as e:
            logger.error("Failed to predict criticality for variant. Error: {}", e)
            raise AlgorithmError("Failed to predict criticality for variant.") from e

    @staticmethod
    def _find_affected_exon_position(var_pos: int, exons: List[Exon]) -> Tuple[int, int]:
        """
        Find start and end positions of the affected exon.

        Args:
            var_pos: The position of the variant in the coding sequence.
            exons: A list of exons of the gene where the variant occurs.

        Returns:
            Tuple[int, int]: The start and end positions of the affected exon.

        Raises:
            AlgorithmError: If the affected exon is not found.
        """
        logger.debug("Finding the affected exon position.")
        for exon in exons:
            if exon.altStartI <= var_pos <= exon.altEndI:
                return exon.altStartI, exon.altEndI
        logger.debug("Affected exon not found. Variant position: {}. Exons: {}", var_pos, exons)
        raise AlgorithmError("Affected exon not found.")

    def _lof_is_frequent_in_population(
        self,
        seqvar: SeqVar,
        exons: List[Exon],
        strand: Optional[GenomicStrand],
    ) -> bool:
        """
        Checks if the Loss-of-Function (LoF) variants in the exon are frequent in the general
        population.

        This function determines the frequency of LoF variants within a specified genomic region and
        evaluates whether this frequency exceeds a defined threshold indicative of common occurrence
        in the general population.

        Implementation of the rule:
        - Calculating the range of the altered region (coding sequence of the transcript).
        - Counting the number of LoF variants and frequent LoF variants in that region.
        - Considering the LoF variants frequent in the general population if the frequency
        of "frequent" LoF variants exceeds 10%.

        Note:
        A LoF variant is considered frequent if its occurrence in the general population
        exceeds some threshold. We use a threshold of 0.1% to determine if the LoF variant is
        frequent.

        Args:
            seqvar: The sequence variant being analyzed.
            exons: A list of exons of the gene where the variant occurs.
            strand: The genomic strand of the gene.

        Returns:
            bool: True if the LoF variant frequency is greater than 0.1%, False otherwise.

        Raises:
            InvalidAPIResponseError: If the API response is invalid or cannot be processed.
        """
        logger.debug("Checking if LoF variants are frequent in the general population.")
        if not strand:
            logger.error("Genomic strand is not available. Cannot determine LoF frequency.")
            raise MissingDataError(
                "Genomic strand position is not available. Cannot determine LoF frequency."
            )

        try:
            start_pos, end_pos = self._find_affected_exon_position(seqvar.pos, exons)
            frequent_lof_variants, lof_variants = self._count_lof_variants(
                seqvar, start_pos, end_pos
            )
            if lof_variants == 0:  # Avoid division by zero
                return False
            if frequent_lof_variants / lof_variants > 0.1:
                return True
            else:
                return False
        except AutoAcmgBaseException as e:
            logger.error("Failed to predict LoF frequency for variant. Error: {}", e)
            raise AlgorithmError("Failed to predict LoF frequency for variant.") from e

    @staticmethod
    def _lof_removes_more_then_10_percent_of_protein(prot_pos: int, prot_length: int) -> bool:
        """
        Check if the LoF variant removes more than 10% of the protein.

        The method checks if the LoF variant removes more than 10% of the protein based on the
        position of the variant in the protein and the length of the protein.

        Note:
            Rule:
                A LoF variant is considered to remove more than 10% of the protein if the variant
                removes more than 10% of the protein.

        Args:
            prot_pos: The position of the variant in the protein.
            prot_length: The length of the protein.

        Returns:
            bool: True if the LoF variant removes more than 10% of the protein, False otherwise.
        """
        logger.debug("Checking if the LoF variant removes more than 10% of the protein.")
        return prot_pos / prot_length > 0.1

    @staticmethod
    def _skipping_exon_pos(seqvar: SeqVar, exons: List[Exon]) -> Tuple[int, int]:
        """
        Calculate the length of the closest to the seqvar exon.

        The method calculates the length of the exon, which can be skipped due to the variant
        consequences.

        Args:
            seqvar: The sequence variant being analyzed.
            exons: A list of exons of the gene where the variant occurs.

        Returns:
            Tuple[int, int]: The start and end positions of the exon skipping region.
        """
        logger.debug("Calculating the length of the exon skipping region.")
        start_pos, end_pos = None, None
        for exon in exons:
            # Include 9 nucleotides upstream and 23 nucleotides downstream of the exon
            if exon.altStartI - 9 <= seqvar.pos <= exon.altEndI + 23:
                start_pos = exon.altStartI
                end_pos = exon.altEndI
                break
        if not start_pos or not end_pos:
            logger.error("Exon not found. Variant position: {}. Exons: {}", seqvar.pos, exons)
            raise AlgorithmError("Exon not found.")
        return start_pos, end_pos

    def _exon_skipping_or_cryptic_ss_disruption(
        self,
        seqvar: SeqVar,
        exons: List[Exon],
        consequences: List[str],
        strand: Optional[GenomicStrand],
    ) -> bool:
        """
        Check if the variant causes exon skipping or cryptic splice site disruption.

        The method checks if the variant causes exon skipping or cryptic splice site disruption
        based on the position of the variant in the coding sequence and the exons of the gene.

        Implementation of the rule:
        - If the exon length is not a multiple of 3, the variant is predicted to cause exon
        skipping.
        - If the variant is a splice acceptor or donor variant, the method predicts cryptic splice
        site disruption.

        Note:
            Rule:
            If the variant causes exon skipping or cryptic splice site disruption, it is considered
            to be pathogenic.

        Args:
            seqvar: The sequence variant being analyzed.
            exons: A list of exons of the gene where the variant occurs.
            consequences: A list of VEP consequences of the sequence variant.

        Returns:
            bool: True if the variant causes exon skipping or cryptic splice site disruption,
                False if preserves reading frame.
        """
        logger.debug(
            "Checking if the variant causes exon skipping or cryptic splice site disruption."
        )
        if not strand:
            logger.error("Strand is not available. Cannot determine exon skipping.")
            raise MissingDataError("Strand is not available. Cannot determine exon skipping.")
        start_pos, end_pos = self._skipping_exon_pos(seqvar, exons)
        if (end_pos - start_pos) % 3 != 0:
            logger.debug("Exon length is not a multiple of 3. Predicted to cause exon skipping.")
            return True

        # Cryptic splice site disruption
        sp = SplicingPrediction(seqvar, consequences=consequences, strand=strand)
        refseq = sp.get_sequence(seqvar.pos - 20, seqvar.pos + 20)
        splice_type = sp.determine_splice_type(consequences)
        cryptic_sites = sp.get_cryptic_ss(refseq, splice_type)
        if len(cryptic_sites) > 0:
            for site in cryptic_sites:
                if abs(site[0] - seqvar.pos) % 3 != 0:
                    logger.debug("Cryptic splice site disruption predicted.")
                    return True
        return False

    @staticmethod
    def _closest_alternative_start_codon(cds_info: Dict[str, CdsInfo], hgvs: str) -> Optional[int]:
        """
        Calculate the closest potential start codon.

        The method calculates the closest potential start codon based on the position of the variant
        in the coding sequence and the CDS information of the gene.

        Args:
            cds_info: A dictionary containing the CDS information for all transcripts.
            hgvs: The main transcript ID.

        Returns:
            Optional[int]: The position of the closest potential start codon, or None if not found.
        """
        logger.debug("Checking if the variant introduces an alternative start codon.")
        if hgvs not in cds_info:  # Should never happen
            logger.error("Main transcript ID {} not found in the transcripts data.", hgvs)
            raise MissingDataError(f"Main transcript ID {hgvs} not found in the transcripts data.")

        main_cds_start = (
            cds_info[hgvs].cds_start
            if cds_info[hgvs].cds_strand == GenomicStrand.Plus
            else cds_info[hgvs].cds_end
        )
        main_strand = cds_info[hgvs].cds_strand

        closest_alternative_start = None
        for transcript_id, info in cds_info.items():
            if transcript_id == hgvs:
                continue
            alt_cds_start = (
                info.cds_start if info.cds_strand == GenomicStrand.Plus else info.cds_end
            )
            if alt_cds_start != main_cds_start and main_strand == info.cds_strand:
                if (
                    not closest_alternative_start  # New start codon found
                    or alt_cds_start - main_cds_start
                    < closest_alternative_start - main_cds_start  # Closer start codon found
                ):
                    closest_alternative_start = alt_cds_start
        return closest_alternative_start

    def _alternative_start_codon(self, cds_info: Dict[str, CdsInfo], hgvs: str) -> bool:
        """
        Check if the variant introduces an alternative start codon in other transcripts.

        Implementation of the rule:
            - Iterating through all transcripts and checking if the coding sequence start
                differs from the main transcript.
            - If the start codon differs, the rule is met.

        Note:
            Rule:
                If the variant introduces an alternative start codon in other transcripts, it is
                considered to be non-pathogenic.

        Args:
            hgvs: The main transcript ID.
            cds_info: A dictionary containing the CDS information for all transcripts.

        Returns:
            bool: True if the variant introduces an alternative start codon in other transcripts,
                False otherwise.
        """
        logger.debug("Checking if the variant introduces an alternative start codon.")
        alt_start = self._closest_alternative_start_codon(cds_info, hgvs)
        return alt_start is not None

    def _upstream_pathogenic_variants(
        self,
        seqvar: SeqVar,
        exons: List[Exon],
        strand: Optional[GenomicStrand],
        cds_info: Dict[str, CdsInfo],
        hgvs: str,
    ) -> bool:
        """
        Look for pathogenic variants upstream of the closest potential in-frame start codon.

        The method checks for pathogenic variants upstream of the closest potential in-frame start
        codon. The method is implemented as follows:
        - Find the closest potential in-frame start codon.
        - Fetch and count pathogenic variants in the specified range.
        - Return True if pathogenic variants are found, otherwise False.

        Args:
            seqvar: The sequence variant being analyzed.
            cds_pos: The position of the variant in the coding sequence.
            exons: A list of exons of the gene where the variant occurs.
            strand: The genomic strand of the gene.
            cds_info: A dictionary containing the CDS information for all transcripts.
            hgvs: The main transcript ID.

        Returns:
            bool: True if pathogenic variants are found upstream of the closest potential in-frame
                start codon, False otherwise.
        """
        logger.debug(
            "Checking for pathogenic variants upstream of the closest in-frame start codon."
        )

        if not strand:
            logger.error("Strand is not available. Cannot determine upstream pathogenic variants.")
            raise MissingDataError(
                "Strand is not available. Cannot determine upstream pathogenic variants."
            )
        start_pos, end_pos = None, None
        if strand == GenomicStrand.Plus:
            start_pos = exons[0].altStartI
            end_pos = self._closest_alternative_start_codon(cds_info, hgvs)
        elif strand == GenomicStrand.Minus:
            start_pos = self._closest_alternative_start_codon(cds_info, hgvs)
            end_pos = exons[-1].altEndI

        if not start_pos or not end_pos:
            logger.debug("No alternative start codon found. Skipping upstream pathogenic variants.")
            raise AlgorithmError(
                "No alternative start codon found. Skipping upstream pathogenic variants."
            )

        # Fetch and count pathogenic variants in the specified range
        try:
            pathogenic_variants, _ = self._count_pathogenic_variants(seqvar, start_pos, end_pos)
            return pathogenic_variants > 0
        except AutoAcmgBaseException as e:
            logger.error("Failed to check upstream pathogenic variants. Error: {}", e)
            raise AlgorithmError("Failed to check upstream pathogenic variants.") from e


class SeqVarPVS1(SeqVarPVS1Helper):
    """Handles the PVS1 criteria assessment for sequence variants."""

    def __init__(self, seqvar: SeqVar, *, config: Optional[Config] = None):
        #: Configuration to use.
        self.config: Config = config or Config()
        self.annonars_client = AnnonarsClient(api_base_url=self.config.api_base_url_annonars)
        # Attributes to be set
        self.seqvar = seqvar
        # Attributes to be computed
        self._seqvar_transcript: TranscriptSeqvar | None = None
        self._gene_transcript: TranscriptGene | None = None
        self._all_seqvar_ts: List[TranscriptSeqvar] = []
        self._all_gene_ts: List[TranscriptGene] = []
        self._consequence: SeqVarConsequence = SeqVarConsequence.NotSet
        self.HGVS: str = ""
        self.consequence: str = ""
        # self.pHGVS: str = ""
        # self.tHGVS: str = ""
        self.HGNC_id: str = ""
        self.transcript_tags: List[str] = []
        self.exons: List[Exon] = []
        # self.cds_pos: int = -1
        self.tx_pos_utr: int = -1
        self.prot_pos: int = -1
        self.prot_length: int = -1
        self.cds_info: Dict[str, CdsInfo] = {}
        # self.cds_start: int = -1
        # self.cds_end: int = -1
        # self.cds_length: int = -1
        self.strand: Optional[GenomicStrand] = None
        # Prediction attributes
        self.prediction: PVS1Prediction = PVS1Prediction.NotPVS1
        self.prediction_path: PVS1PredictionSeqVarPath = PVS1PredictionSeqVarPath.NotSet

    def initialize(self):
        """Setup the PVS1 class.

        Fetches the transcript data and sets the attributes. Use this method before making
        predictions.

        Raises:
            MissingDataError: If the transcript data is not fully set or if some attributes could
            not be set.
        """
        logger.debug("Setting up the SeqVarPVS1 class.")
        # Fetch transcript data
        seqvar_transcript_helper = SeqVarTranscriptsHelper(self.seqvar, config=self.config)
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
            raise MissingDataError(
                "Transcript data is not fully set. Cannot initialize the PVS1 class."
            )

        # Set attributes
        logger.debug("Setting up the attributes for the PVS1 class.")
        self.consequences = self._seqvar_transcript.consequences
        self.HGVS = self._gene_transcript.id
        # self.pHGVS = self._choose_hgvs_p(self.HGVS, self._seqvar_transcript, self._all_seqvar_ts)
        # self.tHGVS = self._seqvar_transcript.hgvs_t if self._seqvar_transcript.hgvs_t else ""
        self.HGNC_id = self._seqvar_transcript.gene_id
        self.transcript_tags = self._seqvar_transcript.feature_tag
        self.exons = self._gene_transcript.genomeAlignments[0].exons
        # self.cds_pos = (
        #     self._seqvar_transcript.cds_pos.ord
        #     if isinstance(self._seqvar_transcript.cds_pos, CdsPos)
        #     else -1
        # )
        self.tx_pos_utr = (
            self._seqvar_transcript.tx_pos.ord
            if isinstance(self._seqvar_transcript.tx_pos, TxPos)
            else -1
        )
        self.prot_pos = (
            self._seqvar_transcript.protein_pos.ord
            if isinstance(self._seqvar_transcript.protein_pos, ProteinPos)
            else -1
        )
        self.prot_length = (
            self._seqvar_transcript.protein_pos.total
            if isinstance(self._seqvar_transcript.protein_pos, ProteinPos)
            else -1
        )
        self.cds_info = {
            ts.id: CdsInfo(
                start_codon=ts.startCodon,
                stop_codon=ts.stopCodon,
                cds_start=ts.genomeAlignments[0].cdsStart,
                cds_end=ts.genomeAlignments[0].cdsEnd,
                cds_strand=GenomicStrand.from_string(ts.genomeAlignments[0].strand),
                exons=ts.genomeAlignments[0].exons,
            )
            for ts in self._all_gene_ts
        }
        # self.cds_start = self._gene_transcript.genomeAlignments[0].cdsStart
        # self.cds_end = self._gene_transcript.genomeAlignments[0].cdsEnd
        # self.cds_length = self._gene_transcript.stopCodon - self._gene_transcript.startCodon
        self.strand = GenomicStrand.from_string(self._gene_transcript.genomeAlignments[0].strand)
        if (
            not self.exons
            # # or self.cds_pos == -1
            # or self.tx_pos_utr == -1
            # or self.prot_pos == -1
            # or self.prot_length == -1
            # # or self.cds_start == -1
            # # or self.cds_end == -1
            # # or self.cds_length == -1
            # or not self.strand
        ):
            logger.error("Some attributes are not set. Cannot initialize the PVS1 class.")
            logger.debug(
                "Exons: {}, CDS position: {}, TX position UTR: {}, Protein position: {}, "
                "Protein length: {}, CDS start: {}, CDS end: {}, CDS length: {}, Strand: {}.",
                self.exons,
                # self.cds_pos,
                self.tx_pos_utr,
                self.prot_pos,
                self.prot_length,
                # self.cds_start,
                # self.cds_end,
                # self.cds_length,
                self.strand,
            )
            raise MissingDataError("Some attributes are not set. Cannot initialize the PVS1 class.")
        logger.debug("SeqVarPVS1 initialized successfully.")

    def verify_PVS1(self) -> Tuple[PVS1Prediction, PVS1PredictionSeqVarPath]:
        """Make the PVS1 prediction.

        The prediction is based on the PVS1 criteria for sequence variants. The prediction
        and prediction path is stored in the prediction and prediction_path attributes.

        Returns:
            Tuple[PVS1Prediction, PVS1PredictionSeqVarPath]: The PVS1 prediction and the prediction
            path.
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
                if self.prot_pos < 374:
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.PTEN
                    return self.prediction, self.prediction_path

            if self._undergo_nmd(self.tx_pos_utr, self.HGNC_id, self.strand, self.exons):
                if self._in_biologically_relevant_transcript(self.transcript_tags):
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.NF1
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.NF2
            else:
                if self._critical4protein_function(self.seqvar, self.exons, self.strand):
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionSeqVarPath.NF3
                else:
                    if self._lof_is_frequent_in_population(
                        self.seqvar, self.exons, self.strand
                    ) or not self._in_biologically_relevant_transcript(self.transcript_tags):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionSeqVarPath.NF4
                    else:
                        if self._lof_removes_more_then_10_percent_of_protein(
                            self.prot_pos, self.prot_length
                        ):
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionSeqVarPath.NF5
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionSeqVarPath.NF6

        elif self._consequence == SeqVarConsequence.SpliceSites:
            if self._exon_skipping_or_cryptic_ss_disruption(
                self.seqvar, self.exons, self.consequences, self.strand
            ) and self._undergo_nmd(self.tx_pos_utr, self.HGNC_id, self.strand, self.exons):
                if self._in_biologically_relevant_transcript(self.transcript_tags):
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.SS1
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.SS2
            elif self._exon_skipping_or_cryptic_ss_disruption(
                self.seqvar, self.exons, self.consequences, self.strand
            ) and not self._undergo_nmd(self.tx_pos_utr, self.HGNC_id, self.strand, self.exons):
                if self._critical4protein_function(self.seqvar, self.exons, self.strand):
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionSeqVarPath.SS3
                else:
                    if self._lof_is_frequent_in_population(
                        self.seqvar, self.exons, self.strand
                    ) or not self._in_biologically_relevant_transcript(self.transcript_tags):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionSeqVarPath.SS4
                    else:
                        if self._lof_removes_more_then_10_percent_of_protein(
                            self.prot_pos, self.prot_length
                        ):
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionSeqVarPath.SS5
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionSeqVarPath.SS6
            else:
                if self._critical4protein_function(self.seqvar, self.exons, self.strand):
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionSeqVarPath.SS10
                else:
                    if self._lof_is_frequent_in_population(
                        self.seqvar, self.exons, self.strand
                    ) or not self._in_biologically_relevant_transcript(self.transcript_tags):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionSeqVarPath.SS7
                    else:
                        if self._lof_removes_more_then_10_percent_of_protein(
                            self.prot_pos, self.prot_length
                        ):
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionSeqVarPath.SS8
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionSeqVarPath.SS9

        elif self._consequence == SeqVarConsequence.InitiationCodon:
            if self._alternative_start_codon(self.cds_info, self.HGVS):
                self.prediction = PVS1Prediction.NotPVS1
                self.prediction_path = PVS1PredictionSeqVarPath.IC3
            else:
                if self._upstream_pathogenic_variants(
                    self.seqvar, self.exons, self.strand, self.cds_info, self.HGVS
                ):
                    self.prediction = PVS1Prediction.PVS1_Moderate
                    self.prediction_path = PVS1PredictionSeqVarPath.IC1
                else:
                    self.prediction = PVS1Prediction.PVS1_Supporting
                    self.prediction_path = PVS1PredictionSeqVarPath.IC2

        else:
            self.prediction = PVS1Prediction.NotPVS1
            self.prediction_path = PVS1PredictionSeqVarPath.NotSet
            logger.info(
                "PVS1 criteria is not met, since the variant consequence is {}.", self._consequence
            )

        return self.prediction, self.prediction_path
