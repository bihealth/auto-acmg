"""PVS1 criteria for Sequence Variants (SeqVar)."""

import re
from typing import Dict, List, Optional, Tuple

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
    SpliceType,
)
from src.defs.auto_pvs1 import (
    CdsInfo,
    GenomicStrand,
    PVS1Prediction,
    PVS1PredictionPathMapping,
    PVS1PredictionSeqVarPath,
    SeqvarConsequenceMapping,
    SeqVarPVS1Consequence,
)
from src.defs.exceptions import (
    AlgorithmError,
    AutoAcmgBaseException,
    InvalidAPIResposeError,
    MissingDataError,
)
from src.defs.mehari import Exon, ProteinPos, TranscriptGene, TranscriptSeqvar, TxPos
from src.defs.seqvar import SeqVar
from src.utils import AutoACMGHelper, SeqVarTranscriptsHelper, SplicingPrediction


class SeqVarPVS1Helper(AutoACMGHelper):
    """Helper methods for PVS1 criteria for sequence variants."""

    def __init__(self):
        super().__init__()
        #: Comment to store the prediction explanation.
        self.comment_pvs1: str = ""

    def _calc_alt_reg(
        self, var_pos: int, exons: List[Exon], strand: GenomicStrand
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

    def _count_pathogenic_vars(
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
        if response and response.clinvar:
            pathogenic_variants = [
                v
                for v in response.clinvar
                if v.records
                and v.records[0].classifications
                and v.records[0].classifications.germlineClassification
                and v.records[0].classifications.germlineClassification.description
                in [
                    "Pathogenic",
                    "Likely pathogenic",
                ]
            ]
            logger.debug(
                "Pathogenic variants: {}, Total variants: {}",
                len(pathogenic_variants),
                len(response.clinvar),
            )
            return len(pathogenic_variants), len(response.clinvar)
        else:
            logger.error("Failed to get variant from range. No ClinVar data.")
            raise InvalidAPIResposeError("Failed to get variant from range. No ClinVar data.")

    def _find_aff_exon_pos(self, var_pos: int, exons: List[Exon]) -> Tuple[int, int]:
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

    def _get_conseq(self, val: SeqVarPVS1Consequence) -> List[str]:
        """
        Get the VEP consequence of the sequence variant by value.

        Args:
            val: The value of the consequence.

        Returns:
            List[str]: The VEP consequences of the sequence variant.
        """
        return [key for key, value in SeqvarConsequenceMapping.items() if value == val]

    def _count_lof_vars(self, seqvar: SeqVar, start_pos: int, end_pos: int) -> Tuple[int, int]:
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
        if response and response.gnomad_genomes:
            frequent_lof_variants = 0
            lof_variants = 0
            for variant in response.gnomad_genomes:
                if not variant.vep:
                    continue
                for vep in variant.vep:
                    if vep.consequence in self._get_conseq(
                        SeqVarPVS1Consequence.NonsenseFrameshift
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

    def _closest_alt_start_cdn(self, cds_info: Dict[str, CdsInfo], hgvs: str) -> Optional[int]:
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

    def _skipping_exon_pos(self, seqvar: SeqVar, exons: List[Exon]) -> Tuple[int, int]:
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

    def undergo_nmd(
        self,
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
            self.comment_pvs1 += (
                f"Variant is in the GJB2 gene (hgnc: {hgnc_id}). "
                "Always predicted to undergo NMD."
            )
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
            self.comment_pvs1 += "The only exon found. Predicted to undergo NMD."
            return False

        nmd_cutoff = sum(tx_sizes[:-1]) - min(50, tx_sizes[-2])
        logger.debug(
            "New stop codon: {}, NMD cutoff: {}.",
            var_pos,
            nmd_cutoff,
        )
        self.comment_pvs1 += (
            f"New stop codon: {var_pos}, NMD cutoff: {nmd_cutoff}."
            f"{'Predicted to undergo NMD.' if var_pos <= nmd_cutoff else 'Predicted to escape NMD.'}"
        )
        return var_pos <= nmd_cutoff

    def in_bio_relevant_tsx(self, transcript_tags: List[str]) -> bool:
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
        self.comment_pvs1 += (
            "Variant is in a biologically relevant transcript. "
            f"Transcript tags: {', '.join(transcript_tags)}."
        )
        return "ManeSelect" in transcript_tags

    def crit4prot_func(
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
            start_pos, end_pos = self._calc_alt_reg(seqvar.pos, exons, strand)
            pathogenic_variants, total_variants = self._count_pathogenic_vars(
                seqvar, start_pos, end_pos
            )
            self.comment_pvs1 += (
                f"Found {pathogenic_variants} pathogenic variants from {total_variants} total "
                f"variants in the range {start_pos} - {end_pos}. "
            )
            if total_variants == 0:  # Avoid division by zero
                self.comment_pvs1 += "No variants found. Predicted to be non-critical."
                return False
            if pathogenic_variants / total_variants > 0.05:
                self.comment_pvs1 += (
                    "Frequency of pathogenic variants "
                    f"{pathogenic_variants/total_variants} exceeds 5%. "
                    "Predicted to be critical."
                )
                return True
            else:
                self.comment_pvs1 += (
                    "Frequency of pathogenic variants "
                    f"{pathogenic_variants/total_variants} does not exceed 5%. "
                    "Predicted to be non-critical."
                )
                return False
        except AutoAcmgBaseException as e:
            logger.error("Failed to predict criticality for variant. Error: {}", e)
            raise AlgorithmError("Failed to predict criticality for variant.") from e

    def lof_freq_in_pop(
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
            start_pos, end_pos = self._find_aff_exon_pos(seqvar.pos, exons)
            frequent_lof_variants, lof_variants = self._count_lof_vars(seqvar, start_pos, end_pos)
            self.comment_pvs1 += (
                f"Found {frequent_lof_variants} frequent LoF variants from {lof_variants} total "
                f"LoF variants in the range {start_pos} - {end_pos}. "
            )
            if lof_variants == 0:  # Avoid division by zero
                self.comment_pvs1 += "No LoF variants found. Predicted to be non-frequent."
                return False
            if frequent_lof_variants / lof_variants > 0.1:
                self.comment_pvs1 += (
                    "Frequency of frequent LoF variants "
                    f"{frequent_lof_variants/lof_variants} exceeds 0.1%. "
                    "Predicted to be frequent."
                )
                return True
            else:
                self.comment_pvs1 += (
                    "Frequency of frequent LoF variants "
                    f"{frequent_lof_variants/lof_variants} does not exceed 0.1%. "
                    "Predicted to be non-frequent."
                )
                return False
        except AutoAcmgBaseException as e:
            logger.error("Failed to predict LoF frequency for variant. Error: {}", e)
            raise AlgorithmError("Failed to predict LoF frequency for variant.") from e

    def lof_rm_gt_10pct_of_prot(self, prot_pos: int, prot_length: int) -> bool:
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
        self.comment_pvs1 += (
            f"Variant removes {prot_pos} amino acids from the protein of length {prot_length}. "
            f"{'Predicted to remove more than 10% of the protein.' if prot_pos / prot_length > 0.1 else 'Predicted to remove less than 10% of the protein.'}"
        )
        return prot_pos / prot_length > 0.1

    def exon_skip_or_cryptic_ss_disrupt(
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
        self.comment_pvs1 += f"Variant's exon position: {start_pos} - {end_pos}."
        if (end_pos - start_pos) % 3 != 0:
            logger.debug("Exon length is not a multiple of 3. Predicted to cause exon skipping.")
            self.comment_pvs1 += (
                "Exon length is not a multiple of 3. Predicted to cause exon skipping."
            )
            return True
        else:
            self.comment_pvs1 += (
                "Exon length is a multiple of 3. Predicted to preserve reading frame "
                "for exon skipping."
            )

        # Cryptic splice site disruption
        sp = SplicingPrediction(seqvar, consequences=consequences, strand=strand, exons=exons)
        refseq = sp.get_sequence(seqvar.pos - 20, seqvar.pos + 20)
        splice_type = sp.determine_splice_type(consequences)
        cryptic_sites = sp.get_cryptic_ss(refseq, splice_type)
        if len(cryptic_sites) > 0:
            for site in cryptic_sites:
                if abs(site[0] - seqvar.pos) % 3 != 0:
                    logger.debug("Cryptic splice site disruption predicted.")
                    self.comment_pvs1 += (
                        "Cryptic splice site disruption predicted. "
                        f"Cryptic splice site: position {site[0]}, splice context {site[1]}, "
                        f"maximnum entropy score {site[2]}. "
                        f"Cryptic splice site - variant position ({seqvar.pos}) = "
                        f"{abs(site[0] - seqvar.pos)} is not devisible by 3."
                    )
                    return True
            logger.debug("Cryptic splice site disruption not predicted.")
            self.comment_pvs1 += "All cryptic splice sites preserve reading frame."
            for i, site in enumerate(cryptic_sites):
                self.comment_pvs1 += (
                    f"Cryptic splice site {i}: position {site[0]}, splice context {site[1]}, "
                    f"maximnum entropy score {site[2]}. "
                    f"Cryptic splice site - variant position ({seqvar.pos}) = "
                    f"{abs(site[0] - seqvar.pos)} is devisible by 3."
                )
        else:
            logger.debug("No cryptic splice site found. Predicted to preserve reading frame.")
            self.comment_pvs1 += (
                "No cryptic splice site found. Predicted to preserve reading frame."
            )
        return False

    def alt_start_cdn(self, cds_info: Dict[str, CdsInfo], hgvs: str) -> bool:
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
        alt_start = self._closest_alt_start_cdn(cds_info, hgvs)
        if alt_start is not None:
            self.comment_pvs1 += (
                f"Alternative start codon found at position {alt_start}. "
                "Predicted to be non-pathogenic."
            )
            return True
        self.comment_pvs1 += "No alternative start codon found."
        return False

    def up_pathogenic_vars(
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
        if strand == GenomicStrand.Plus:
            start_pos = exons[0].altStartI
            end_pos = exons[0].altEndI
        elif strand == GenomicStrand.Minus:
            start_pos = exons[-1].altStartI
            end_pos = exons[-1].altEndI

        # Fetch and count pathogenic variants in the specified range
        try:
            pathogenic_variants, _ = self._count_pathogenic_vars(seqvar, start_pos, end_pos)
            self.comment_pvs1 += (
                f"Found {pathogenic_variants} pathogenic variants upstream of the closest "
                f"in-frame start codon. The search range: {start_pos} - {end_pos} "
                f"with {'+' if strand == GenomicStrand.Plus else '-'} strand."
            )
            return pathogenic_variants > 0
        except AutoAcmgBaseException as e:
            logger.error("Failed to check upstream pathogenic variants. Error: {}", e)
            raise AlgorithmError("Failed to check upstream pathogenic variants.") from e


class AutoPVS1(SeqVarPVS1Helper):
    """Handles the PVS1 criteria assessment for sequence variants."""

    def __init__(self):
        super().__init__()
        self.prediction: PVS1Prediction = PVS1Prediction.NotPVS1
        self.prediction_path: PVS1PredictionSeqVarPath = PVS1PredictionSeqVarPath.NotSet

    def _convert_consequence(self, var_data: AutoACMGData) -> SeqVarPVS1Consequence:
        """
        Convert the VEP consequence of the sequence variant to the internal representation.

        Args:
            seqvar: The sequence variant being analyzed.

        Returns:
            SeqVarConsequence: The internal representation of the VEP consequence.
        """
        for conseq in var_data.consequence.mehari:
            if conseq in SeqvarConsequenceMapping:
                return SeqvarConsequenceMapping[conseq]
        if var_data.consequence.cadd in SeqvarConsequenceMapping:
            return SeqvarConsequenceMapping[var_data.consequence.cadd]
        return SeqVarPVS1Consequence.NotSet

    def verify_pvs1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[PVS1Prediction, PVS1PredictionSeqVarPath, str]:
        """Make the PVS1 prediction.

        The prediction is based on the PVS1 criteria for sequence variants. The prediction
        and prediction path is stored in the prediction and prediction_path attributes.

        Returns:
            Tuple[PVS1Prediction, PVS1PredictionSeqVarPath, str]: The prediction, prediction path,
                and the comment.
        """
        cons = self._convert_consequence(var_data)
        if cons == SeqVarPVS1Consequence.NonsenseFrameshift:
            self.comment_pvs1 = "Analysing as nonsense or frameshift variant. =>\n"
            if var_data.hgnc_id == "HGNC:9588":  # Follow guidelines for PTEN
                if var_data.prot_pos < 374:
                    self.comment_pvs1 = (
                        f"The variant is related to PTEN gene (hgnc: {var_data.hgnc_id})."
                        "The alteration results in a premature stop codon at protein "
                        f"position: {var_data.prot_pos} which is less than 374."
                    )
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.PTEN
                    return self.prediction, self.prediction_path, self.comment_pvs1

            if self.undergo_nmd(
                var_data.tx_pos_utr, var_data.hgnc_id, var_data.strand, var_data.exons
            ):
                self.comment_pvs1 += " =>\n"
                if self.in_bio_relevant_tsx(var_data.transcript_tags):
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.NF1
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.NF2
            else:
                self.comment_pvs1 += " =>\n"
                if self.crit4prot_func(seqvar, var_data.exons, var_data.strand):
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionSeqVarPath.NF3
                else:
                    self.comment_pvs1 += " =>\n"
                    if self.lof_freq_in_pop(
                        seqvar, var_data.exons, var_data.strand
                    ) or not self.in_bio_relevant_tsx(var_data.transcript_tags):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionSeqVarPath.NF4
                    else:
                        self.comment_pvs1 += " =>\n"
                        if self.lof_rm_gt_10pct_of_prot(var_data.prot_pos, var_data.prot_length):
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionSeqVarPath.NF5
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionSeqVarPath.NF6

        elif cons == SeqVarPVS1Consequence.SpliceSites:
            self.comment_pvs1 = "Analysing as splice site variant. =>\n"
            if self.exon_skip_or_cryptic_ss_disrupt(
                seqvar, var_data.exons, var_data.consequence.mehari, var_data.strand
            ) and self.undergo_nmd(
                var_data.tx_pos_utr, var_data.hgnc_id, var_data.strand, var_data.exons
            ):
                self.comment_pvs1 += " =>\n"
                if self.in_bio_relevant_tsx(var_data.transcript_tags):
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.SS1
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.SS2
            elif self.exon_skip_or_cryptic_ss_disrupt(
                seqvar, var_data.exons, var_data.consequence.mehari, var_data.strand
            ) and not self.undergo_nmd(
                var_data.tx_pos_utr, var_data.hgnc_id, var_data.strand, var_data.exons
            ):
                self.comment_pvs1 += " =>\n"
                if self.crit4prot_func(seqvar, var_data.exons, var_data.strand):
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionSeqVarPath.SS3
                else:
                    self.comment_pvs1 += " =>\n"
                    if self.lof_freq_in_pop(
                        seqvar, var_data.exons, var_data.strand
                    ) or not self.in_bio_relevant_tsx(var_data.transcript_tags):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionSeqVarPath.SS4
                    else:
                        self.comment_pvs1 += " =>\n"
                        if self.lof_rm_gt_10pct_of_prot(var_data.prot_pos, var_data.prot_length):
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionSeqVarPath.SS5
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionSeqVarPath.SS6
            else:
                self.comment_pvs1 += " =>\n"
                if self.crit4prot_func(seqvar, var_data.exons, var_data.strand):
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionSeqVarPath.SS10
                else:
                    self.comment_pvs1 += " =>\n"
                    if self.lof_freq_in_pop(
                        seqvar, var_data.exons, var_data.strand
                    ) or not self.in_bio_relevant_tsx(var_data.transcript_tags):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionSeqVarPath.SS7
                    else:
                        self.comment_pvs1 += " =>\n"
                        if self.lof_rm_gt_10pct_of_prot(var_data.prot_pos, var_data.prot_length):
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionSeqVarPath.SS8
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionSeqVarPath.SS9

        elif cons == SeqVarPVS1Consequence.InitiationCodon:
            self.comment_pvs1 = "Analysing as initiation codon variant. =>\n"
            if self.alt_start_cdn(var_data.cds_info, var_data.transcript_id):
                self.prediction = PVS1Prediction.NotPVS1
                self.prediction_path = PVS1PredictionSeqVarPath.IC3
            else:
                self.comment_pvs1 += " =>\n"
                if self.up_pathogenic_vars(
                    seqvar,
                    var_data.exons,
                    var_data.strand,
                    var_data.cds_info,
                    var_data.transcript_id,
                ):
                    self.prediction = PVS1Prediction.PVS1_Moderate
                    self.prediction_path = PVS1PredictionSeqVarPath.IC1
                else:
                    self.prediction = PVS1Prediction.PVS1_Supporting
                    self.prediction_path = PVS1PredictionSeqVarPath.IC2

        else:
            self.prediction = PVS1Prediction.UnsupportedConsequence
            self.prediction_path = PVS1PredictionSeqVarPath.NotSet
            logger.info(
                "PVS1 criteria is not met, since the variant consequence is {}.",
                cons.name,
            )
            self.comment_pvs1 = (
                f"Variant consequence is {cons.name}. " "PVS1 criteria cannot be applied."
            )

        return self.prediction, self.prediction_path, self.comment_pvs1

    def predict_pvs1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Predict the PVS1 criteria for the sequence variant.

        Args:
            seqvar: The sequence variant being analyzed.

        Returns:
            Tuple[PVS1Prediction, PVS1PredictionSeqVarPath, str]: The prediction, prediction path,
                and the comment.
        """
        pred, path, comment = self.verify_pvs1(seqvar, var_data)
        evidence_strength_mapping: Dict[PVS1Prediction, AutoACMGStrength] = {
            PVS1Prediction.PVS1: AutoACMGStrength.PathogenicVeryStrong,
            PVS1Prediction.PVS1_Strong: AutoACMGStrength.PathogenicStrong,
            PVS1Prediction.PVS1_Moderate: AutoACMGStrength.PathogenicModerate,
            PVS1Prediction.PVS1_Supporting: AutoACMGStrength.PathogenicSupporting,
            PVS1Prediction.NotPVS1: AutoACMGStrength.PathogenicVeryStrong,
        }
        return AutoACMGCriteria(
            name="PVS1",
            prediction=(
                AutoACMGPrediction.Met
                if pred
                in [
                    PVS1Prediction.PVS1,
                    PVS1Prediction.PVS1_Strong,
                    PVS1Prediction.PVS1_Moderate,
                    PVS1Prediction.PVS1_Supporting,
                ]
                else AutoACMGPrediction.NotMet
            ),
            strength=evidence_strength_mapping[pred],
            summary=comment,
            description=f"Prediction path: {PVS1PredictionPathMapping[path]}",
        )
