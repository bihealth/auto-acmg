"""PVS1 criteria for Structural Variants (StrucVar)."""

from typing import Dict, List, Tuple

from loguru import logger

from src.core.config import settings
from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGStrength,
    AutoACMGStrucVarData,
    GenomicStrand,
)
from src.defs.auto_pvs1 import PVS1Prediction, PVS1PredictionPathMapping, PVS1PredictionStrucVarPath
from src.defs.exceptions import (
    AlgorithmError,
    AutoAcmgBaseException,
    InvalidAPIResposeError,
    MissingDataError,
)
from src.defs.mehari import Exon
from src.defs.strucvar import StrucVar, StrucVarType
from src.utils import AutoACMGHelper


class StrucVarHelper(AutoACMGHelper):
    """Helper methods for PVS1 criteria for Structural Variants (StrucVar)."""

    def __init__(self):
        super().__init__()
        #: Comment to store the prediction explanation.
        self.comment_pvs1: str = ""

    def _minimal_deletion(self, strucvar: StrucVar, exons: List[Exon]) -> bool:
        """
        Check if the variant is a minimal deletion. A minimal deletion affects at least one full
        exon.

        Args:
            strucvar: The structural variant.
            exons: The exons of the gene.

        Returns:
            bool: True if the deletion affects at least one full exon, False otherwise.

        Raises:
            AlgorithmError: If the variant is not a deletion.
            MissingDataError: If exons are not available.
        """
        if strucvar.sv_type != StrucVarType.DEL:
            raise AlgorithmError("Variant is not a deletion.")
        if not exons:
            raise MissingDataError(
                "Exons are not available. Cannot determine if the variant is a minimal deletion."
            )
        for exon in exons:
            if strucvar.start <= exon.altStartI and strucvar.stop >= exon.altEndI:
                return True
        return False

    def _count_pathogenic_vars(self, strucvar: StrucVar) -> Tuple[int, int]:
        """
        Counts pathogenic variants in the range specified by the structural variant.

        The method retrieves variants from the range defined by the structural variant's start and
        stop positions and iterates through the ClinVar data of each variant to count the number of
        pathogenic variants and the total number of variants. The method considers a variant
        pathogenic if its classification is "Pathogenic" or "Likely pathogenic".

        Args:
            strucvar: The structural variant being analyzed.

        Returns:
            Tuple[int, int]: The number of pathogenic variants and the total number of variants.

        Raises:
            InvalidAPIResposeError: If the API response is invalid or cannot be processed.
        """
        logger.debug(
            "Counting pathogenic variants from position {} to {}.",
            strucvar.start,
            strucvar.stop,
        )
        if strucvar.stop < strucvar.start:
            raise AlgorithmError("End position is less than the start position.")

        response = self.annonars_client.get_variant_from_range(
            strucvar, strucvar.start, strucvar.stop
        )
        if response and response.clinvar:
            pathogenic_variants = [
                v
                for v in response.clinvar
                if v.records
                and (clf := v.records[0].classifications)
                and (gc := clf.germlineClassification)
                and gc.description in ["Pathogenic", "Likely pathogenic"]
            ]
            logger.debug(
                "Pathogenic variants: {}, Total variants: {}",
                len(pathogenic_variants),
                len(response.clinvar),
            )
            return len(pathogenic_variants), len(response.clinvar)
        else:
            raise InvalidAPIResposeError("Failed to get variant from range. No ClinVar data.")

    def _count_lof_vars(self, strucvar: StrucVar) -> Tuple[int, int]:
        """
        Counts Loss-of-Function (LoF) variants within the range of a structural variant.

        The method retrieves variants from the range defined by the structural variant's start and
        stop positions and iterates through the available data of each variant to count the number
        of LoF variants and the number of frequent LoF variants, based on the gnomAD genomes data (
        for the consequences of Nonsense and Frameshift variants) and the allele frequency (for the
        frequency of the LoF variants in the general population).

        Note:
            A LoF variant is considered frequent if its occurrence in the general population exceeds
            a threshold of 0.1%.

        Args:
            strucvar: The structural variant being analyzed.

        Returns:
            Tuple[int, int]: The number of frequent LoF variants and the total number of LoF
            variants.

        Raises:
            AlgorithmError: If the end position is less than the start position.
            InvalidAPIResposeError: If the API response is invalid or cannot be processed.
        """
        logger.debug(
            "Counting LoF variants from position {} to {}.",
            strucvar.start,
            strucvar.stop,
        )
        if strucvar.stop < strucvar.start:
            raise AlgorithmError("End position is less than the start position.")

        response = self.annonars_client.get_variant_from_range(
            strucvar, strucvar.start, strucvar.stop
        )
        if response and response.gnomad_genomes:
            frequent_lof_variants = 0
            lof_variants = 0
            for variant in response.gnomad_genomes:
                if not variant.vep:
                    continue
                for vep in variant.vep:
                    if vep.consequence in ["Nonsense", "Frameshift"]:
                        lof_variants += 1
                        if variant.alleleCounts:
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
            raise InvalidAPIResposeError(
                "Failed to get variant from range. No gnomAD genomes data."
            )

    def _calc_cds(
        self,
        exons: List[Exon],
        strand: GenomicStrand,
        start_codon: int,
        stop_codon: int,
    ) -> List[Exon]:
        """
        Remove UTRs from exons.

        Args:
            exons: List of exons for the gene.
            strand: The genomic strand of the gene.
            start_codon: Position of the start codon.
            stop_codon: Position of the stop codon.

        Returns:
            List[Exon]: List of exons without UTRs.

        Raises:
            MissingDataError: If the genomic strand is not set.
        """
        if strand == GenomicStrand.NotSet:
            raise MissingDataError("Genomic strand is not set. Cannot remove UTRs.")

        if strand == GenomicStrand.Plus:
            # Remove 5' UTR
            exons_remove = []
            for i, exon in enumerate(exons):
                if exon.altEndI - exon.altStartI + 1 > start_codon:
                    exon.altStartI += start_codon
                    break
                else:
                    exons_remove.append(i)
                    start_codon -= exon.altEndI - exon.altStartI + 1

            # Remove exons
            for i in exons_remove[::-1]:
                exons.pop(i)

            exons_remove = []
            # Remove 3' UTR
            for i, exon in enumerate(exons[::-1]):
                if exon.altEndI - exon.altStartI + 1 > stop_codon:
                    exon.altEndI -= stop_codon
                    break
                else:
                    exons_remove.append(i + 1)
                    stop_codon -= exon.altEndI - exon.altStartI + 1

            # Remove exons
            for i in exons_remove[::-1]:
                exons.pop(-i)

        elif strand == GenomicStrand.Minus:
            exons_remove = []
            # Remove 3' UTR
            for i, exon in enumerate(exons):
                if exon.altEndI - exon.altStartI + 1 > stop_codon:
                    exon.altStartI += stop_codon
                    break
                else:
                    exons_remove.append(i)
                    stop_codon -= exon.altEndI - exon.altStartI + 1

            # Remove exons
            for i in exons_remove[::-1]:
                exons.pop(i)

            exons_remove = []
            # Remove 5' UTR
            for i, exon in enumerate(exons[::-1]):
                if exon.altEndI - exon.altStartI + 1 > start_codon:
                    exon.altEndI -= start_codon
                    break
                else:
                    exons_remove.append(i + 1)
                    start_codon -= exon.altEndI - exon.altStartI + 1

            # Remove exons
            for i in exons_remove[::-1]:
                exons.pop(-i)

        return exons

    def full_gene_del(self, strucvar: StrucVar, exons: List[Exon]) -> bool:
        """
        Check if the variant is a full gene deletion. The deletion affects the whole gene if the
        start position of the deletion is less than or equal to the start of the first exon and the
        stop position of the deletion is greater than or equal to the end of the last exon.

        Args:
            strucvar: The structural variant.
            exons: The exons of the gene.

        Returns:
            bool: True if the variant is a full gene deletion, False otherwise.

        Raises:
            MissingDataError: If exons are not available.
        """
        logger.info("Checking if the variant is a full gene deletion.")
        if not exons:
            raise MissingDataError(
                "Exons are not available. Cannot determine if the variant is a full gene deletion."
            )

        gene_start = exons[0].altStartI
        gene_end = exons[-1].altEndI
        self.comment_pvs1 += f"Gene start: {gene_start}, gene end: {gene_end}."
        return strucvar.start <= gene_start and strucvar.stop >= gene_end

    def del_disrupt_rf(self, strucvar: StrucVar, exons: List[Exon], strand: GenomicStrand) -> bool:
        """
        Check if the single or multiple exon deletion disrupts the reading frame.

        Find the start and end positions of alteration based on the affected exon(s). If the
        positions lie within the intron(s) of the affected exon(s), the deletion does not disrupt
        the reading frame. Otherwise, there're two cases:
        - Check if the deletion starts within an exon. If so, check if the offset from the start
        of the exon to the start of the deletion is a multiple of 3. If so, the deletion does not
        disrupt the reading frame.
        - Check if the deletion stops within an exon. If so, check if the offset from the start
        of the last affected exon to the stop of the deletion is a multiple of 3. If so, the
        deletion does not disrupt the reading frame.

        Args:
            strucvar: The structural variant.
            exons: The exons of the gene.
            strand: The genomic strand of the variant.

        Returns:
            True if the deletion disrupts the reading frame, False otherwise.

        Raises:
            MissingDataError: If exons or strand are not available.
            AlgorithmError: Less than 1 full exon affected.
        """
        logger.info("Checking if the deletion disrupts the reading frame.")
        if not exons or strand == GenomicStrand.NotSet:
            raise MissingDataError(
                "Exons or strand are not available. Cannot determine if the deletion disrupts the "
                "reading frame."
            )

        # Find affected exons
        affected_exons = [
            exon
            for exon in exons
            if
            (
                # The deletion affects the whole exon
                (strucvar.start <= exon.altStartI and strucvar.stop >= exon.altEndI)
                # The deletion starts within the exon
                or (strucvar.start > exon.altStartI and strucvar.start < exon.altEndI)
                # The deletion ends within the exon
                or (strucvar.stop > exon.altStartI and strucvar.stop < exon.altEndI)
            )
        ]

        if (
            # No affected exons
            len(affected_exons) < 1
            # Deletion doesn't affect the full exon
            or (
                strucvar.start > affected_exons[0].altStartI
                and strucvar.stop < affected_exons[0].altEndI
            )
        ):
            raise AlgorithmError("The deletion affects less than one full exon.")

        first_affected_exon = affected_exons[0]
        last_affected_exon = affected_exons[-1]

        # Check if deletion is entirely within introns
        if (
            strucvar.start <= first_affected_exon.altStartI
            and strucvar.stop >= last_affected_exon.altEndI
        ):
            self.comment_pvs1 += "The deletion is entirely within introns. "
            return False

        if strand == GenomicStrand.Plus:
            # Case 1: Check if deletion starts within an exon
            if (
                strucvar.start > first_affected_exon.altStartI
                and strucvar.start <= first_affected_exon.altEndI
            ):
                self.comment_pvs1 += "The deletion starts within an exon. "
                return (strucvar.start - first_affected_exon.altStartI + 1) % 3 != 0

            # Case 2: Check if deletion stops within an exon
            if (
                strucvar.stop >= last_affected_exon.altStartI
                and strucvar.stop < last_affected_exon.altEndI
            ):
                self.comment_pvs1 += "The deletion stops within an exon. "
                return (strucvar.stop - last_affected_exon.altStartI + 1) % 3 != 0

        elif strand == GenomicStrand.Minus:
            # Case 1: Check if deletion starts within an exon
            if (
                strucvar.stop < last_affected_exon.altEndI
                and strucvar.stop >= last_affected_exon.altStartI
            ):
                self.comment_pvs1 += "The deletion starts within an exon. "
                return (last_affected_exon.altEndI - strucvar.stop + 1) % 3 != 0

            # Case 2: Check if deletion stops within an exon
            if (
                strucvar.start <= first_affected_exon.altEndI
                and strucvar.start > first_affected_exon.altStartI
            ):
                self.comment_pvs1 += "The deletion stops within an exon. "
                return (first_affected_exon.altEndI - strucvar.start + 1) % 3 != 0

        # If none of the above cases apply, the deletion doesn't disrupt the reading frame
        return False

    def undergo_nmd(self, strucvar: StrucVar, exons: List[Exon], strand: GenomicStrand) -> bool:
        """
        Check if the variant undergoes NMD.

        Check if the whole deletion affects only the last exon and 50 base pairs of the penultimate
        exon. If so, the variant does not undergo NMD.

        Args:
            strucvar: The structural variant.
            exons: The exons of the gene.
            strand: The genomic strand of the variant.

        Returns:
            True if the variant undergoes NMD, False otherwise.

        Raises:
            MissingDataError: If exons or strand are not available.
            AlgorithmError: If less than 2 exons are available.
        """
        logger.info("Checking if the variant undergoes NMD.")
        if not exons or strand == GenomicStrand.NotSet:
            raise MissingDataError(
                "Exons or strand are not available. Cannot determine if the variant undergoes NMD."
            )
        if len(exons) < 2:
            raise AlgorithmError(
                "Exons are not available. Cannot determine if the variant undergoes NMD."
            )

        if strand == GenomicStrand.Plus:
            nmd_cutoff = max(exons[-2].altEndI - 50, exons[-2].altStartI)
            self.comment_pvs1 += f"NMD cutoff: {nmd_cutoff}."
            if strucvar.start >= nmd_cutoff:
                self.comment_pvs1 += "The variant undergoes NMD."
                return False
        elif strand == GenomicStrand.Minus:
            nmd_cutoff = min(exons[1].altStartI + 50, exons[1].altEndI)
            self.comment_pvs1 += f"NMD cutoff: {nmd_cutoff}."
            if strucvar.stop <= nmd_cutoff:
                self.comment_pvs1 += "The variant undergoes NMD."
                return False

        self.comment_pvs1 += "The variant does not undergo NMD."
        return True

    def in_bio_relevant_tsx(self, transcript_tags: List[str]) -> bool:
        """
        Check if the deletion is in a biologically relevant transcript. Check if the transcript has
        a MANE Select tag.

        Args:
            transcript_tags: The tags of the transcript.

        Returns:
            bool: True if the deletion is in a biologically relevant transcript, False otherwise.
        """
        logger.info("Checking if the deletion is in a biologically relevant transcript.")
        self.comment_pvs1 += f"Transcript tags: {', '.join(transcript_tags)}."
        return "TRANSCRIPT_TAG_MANE_SELECT" in transcript_tags

    def crit4prot_func(self, strucvar: StrucVar) -> bool:
        """
        Check if the deletion is critical for protein function.

        This method is implemented by fetching variants from the start to the end of the structural
        variant, then counting the number of pathogenic variants in that region, by iterating
        through the clinvar data of each variant. Consider the region critical if the frequency of
        pathogenic variants exceeds 5%.

        Args:
            strucvar: The structural variant being analyzed.

        Returns:
            bool: True if the deletion is critical for protein function, False otherwise.

        Raises:
            AlgorithmError: If the API response is invalid or cannot be processed.
        """
        logger.debug("Checking if the deletion is critical for protein function.")

        try:
            pathogenic_variants, total_variants = self._count_pathogenic_vars(strucvar)
            self.comment_pvs1 += (
                f"Found {pathogenic_variants} pathogenic variants from {total_variants} total "
                f"variants in the range {strucvar.start} - {strucvar.stop}. "
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
            raise AlgorithmError("Failed to predict criticality for variant.") from e

    def lof_freq_in_pop(self, strucvar: StrucVar) -> bool:
        """
        Checks if the Loss-of-Function (LoF) variants within the structural variant are frequent in
        the general population.

        This function determines the frequency of LoF variants within the range specified by the
        structural variant and evaluates whether this frequency exceeds a defined threshold
        indicative of common occurrence in the general population.

        Implementation of the rule:
        - Retrieving the number of LoF variants and frequent LoF variants in the range defined by
        the structural variant.
        - Considering the LoF variants frequent in the general population if the frequency of
        "frequent" LoF variants exceeds 10%.

        Note:
        A LoF variant is considered frequent if its occurrence in the general population exceeds
        some threshold. We use a threshold of 10% to determine if the LoF variant is frequent.

        Args:
            strucvar: The structural variant being analyzed.

        Returns:
            bool: True if the LoF variant frequency is greater than 10%, False otherwise.

        Raises:
            AlgoritmError: If the API response is invalid or cannot be processed.
        """
        logger.debug(
            "Checking if LoF variants are frequent in the general population for "
            "structural variants."
        )

        try:
            frequent_lof_variants, lof_variants = self._count_lof_vars(strucvar)
            self.comment_pvs1 += (
                f"Found {frequent_lof_variants} frequent LoF variants from {lof_variants} total "
                f"LoF variants in the range {strucvar.start} - {strucvar.stop}. "
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
            raise AlgorithmError("Failed to predict LoF frequency for structural variant.") from e

    def lof_rm_gt_10pct_of_prot(
        self,
        strucvar: StrucVar,
        exons: List[Exon],
        strand: GenomicStrand,
        start_codon: int,
        stop_codon: int,
    ) -> bool:
        """
        Determine if the deletion removes more than 10% of the protein-coding sequence.

        First remove the UTRs from the exons. Then iterate through the CDS exons and calculate the
        total CDS length and the length of the deleted region. Return True if the deletion removes
        more than 10% of the protein-coding sequence, False otherwise.

        Args:
            strucvar: The structural variant being analyzed.
            exons: List of exons for the gene.
            strand: The genomic strand of the gene.
            start_codon: Position of the start codon.
            stop_codon: Position of the stop codon.

        Returns:
            bool: True if the deletion removes more than 10% of the protein, False otherwise.

        Raises:
            AlgorithmError: If the total CDS length is zero.
        """
        total_cds_length = 0
        deleted_length = 0

        cds = self._calc_cds(exons, strand, start_codon, stop_codon)

        for exon in cds:
            total_cds_length += exon.altEndI - exon.altStartI + 1

            overlap_start = max(strucvar.start, exon.altStartI)
            overlap_end = min(strucvar.stop, exon.altEndI)

            if overlap_start <= overlap_end:  # Check if there is an overlap
                overlap_length = overlap_end - overlap_start + 1
                deleted_length += overlap_length

        if total_cds_length == 0:
            raise AlgorithmError(
                "Total CDS length is zero. Cannot determine if deletion removes "
                "more than 10% of the protein."
            )
        deletion_percentage = (deleted_length // 3) / (total_cds_length // 3)
        return deletion_percentage > 0.1

    @staticmethod
    def dup_disrupt_rf() -> bool:
        """
        Check if the duplication disrupts the reading frame.
        NOT IMPLEMENTED!
        """
        return False

    @staticmethod
    def proven_in_tandem() -> bool:
        """
        Check if the duplication is proven in tandem.
        NOT IMPLEMENTED!
        """
        return False

    @staticmethod
    def presumed_in_tandem() -> bool:
        """
        Check if the duplication is presumed in tandem.
        NOT IMPLEMENTED!
        """
        return False


class AutoPVS1(StrucVarHelper):
    """Handles the PVS1 criteria assesment for structural variants."""

    def __init__(self):
        super().__init__()
        self.prediction: PVS1Prediction = PVS1Prediction.NotPVS1
        self.prediction_path: PVS1PredictionStrucVarPath = PVS1PredictionStrucVarPath.NotSet

    def verify_pvs1(  # pragma: no cover
        self, strucvar: StrucVar, var_data: AutoACMGStrucVarData
    ) -> Tuple[PVS1Prediction, PVS1PredictionStrucVarPath, str]:
        """Make the PVS1 prediction.

        The prediction is based on the PVS1 decision tree for structural variants.

        Args:
            strucvar: The structural variant.
            var_data: The variant information.

        Returns:
            Tuple[PVS1Prediction, PVS1PredictionStrucVarPath, str]: The prediction, prediction path,
            and the comment.
        """
        if strucvar.sv_type == StrucVarType.DEL:
            self.comment_pvs1 = "Analysing the deletion variant. => "
            if not self._minimal_deletion(strucvar, var_data.exons):
                return (
                    PVS1Prediction.NotPVS1,
                    PVS1PredictionStrucVarPath.NotSet,
                    "Variant is not a minimal deletion. Must affect at least one full exon.",
                )
            if self.full_gene_del(strucvar, var_data.exons):
                self.prediction = PVS1Prediction.PVS1
                self.prediction_path = PVS1PredictionStrucVarPath.DEL1
            elif self.del_disrupt_rf(
                strucvar, var_data.exons, var_data.strand
            ) and self.undergo_nmd(strucvar, var_data.exons, var_data.strand):
                self.comment_pvs1 += " =>"
                if self.in_bio_relevant_tsx(var_data.transcript_tags):
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL2
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL3
            elif self.del_disrupt_rf(
                strucvar, var_data.exons, var_data.strand
            ) and not self.undergo_nmd(strucvar, var_data.exons, var_data.strand):
                self.comment_pvs1 += " =>"
                if self.crit4prot_func(strucvar):
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL4
                else:
                    self.comment_pvs1 += " =>"
                    if self.lof_freq_in_pop(strucvar) or not self.in_bio_relevant_tsx(
                        var_data.transcript_tags
                    ):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionStrucVarPath.DEL5_1
                    else:
                        self.comment_pvs1 += " =>"
                        if self.lof_rm_gt_10pct_of_prot(
                            strucvar,
                            var_data.exons,
                            var_data.strand,
                            var_data.start_cdn,
                            var_data.stop_cdn,
                        ):
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL6_1
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL7_1
            else:
                self.comment_pvs1 += " =>"
                if self.crit4prot_func(strucvar):
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL8
                else:
                    self.comment_pvs1 += " =>"
                    if self.lof_freq_in_pop(strucvar) or not self.in_bio_relevant_tsx(
                        var_data.transcript_tags
                    ):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionStrucVarPath.DEL5_2
                    else:
                        self.comment_pvs1 += " =>"
                        if self.lof_rm_gt_10pct_of_prot(
                            strucvar,
                            var_data.exons,
                            var_data.strand,
                            var_data.start_cdn,
                            var_data.stop_cdn,
                        ):
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL6_2
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL7_2

        elif strucvar.sv_type == StrucVarType.DUP:
            self.comment_pvs1 = (
                "THE CRITERIA RELY ON THE DUPLICATION_TANDEM SETTING. "
                "Please, specify the entry in CLI or API. Per default, the criteria is set to "
                "False (Not PVS1). "
            )
            self.comment_pvs1 += "Analysing the duplication variant. => "

            if settings.DUPLICATION_TANDEM:
                self.comment_pvs1 += (
                    "The duplication is in tandem AND disrupts reading frame AND undergoes NMD. "
                )
                self.prediction = PVS1Prediction.PVS1
                self.prediction_path = PVS1PredictionStrucVarPath.DUP1
            else:
                self.comment_pvs1 += (
                    "The duplication is not in tandem OR does not disrupt reading frame OR "
                    "does not undergo NMD."
                )
                self.prediction = PVS1Prediction.NotPVS1
                self.prediction_path = PVS1PredictionStrucVarPath.DUP3

            # if self.proven_in_tandem():
            #     self.comment_pvs1 += " =>"
            #     if self.dup_disrupt_rf() and self.undergo_nmd(
            #         strucvar, var_data.exons, var_data.strand
            #     ):
            #         self.prediction = PVS1Prediction.PVS1
            #         self.prediction_path = PVS1PredictionStrucVarPath.DUP1
            #     else:
            #         self.prediction = PVS1Prediction.NotPVS1
            #         self.prediction_path = PVS1PredictionStrucVarPath.DUP2_1
            # elif self.presumed_in_tandem():
            #     self.comment_pvs1 += " =>"
            #     if self.dup_disrupt_rf() and self.undergo_nmd(
            #         strucvar, var_data.exons, var_data.strand
            #     ):
            #         self.prediction = PVS1Prediction.PVS1_Strong
            #         self.prediction_path = PVS1PredictionStrucVarPath.DUP3
            #     else:
            #         self.prediction = PVS1Prediction.NotPVS1
            #         self.prediction_path = PVS1PredictionStrucVarPath.DUP2_2

            # else:
            #     self.prediction = PVS1Prediction.NotPVS1
            #     self.prediction_path = PVS1PredictionStrucVarPath.DUP4

        else:
            self.prediction = PVS1Prediction.NotSet
            self.prediction_path = PVS1PredictionStrucVarPath.NotSet
            logger.error("Unsupported structural variant type: {}", strucvar.sv_type)
            self.comment_pvs1 = "Unsupported structural variant type."

        return self.prediction, self.prediction_path, self.comment_pvs1

    def predict_pvs1(self, strucvar: StrucVar, var_data: AutoACMGStrucVarData) -> AutoACMGCriteria:
        """Predict the PVS1 criteria."""
        pred, path, comment = self.verify_pvs1(strucvar, var_data)
        evidence_strength_mapping: Dict[PVS1Prediction, AutoACMGStrength] = {
            PVS1Prediction.PVS1: AutoACMGStrength.PathogenicVeryStrong,
            PVS1Prediction.PVS1_Strong: AutoACMGStrength.PathogenicStrong,
            PVS1Prediction.PVS1_Moderate: AutoACMGStrength.PathogenicModerate,
            PVS1Prediction.PVS1_Supporting: AutoACMGStrength.PathogenicSupporting,
            PVS1Prediction.NotPVS1: AutoACMGStrength.PathogenicVeryStrong,
            PVS1Prediction.UnsupportedConsequence: AutoACMGStrength.PathogenicVeryStrong,
            PVS1Prediction.NotSet: AutoACMGStrength.PathogenicVeryStrong,
        }
        return AutoACMGCriteria(
            name="PVS1",
            prediction=(
                AutoACMGPrediction.Applicable
                if pred
                in [
                    PVS1Prediction.PVS1,
                    PVS1Prediction.PVS1_Strong,
                    PVS1Prediction.PVS1_Moderate,
                    PVS1Prediction.PVS1_Supporting,
                ]
                else AutoACMGPrediction.NotApplicable
            ),
            strength=evidence_strength_mapping[pred],
            summary=comment,
            description=f"Prediction path: {PVS1PredictionPathMapping[path]}",
        )
