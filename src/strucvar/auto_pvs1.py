"""PVS1 criteria for Structural Variants (StrucVar)."""

from typing import Dict, List, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGStrength,
    AutoACMGStrucVarData,
    GenomicStrand,
)
from src.defs.auto_pvs1 import PVS1Prediction, PVS1PredictionPathMapping, PVS1PredictionStrucVarPath
from src.defs.exceptions import AlgorithmError, MissingDataError
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
        Check if the variant is a minimal deletion.

        Args:
            strucvar: The structural variant.
            exons: The exons of the gene.

        Returns:
            True if the deletion affects at least one full exon, False otherwise.

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

    def full_gene_del(self, strucvar: StrucVar, exons: List[Exon]) -> bool:
        """
        Check if the variant is a full gene deletion.

        Args:
            strucvar: The structural variant.
            exons: The exons of the gene.

        Returns:
            True if the variant is a full gene deletion, False otherwise.

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
            if (
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
        Check if the deletion is in a biologically relevant transcript.

        Check if the transcript has a MANE Select tag.

        Args:
            transcript_tags: The tags of the transcript.

        Returns:
            True if the deletion is in a biologically relevant transcript, False otherwise.
        """
        logger.info("Checking if the deletion is in a biologically relevant transcript.")
        self.comment_pvs1 += f"Transcript tags: {', '.join(transcript_tags)}."
        return "ManeSelect" in transcript_tags

    @staticmethod
    def crit4prot_func() -> bool:
        """Check if the deletion is critical for protein function."""
        return False

    @staticmethod
    def lof_freq_in_pop() -> bool:
        """Check if loss-of-function is frequent in the population."""
        return False

    @staticmethod
    def lof_rm_gt_10pct_of_prot() -> bool:
        """Check if the loss-of-function removes more than 10% of the protein."""
        return False

    @staticmethod
    def dup_disrupt_rf() -> bool:
        """Check if the duplication disrupts the reading frame."""
        return False

    @staticmethod
    def proven_in_tandem() -> bool:
        """Check if the duplication is proven in tandem."""
        return False

    @staticmethod
    def presumed_in_tandem() -> bool:
        """Check if the duplication is presumed in tandem."""
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

        The prediction is based on the PVS1 criteria for structural variants.
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
                if self.crit4prot_func():
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL4
                else:
                    self.comment_pvs1 += " =>"
                    if self.lof_freq_in_pop() or not self.in_bio_relevant_tsx(
                        var_data.transcript_tags
                    ):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionStrucVarPath.DEL5_1
                    else:
                        self.comment_pvs1 += " =>"
                        if self.lof_rm_gt_10pct_of_prot():
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL6_1
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL7_1
            else:
                self.comment_pvs1 += " =>"
                if self.crit4prot_func():
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionStrucVarPath.DEL8
                else:
                    self.comment_pvs1 += " =>"
                    if self.lof_freq_in_pop() or not self.in_bio_relevant_tsx(
                        var_data.transcript_tags
                    ):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionStrucVarPath.DEL5_2
                    else:
                        self.comment_pvs1 += " =>"
                        if self.lof_rm_gt_10pct_of_prot():
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL6_2
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionStrucVarPath.DEL7_2

        elif strucvar.sv_type == StrucVarType.DUP:
            self.comment_pvs1 = "Analysing the duplication variant. => "
            if self.proven_in_tandem():
                self.comment_pvs1 += " =>"
                if self.dup_disrupt_rf() and self.undergo_nmd(
                    strucvar, var_data.exons, var_data.strand
                ):
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DUP1
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DUP2_1
            elif self.presumed_in_tandem():
                self.comment_pvs1 += " =>"
                if self.dup_disrupt_rf() and self.undergo_nmd(
                    strucvar, var_data.exons, var_data.strand
                ):
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionStrucVarPath.DUP3
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionStrucVarPath.DUP2_2

            else:
                self.prediction = PVS1Prediction.NotPVS1
                self.prediction_path = PVS1PredictionStrucVarPath.DUP4

        else:
            self.prediction = PVS1Prediction.NotSet
            self.prediction_path = PVS1PredictionStrucVarPath.NotSet
            logger.error("Unsupported structural variant type: {}", strucvar.sv_type)
            self.comment_pvs1 = "Unsupported structural variant type."

        return self.prediction, self.prediction_path, self.comment_pvs1

    def predict_pvs1(self, strucvar: StrucVar, var_data: AutoACMGStrucVarData) -> AutoACMGCriteria:
        """Predict the PVS1 criteria for structural variants."""
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
