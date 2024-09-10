"""Implementation of BP7 criteria."""

from typing import Optional, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    BP7,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    GenomicStrand,
)
from src.defs.exceptions import AutoAcmgBaseException, MissingDataError
from src.defs.seqvar import SeqVar
from src.utils import AutoACMGHelper


class AutoBP7(AutoACMGHelper):
    """Class for BP7 prediction."""

    def __init__(self):
        super().__init__()
        #: Prediction result.
        self.prediction_bp7: Optional[BP7] = None
        #: Comment to store the prediction explanation.
        self.comment_bp7: str = ""

    @staticmethod
    def _spliceai_impact(var_data: AutoACMGSeqVarData) -> bool:
        """
        Predict splice site alterations using SpliceAI.

        If any of SpliceAI scores are greater than specific thresholds, the variant is considered a
        splice site alteration. The thresholds are defined in the variant data thresholds.

        Args:
            var_data: The data containing variant scores and thresholds.

        Returns:
            bool: True if the variant is a splice site alteration, False otherwise.
        """
        score_checks = {
            "spliceAI_acceptor_gain": var_data.thresholds.spliceAI_acceptor_gain,
            "spliceAI_acceptor_loss": var_data.thresholds.spliceAI_acceptor_loss,
            "spliceAI_donor_gain": var_data.thresholds.spliceAI_donor_gain,
            "spliceAI_donor_loss": var_data.thresholds.spliceAI_donor_loss,
        }
        return any(
            (getattr(var_data.scores.cadd, score_name) or 0) > threshold
            for score_name, threshold in score_checks.items()
        )

    def _is_conserved(self, var_data: AutoACMGSeqVarData) -> bool:
        """
        Predict if the variant is conserved.

        Check if the variant is conserved using the phyloP100 score.

        Args:
            variant_info: The variant information.

        Returns:
            bool: True if the variant is not conserved, False otherwise.
        """
        phylop = var_data.scores.cadd.phyloP100 or var_data.scores.dbnsfp.phyloP100
        if not phylop:
            raise MissingDataError("Missing phyloP100 score.")
        if phylop >= var_data.thresholds.phyloP100:
            return True
        return False

    @staticmethod
    def _is_synonymous(var_data: AutoACMGSeqVarData) -> bool:
        """
        Predict if the variant is synonymous.

        Check if the variant's consequence is synonymous.

        Args:
            variant_info: The variant information.

        Returns:
            bool: True if the variant is synonymous, False otherwise.
        """
        if (
            "synonymous_variant" in var_data.consequence.mehari
            or var_data.consequence.cadd == "synonymous"
        ):
            return True
        return False

    @staticmethod
    def _is_intronic(var_data: AutoACMGSeqVarData) -> bool:
        """
        Predict if the variant is intronic.

        Check if the variant's consequence is intronic.

        Args:
            variant_info: The variant information.

        Returns:
            bool: True if the variant is intronic, False otherwise.
        """
        intronic_consequences = [
            "intron_variant",
            "intergenic_variant",
            "upstream_gene_variant",
            "downstream_gene_variant",
            "5_prime_utr_variant",
            "5_prime_UTR_variant",
            "3_prime_utr_variant",
            "3_prime_UTR_variant",
            "splice_region_variant",
            "splice_donor_5th_base_variant",
            "splice_donor_region_variant",
            "splice_polypyrimidine_tract_variant",
        ]
        if any(
            consequence in var_data.consequence.mehari for consequence in intronic_consequences
        ) or any(consequence in var_data.consequence.cadd for consequence in ["intron", "splice"]):
            return True
        return False

    def _affect_canonical_ss(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """
        Predict if the variant affects canonical splice site.

        Check if the variant position is within +1/-2 of the start/end of an intron.
        Note, that for the minus strand, the donor site is at the start of the intron and the
        acceptor site is at the end of the intron.

        Args:
            seqvar: The sequence variant.
            var_data: The variant data.

        Returns:
            bool: True if the variant affects canonical splice site, False otherwise.

        Raises:
            MissingDataError: If the strand information is missing.
        """
        for exon in var_data.exons:
            if var_data.strand == GenomicStrand.Plus:
                # Check the acceptor site
                if (
                    exon.altStartI - var_data.thresholds.bp7_acceptor
                    <= seqvar.pos
                    <= exon.altStartI
                ):
                    return True
                # Check the donor site
                if exon.altEndI <= seqvar.pos <= exon.altEndI + var_data.thresholds.bp7_donor:
                    return True
            elif var_data.strand == GenomicStrand.Minus:
                # Check the donor site
                if exon.altStartI - var_data.thresholds.bp7_donor <= seqvar.pos <= exon.altStartI:
                    return True
                # Check the acceptor site
                if exon.altEndI <= seqvar.pos <= exon.altEndI + var_data.thresholds.bp7_acceptor:
                    return True
            else:
                raise MissingDataError("Missing strand information.")
        return False

    def _is_bp7_exception(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """
        Help function to check if the variant is an exception.

        Per default there are no exceptions for BP7.

        Args:
            seqvar: The sequence variant.
            var_data: The variant data.

        Returns:
            bool: True if the variant is an exception, False otherwise.
        """
        return False

    def verify_bp7(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> Tuple[Optional[BP7], str]:
        """Predict BP7 criterion."""
        self.prediction_bp7 = BP7()
        self.comment_bp7 = ""
        try:
            if seqvar.chrom == "MT":
                self.comment_bp7 = "Variant is in the mitochondrial genome. BP7 is not met."
                self.prediction_bp7.BP7 = False
                return self.prediction_bp7, self.comment_bp7

            if (
                (
                    # Synonymous variant
                    self._is_synonymous(var_data)
                    # Intronic variant that does not affect canonical splice site (ClinGen modified)
                    or (
                        self._is_intronic(var_data)
                        and not self._affect_canonical_ss(seqvar, var_data)
                    )
                )
                and not self._is_conserved(var_data)
                and not self._spliceai_impact(var_data)
                and not self._is_bp7_exception(seqvar, var_data)
            ):
                self.comment_bp7 += (
                    "Synonymous variant is not conserved and not predicted to affect splicing. "
                    "Or intronic variant that does not affect canonical splice site. "
                    f"PhyloP100 score: {var_data.scores.cadd.phyloP100}. "
                    f"SpliceAI scores: {var_data.scores.cadd.spliceAI_acceptor_gain},"
                    f"{var_data.scores.cadd.spliceAI_acceptor_loss}, "
                    f"{var_data.scores.cadd.spliceAI_donor_gain}, "
                    f"{var_data.scores.cadd.spliceAI_donor_loss}. "
                    "BP7 is met."
                )
                self.prediction_bp7.BP7 = True
            else:
                self.comment_bp7 += (
                    "Variant is not synonymous, or is conserved, or predicted to affect splicing. "
                    "Or intronic variant that affects canonical splice site. "
                    f"PhyloP100 score: {var_data.scores.cadd.phyloP100}. "
                    f"SpliceAI scores: {var_data.scores.cadd.spliceAI_acceptor_gain},"
                    f"{var_data.scores.cadd.spliceAI_acceptor_loss}, "
                    f"{var_data.scores.cadd.spliceAI_donor_gain}, "
                    f"{var_data.scores.cadd.spliceAI_donor_loss}. "
                    "BP7 is not met."
                )
                self.prediction_bp7.BP7 = False

        except AutoAcmgBaseException as e:
            self.comment_bp7 = f"Failed to predict BP7 criterion. Error: {e}"
            self.prediction_bp7 = None

        return self.prediction_bp7, self.comment_bp7

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Predict BP7 criterion."""
        logger.info("Predict BP7")
        pred, comment = self.verify_bp7(seqvar, var_data)
        if pred:
            pred_bp7 = (
                AutoACMGPrediction.Applicable
                if pred.BP7
                else (
                    AutoACMGPrediction.NotApplicable
                    if pred.BP7 is False
                    else AutoACMGPrediction.Failed
                )
            )
            strength_bp7 = pred.BP7_strength
        else:
            pred_bp7 = AutoACMGPrediction.Failed
            strength_bp7 = AutoACMGStrength.BenignSupporting
        return AutoACMGCriteria(
            name="BP7",
            prediction=pred_bp7,
            strength=strength_bp7,
            summary=comment,
        )
