"""Implementation of PS1 and PM5 prediction for sequence variants."""

import re
from typing import Optional, Tuple

from loguru import logger

from src.defs.annonars_variant import AnnonarsVariantResponse, VariantResult
from src.defs.auto_acmg import (
    PS1PM5,
    AminoAcid,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
)
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException
from src.defs.seqvar import SeqVar
from src.utils import AutoACMGHelper

#: DNA bases
DNA_BASES = ["A", "C", "G", "T"]

#: Regular expression for parsing pHGVSp
REGEX_HGVSP = re.compile(r"p\.(\D+)(\d+)(\D+)")


class AutoPS1PM5(AutoACMGHelper):
    """Class for PS1 and PM5 prediction."""

    def __init__(self):
        super().__init__()
        #: Prediction result.
        self.prediction_ps1pm5: Optional[PS1PM5] = None
        #: Comment to store the prediction explanation.
        self.comment_ps1pm5: str = ""

    def _get_var_info(self, seqvar: SeqVar) -> Optional[AnnonarsVariantResponse]:
        """Get variant information from Annonars.

        Args:
            seqvar: The sequence variant.

        Returns:
            Optional[AnnonarsVariantResponse]: Annonars response if the variant is found, None
            otherwise.
        """
        try:
            return self.annonars_client.get_variant_info(seqvar)
        except AutoAcmgBaseException:
            return None

    def _parse_HGVSp(self, pHGVSp: str) -> Optional[AminoAcid]:
        """Parse the pHGVSp from VEP into its components.

        Args:
            pHGVSp: The protein change in HGVS format.

        Returns:
            Optional[AminoAcid]: The amino acid change if the pHGVSp is valid, None otherwise.
        """
        try:
            # If multiple pHGVSp values are present, take the first one
            if ";" in pHGVSp:
                pHGVSp = pHGVSp.split(";")[0]
            match = re.match(REGEX_HGVSP, pHGVSp)
            if match:
                return AminoAcid[match.group(3)]
            else:
                return None
        except AutoAcmgBaseException:
            logger.debug("Invalid pHGVSp: {}", pHGVSp)
            return None

    def _is_pathogenic(self, variant_info: VariantResult) -> bool:
        """Check if the variant is pathogenic based on ClinVar data.

        Args:
            variant_info: Annonars variant information

        Returns:
            bool: True if the variant is pathogenic
        """
        if variant_info.clinvar and variant_info.clinvar.records:
            for rec in variant_info.clinvar.records:
                if (
                    (r := rec.classifications)
                    and (g := r.germlineClassification)
                    and g.description in ["Pathogenic", "Likely pathogenic"]
                ):
                    return True
        return False

    def _is_missense(self, var_data: AutoACMGSeqVarData) -> bool:
        """
        Check if the variant's consequence is missense.

        Args:
            var_data: The variant information.

        Returns:
            bool: True if the variant is a missense variant, False otherwise.
        """
        if "missense" in var_data.consequence.cadd:
            return True
        if any("missense" in cons for cons in var_data.consequence.mehari):
            return True
        return False

    def _is_splice_affecting(self, var_data: AutoACMGSeqVarData) -> bool:
        """
        Check if the variant is a splice-affecting variant.

        Args:
            var_data: The variant information.

        Returns:
            bool: True if the variant is a splice-affecting variant, False otherwise.
        """
        if "splice" in var_data.consequence.cadd:
            return True
        if any("splice" in cons for cons in var_data.consequence.mehari):
            return True
        return False

    def _affect_splicing(self, var_data: AutoACMGSeqVarData) -> bool:
        """
        Check if the variant affects splicing. If any of spliceAI scores are above the threshold,
        the variant is considered to affect splicing.

        Args:
            var_data: The variant information.

        Returns:
            bool: True if the variant affects splicing, False otherwise.
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

    def verify_ps1pm5(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PS1PM5], str]:
        """
        Predicts the criteria PS1 and PM5 for the provided sequence variant.

        Implementation of the rule PS1 and PM5:
            The method implements the rule by:
            - Getting the primary variant information & parsing the primary amino acid change.
            - Iterating over all possible alternative bases & getting the alternative variant
            information.
            - Parsing the alternative amino acid change & checking if the alternative variant is
            pathogenic.
            - If the alternative variant is pathogenic and the amino acid change is the same as the
            primary variant, then PS1 is set to True.
            - If the alternative variant is pathogenic and the amino acid change is different from
            the primary variant, then PM5 is set to True.

        Note:
            Rules:
            PS1: Same amino acid change as a previously established pathogenic variant regardless of
            nucleotide change.
            PM5: Novel missense change at an amino acid residue where a different missense change
            determined to be pathogenic has been seen before.

        Returns:
            Tuple[Optional[PS1PM5], str]: Prediction result and the comment with the explanation.
        """
        self.prediction_ps1pm5 = PS1PM5()
        self.comment_ps1pm5 = ""
        try:
            if not self._is_missense(var_data) and not self._is_splice_affecting(var_data):
                raise AlgorithmError(
                    "Variant is not a missense or splice-affecting variant. PS1/PM5 not applicable."
                )

            primary_aa_change = self._parse_HGVSp(var_data.pHGVS)
            if not primary_aa_change:
                raise AlgorithmError("No valid primary amino acid change for PS1/PM5 prediction.")

            self.comment_ps1pm5 = f"Primary amino acid change: {primary_aa_change.name}. "
            for alt_base in DNA_BASES:
                # Skip the same base insert
                if alt_base == seqvar.insert:
                    continue

                self.comment_ps1pm5 += f"Analysing alternative variant with base: {alt_base}. "
                alt_seqvar = SeqVar(
                    genome_release=seqvar.genome_release,
                    chrom=seqvar.chrom,
                    pos=seqvar.pos,
                    delete=seqvar.delete,
                    insert=alt_base,
                )
                alt_info = self._get_var_info(alt_seqvar)

                if alt_info and alt_info.result.dbnsfp and alt_info.result.dbnsfp.HGVSp_VEP:
                    alt_aa_change = self._parse_HGVSp(alt_info.result.dbnsfp.HGVSp_VEP)
                    self.comment_ps1pm5 += (
                        "Alternative amino acid change: "
                        f"{alt_aa_change.name if alt_aa_change else "N/A"}. "
                    )
                    if alt_aa_change and self._is_pathogenic(alt_info.result):
                        if primary_aa_change == alt_aa_change:
                            self.comment_ps1pm5 += (
                                "Alternative variant is pathogenic "
                                "and amino acid change is the same as primary variant."
                                "Result: PS1 is met."
                            )
                            self.prediction_ps1pm5.PS1 = True
                        if primary_aa_change != alt_aa_change:
                            self.comment_ps1pm5 += (
                                "Alternative variant is pathogenic "
                                "and amino acid change is different from primary variant."
                                "Result: PM5 is met."
                            )
                            self.prediction_ps1pm5.PM5 = True
                else:
                    self.comment_ps1pm5 += "Missing dbnsfp data."
                    self.comment_ps1pm5 += f"Failed to get variant information for {alt_seqvar}. "

        except AutoAcmgBaseException as e:
            self.comment_ps1pm5 = f"Error occurred during PS1/PM5 prediction. Error: {e}"
            self.prediction_ps1pm5 = None
        return self.prediction_ps1pm5, self.comment_ps1pm5

    def predict_ps1pm5(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Predict PS1 and PM5 criteria."""
        logger.info("Predict PS1 and PM5")
        pred, comment = self.verify_ps1pm5(seqvar, var_data)
        if pred:
            ps1_pred = (
                AutoACMGPrediction.Applicable
                if pred.PS1
                else (
                    AutoACMGPrediction.NotApplicable
                    if pred.PS1 is False
                    else AutoACMGPrediction.Failed
                )
            )
            pm5_pred = (
                AutoACMGPrediction.Applicable
                if pred.PM5
                else (
                    AutoACMGPrediction.NotApplicable
                    if pred.PM5 is False
                    else AutoACMGPrediction.Failed
                )
            )
            ps1_strength = pred.PS1_strength
            pm5_strength = pred.PM5_strength
        else:
            ps1_pred = AutoACMGPrediction.Failed
            pm5_pred = AutoACMGPrediction.Failed
            ps1_strength = AutoACMGStrength.PathogenicStrong
            pm5_strength = AutoACMGStrength.PathogenicModerate
        return (
            AutoACMGCriteria(
                name="PS1",
                prediction=ps1_pred,
                strength=ps1_strength,
                summary=comment,
            ),
            AutoACMGCriteria(
                name="PM5",
                prediction=pm5_pred,
                strength=pm5_strength,
                summary=comment,
            ),
        )
