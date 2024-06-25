"""Implementation of PS1 and PM5 prediction for sequence variants."""

import re
from typing import Optional, Tuple

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.annonars_variant import AnnonarsVariantResponse, VariantResult
from src.defs.auto_acmg import PS1PM5, AminoAcid
from src.defs.auto_pvs1 import SeqVarConsequence
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.utils import SeqVarTranscriptsHelper

#: DNA bases
DNA_BASES = ["A", "C", "G", "T"]

#: Regular expression for parsing pHGVSp
REGEX_HGVSP = re.compile(r"p\.(\D+)(\d+)(\D+)")


class AutoPS1PM5:
    """Predicts PS1 and PM5 criteria for sequence variants."""

    def __init__(
        self,
        seqvar: SeqVar,
        variant_info: VariantResult,
        *,
        config: Optional[Config] = None,
    ):
        #: Configuration to use.
        self.config: Config = config or Config()
        #: Sequence variant to predict.
        self.seqvar: SeqVar = seqvar
        #: Variant information.
        self.variant_info: VariantResult = variant_info
        #: Annonars client.
        self.annonars_client: AnnonarsClient = AnnonarsClient(
            api_base_url=self.config.api_base_url_annonars
        )
        #: Prediction result.
        self.prediction: Optional[PS1PM5] = None
        #: Comment to store the prediction explanation.
        self.comment: str = ""

    def _get_variant_info(self, seqvar: SeqVar) -> Optional[AnnonarsVariantResponse]:
        """Get variant information from Annonars.

        Returns:
            AnnonarsVariantResponse: Annonars response.
        """
        try:
            logger.debug("Getting variant information for {}.", seqvar)
            return self.annonars_client.get_variant_info(seqvar)
        except AutoAcmgBaseException as e:
            logger.error("Failed to get variant information. Error: {}", e)
            return None

    def _parse_HGVSp(self, pHGVSp: str) -> Optional[AminoAcid]:
        """Parse the pHGVSp from VEP into its components.

        Args:
            pHGVSp (str): The protein change in HGVS format.

        Returns:
            AminoAcid: The new amino acid if the pHGVSp is valid, None otherwise.
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
        """Check if the variant is pathogenic.

        Args:
            variant_info (dict): Annonars variant information

        Returns:
            bool: True if the variant is pathogenic
        """
        if variant_info.clinvar and variant_info.clinvar.records:
            for rec in variant_info.clinvar.records:
                if (
                    rec.classifications
                    and rec.classifications.germlineClassification
                    and rec.classifications.germlineClassification.description
                    in [
                        "Pathogenic",
                        "Likely pathogenic",
                    ]
                ):
                    return True
        return False

    def _is_missense(self) -> bool:
        """
        Check if the variant is a missense variant.

        Fetches the transcript information for the sequence variant and checks if the variant
        consequence is missense.

        Returns:
            bool: True if the variant is a missense variant.
        """
        seqvar_transcript_helper = SeqVarTranscriptsHelper(self.seqvar, config=self.config)
        seqvar_transcript_helper.initialize()
        (_, _, _, _, consequence) = seqvar_transcript_helper.get_ts_info()
        return consequence == SeqVarConsequence.Missense

    def predict(self) -> Tuple[Optional[PS1PM5], str]:
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
        try:
            # Initialize the prediction result
            self.prediction = PS1PM5()

            if (
                not self.variant_info
                or not self.variant_info.dbnsfp
                or not self.variant_info.dbnsfp.HGVSp_VEP
            ):
                raise MissingDataError(
                    "No valid primary variant information for PS1/PM5 prediction."
                )

            if not self._is_missense():
                raise AlgorithmError("Variant is not a missense variant. PS1/PM5 not applicable.")

            self.comment = (
                "Extracting primary amino acid change from pHGVS: "
                f"{self.variant_info.dbnsfp.HGVSp_VEP}.\n"
            )
            primary_aa_change = self._parse_HGVSp(self.variant_info.dbnsfp.HGVSp_VEP)
            if not primary_aa_change:
                raise AlgorithmError("No valid primary amino acid change for PS1/PM5 prediction.")

            self.comment += f"Primary amino acid change: {primary_aa_change.name}. =>\n"
            for alt_base in DNA_BASES:
                # Skip the same base insert
                if alt_base == self.seqvar.insert:
                    continue

                self.comment += f"Analysing alternative variant with base: {alt_base}.\n"
                alt_seqvar = SeqVar(
                    genome_release=self.seqvar.genome_release,
                    chrom=self.seqvar.chrom,
                    pos=self.seqvar.pos,
                    delete=self.seqvar.delete,
                    insert=alt_base,
                )
                alt_info = self._get_variant_info(alt_seqvar)

                if alt_info and alt_info.result.dbnsfp and alt_info.result.dbnsfp.HGVSp_VEP:
                    self.comment += (
                        f"Extracting alternative amino acid change from pHGVS: "
                        f"{alt_info.result.dbnsfp.HGVSp_VEP}.\n"
                    )
                    alt_aa_change = self._parse_HGVSp(alt_info.result.dbnsfp.HGVSp_VEP)
                    self.comment += (
                        "Alternative amino acid change: "
                        f"{alt_aa_change.name if alt_aa_change else "N/A"}. =>\n"
                    )
                    if alt_aa_change and self._is_pathogenic(alt_info.result):
                        if primary_aa_change == alt_aa_change:
                            self.comment += (
                                "Alternative variant is pathogenic "
                                "and amino acid change is the same as primary variant.\n"
                                "Result: PS1 is met.\n"
                            )
                            self.prediction.PS1 = True  # Same amino acid change and pathogenic
                        if primary_aa_change != alt_aa_change:
                            self.comment += (
                                "Alternative variant is pathogenic "
                                "and amino acid change is different from primary variant.\n"
                                "Result: PM5 is met.\n"
                            )
                            self.prediction.PM5 = (
                                True  # Different amino acid change at same residue, pathogenic
                            )
                else:
                    self.comment += f"Failed to get variant information for {alt_seqvar}.\n"

        except AutoAcmgBaseException as e:
            logger.error("Error occurred during PS1/PM5 prediction. Error: {}", e)
            self.comment += f"Error occurred during PS1/PM5 prediction. Error: {e}"
            self.prediction = None

        # Return the prediction result and the comment with the explanation
        return self.prediction, self.comment
