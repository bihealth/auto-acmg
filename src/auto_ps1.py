"""Implementation of PS1 prediction for sequence variants."""

import re
from typing import Optional

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.defs.annonars_variant import AnnonarsVariantResponse, VariantResult
from src.defs.auto_acmg import PS1PM5, AminoAcid
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar

#: DNA bases
DNA_BASES = ["A", "C", "G", "T"]


class AutoPS1PM5:
    """Predicts PS1 criteria for sequence variants."""

    def __init__(self, seqvar: SeqVar, genome_release: GenomeRelease):
        self.seqvar = seqvar
        self.genome_release = genome_release
        self.annonars_client = AnnonarsClient()

    def get_variant_info(self, seqvar: SeqVar) -> Optional[AnnonarsVariantResponse]:
        """Get variant information from Annonars.

        Returns:
            dict: Annonars response.
        """
        try:
            logger.debug("Getting variant information for {}.", seqvar)
            return self.annonars_client.get_variant_info(seqvar)
        except Exception as e:
            logger.error("Failed to get variant information. Error: {}", e)
            return None

    @staticmethod
    def parse_HGVSp(pHGVSp: str) -> Optional[str]:
        """Parse the pHGVSp string into its components.

        Args:
            pHGVSp (str): The protein change in HGVS format.

        Returns:
            str: The new amino acid if the pHGVSp is valid, None otherwise.
        """
        match = re.match(r"p\.(\D+)(\d+)(\D+)", pHGVSp)
        if match:
            # original_aa = AminoAcid[match.group(1)].value
            # position = int(match.group(2))
            new_aa = AminoAcid[match.group(3)].value
            return new_aa
        logger.debug("Invalid pHGVSp: {}", pHGVSp)
        return None

    @staticmethod
    def _is_pathogenic(variant_info: VariantResult) -> bool:
        """Check if the variant is pathogenic.

        Args:
            variant_info (dict): Annonars variant information

        Returns:
            bool: True if the variant is pathogenic
        """
        if (
            variant_info.clinvar
            and variant_info.clinvar.referenceAssertions
            and variant_info.clinvar.referenceAssertions[0].clinicalSignificance
        ):
            return (
                variant_info.clinvar.referenceAssertions[0].clinicalSignificance == "CLINICAL_SIGNIFICANCE_PATHOGENIC"
            )
        return False

    def predict(self) -> Optional[PS1PM5]:
        """Predicts the criteria PS1 and PM5 for the provided sequence variant.

        Returns:
            PS1PM5: The prediction result.
        """
        prediction = PS1PM5()

        primary_info = self.get_variant_info(self.seqvar)
        if not primary_info or not primary_info.result.dbnsfp or not primary_info.result.dbnsfp.HGVSp_VEP:
            logger.error("No primary variant information available for prediction.")
            return prediction

        primary_aa_change = self.parse_HGVSp(primary_info.result.dbnsfp.HGVSp_VEP)
        if not primary_aa_change:
            logger.debug("No valid primary amino acid change for PS1/PM5 prediction.")
            return prediction

        for alt_base in DNA_BASES:
            if alt_base == self.seqvar.insert:
                continue  # Skip the same base insert
            alt_seqvar = SeqVar(
                genome_release=self.genome_release,
                chrom=self.seqvar.chrom,
                pos=self.seqvar.pos,
                delete=self.seqvar.delete,
                insert=alt_base,
            )
            alt_info = self.get_variant_info(alt_seqvar)
            if alt_info and alt_info.result.dbnsfp and alt_info.result.dbnsfp.HGVSp_VEP:
                alt_aa_change = self.parse_HGVSp(alt_info.result.dbnsfp.HGVSp_VEP)
                if alt_aa_change and self._is_pathogenic(alt_info.result):
                    if primary_aa_change == alt_aa_change:
                        prediction.PS1 = True  # Same amino acid change and pathogenic
                    if primary_aa_change != alt_aa_change:
                        prediction.PM5 = True  # Different amino acid change at same residue, pathogenic

        return prediction
