"""Implementation of PS1 prediction for sequence variants."""

import re
from typing import Optional, Tuple

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.defs.annonars_variant import AnnonarsVariantResponse, VariantResult
from src.defs.auto_acmg import AminoAcid
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar

#: DNA bases
DNA_BASES = ["A", "C", "G", "T"]


class AutoPS1:
    """Predicts PS1 criteria for sequence variants."""

    def __init__(self, seqvar: SeqVar, genome_release: GenomeRelease):
        self.seqvar = seqvar
        self.genome_release = genome_release

    def get_variant_info(self, seqvar: SeqVar) -> Optional[AnnonarsVariantResponse]:
        """Get variant information from Annonars.

        Returns:
            dict: Annonars response.
        """
        try:
            logger.debug("Getting variant information for {}.", seqvar)
            annonars_client = AnnonarsClient()
            return annonars_client.get_variant_info(seqvar)
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

    def _is_pathogenic(self, variant_info: VariantResult) -> bool:
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
                variant_info.clinvar.referenceAssertions[0].clinicalSignificance
                == "CLINICAL_SIGNIFICANCE_PATHOGENIC"
            )
        return False

    def predict(self) -> Optional[bool]:
        """Predicts PS1 criteria for sequence variants.

        Returns:
            bool: True if the variant meets PS1 criteria, False if it does not, or None if
            prediction failed.
        """
        try:
            # Predict PS1 criteria
            variant_info = self.get_variant_info(self.seqvar)
            if not variant_info:
                logger.error("No variant info available for PS1 prediction.")
                return None

            # Parse the pHGVSp string
            if variant_info.result.dbnsfp and variant_info.result.dbnsfp.HGVSp_VEP:
                new_aa = self.parse_HGVSp(variant_info.result.dbnsfp.HGVSp_VEP)
            else:
                logger.debug("No pHGVSp available for sequence variant.")
                new_aa = None

            # Check similar variants for pathogenicity
            for alt_base in DNA_BASES:
                if alt_base == self.seqvar.delete:
                    continue
                alt_seqvar = SeqVar(
                    genome_release=self.genome_release,
                    chrom=self.seqvar.chrom,
                    pos=self.seqvar.pos,
                    delete=self.seqvar.delete,
                    insert=alt_base,
                )
                variant_info = self.get_variant_info(alt_seqvar)
                if variant_info:
                    if not variant_info.result.dbnsfp or not variant_info.result.dbnsfp.HGVSp_VEP:
                        logger.debug("No pHGVSp available for sequence variant: {}", alt_seqvar)
                        continue

                    # Parse the pHGVSp string for the alternative variant
                    alt_new_aa = self.parse_HGVSp(variant_info.result.dbnsfp.HGVSp_VEP)
                    # Check if the alternative variant is pathogenic
                    if self._is_pathogenic(variant_info.result) and (
                        not new_aa or new_aa == alt_new_aa
                    ):
                        logger.info(
                            "PS1 criteria met for sequence variant {} with new amino acid {}.",
                            self.seqvar.user_repr,
                            new_aa,
                        )
                        return True
            return False

        except Exception as e:
            logger.error("Failed to predict PS1 criteria. Error: {}", e)
            return None
