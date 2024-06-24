"""Implementation of BA1, BS1, BS2, PM2 prediction for sequence variants."""

import re
from typing import List, Optional, Tuple

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.annonars_variant import AlleleCount, AnnonarsVariantResponse, VariantResult
from src.defs.auto_acmg import BA1BS1BS2PM2
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


class AutoBA1BS1BS2PM2:
    """Predicts BA1, BS1, BS2, PM2 criteria for sequence variants."""

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
        self.prediction: Optional[BA1BS1BS2PM2] = None
        #: Comment to store the prediction explanation.
        self.comment: str = ""

    def _get_control_af(self, variant_data: VariantResult) -> AlleleCount:
        """
        Get the allele frequency for the control population.

        Args:
            variant_data: The variant data.

        Returns:
            The allele frequency for the control population.
        """
        if not variant_data.gnomad_exomes or not variant_data.gnomad_exomes.alleleCounts:
            raise MissingDataError("No allele counts found in variant data")
        for af in variant_data.gnomad_exomes.alleleCounts:
            if af.cohort == "controls":
                return af
        raise MissingDataError("No controls allele counts found in variant data.")

    def _get_af(self, seqvar: SeqVar, variant_data: VariantResult) -> float:
        """
        Get the allele frequency for the sequence variant.

        Args:
            seqvar: The sequence variant.
            variant_data: The variant data.

        Returns:
            The allele frequency.
        """
        if seqvar.chrom.startswith("M"):
            if not variant_data.gnomad_mtdna or not variant_data.gnomad_mtdna.afHet:
                self.comment += "No gnomad data found for mitochondrial variant.\n"
                raise MissingDataError("No gnomad data found for mitochondrial variant.")
            else:
                self.comment += (
                    "Mitochondrial variant with allele frequency: "
                    f"{variant_data.gnomad_mtdna.afHet}.\n"
                )
                return variant_data.gnomad_mtdna.afHet
        else:
            controls_af = self._get_control_af(variant_data)
            if not controls_af.raw or not controls_af.raw.af:
                self.comment += "No allele frequency found in control data.\n"
                raise MissingDataError("No allele frequency found in control data.")
            self.comment += f"Allele frequency: {controls_af.raw.af}.\n"
            return controls_af.raw.af

    def _check_zygosity(self, seqvar: SeqVar, variant_data: VariantResult) -> bool:
        """
        Check the zygosity of the sequence variant.

        Args:
            variant_data: The variant data.

        Returns:
            True if the variant is recessive (homozygous), dominant (heterozygous), or X-linked
            (hemizygous) disorder.
        """
        if seqvar.chrom.startswith("M"):
            self.comment += "Mitochondrial variants are not considered for BS2 criteria."
            return False
        controls_af = self._get_control_af(variant_data)
        if not controls_af.raw or not controls_af.raw.nhomalt:
            self.comment += "No nhomalt found in control data.\n"
            raise MissingDataError("No raw data found in control data.")
        if controls_af.raw.nhomalt > 5:
            self.comment += f"nhomalt is {controls_af.raw.nhomalt} > 5.\n"
            return True
        self.comment += f"nhomalt is {controls_af.raw.nhomalt} <= 5.\n"
        return False

    def predict(self) -> Tuple[Optional[BA1BS1BS2PM2], str]:
        """
        Predicts the BA1, BS1, BS2, PM2 criteria for the sequence variant.

        Note:
            Rules:
            BA1: Allele frequency is >5%.

            BS1: Allele frequency is between 1% and 5%.

            BS2: Observed in a healthy adult individual for a recessive (homozygous), dominant
            (heterozygous), or X-linked (hemizygous) disorder, with full penetrance expected at an
            early age.

            PM2: Absent from controls allele frequency data.

        Returns:
            BA1BS1BS2PM2: The prediction result.
        """
        try:
            self.prediction = BA1BS1BS2PM2()

            self.comment += "Check allele frequency for the control population.\n"
            af = self._get_af(self.seqvar, self.variant_info)
            if af > 0.05:
                self.prediction.BA1 = True
            elif af >= 0.01:
                self.prediction.BS1 = True
            else:
                self.prediction.PM2 = True

            self.comment += "Check zygosity.\n"
            if af >= 0.01 and self._check_zygosity(self.seqvar, self.variant_info):
                self.prediction.BS2 = True

        except AutoAcmgBaseException as e:
            logger.error("Error occurred during BA1, BS1, BS2, PM2 prediction. Error: {}", e)
            self.comment += f"An error occurred while predicting BA1, BS1, BS2, PM2 criteria: {e}"
            self.prediction = None

        # Return the prediction result and explanation
        return self.prediction, self.comment
