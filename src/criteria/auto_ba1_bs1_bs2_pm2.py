"""Implementation of BA1, BS1, BS2, PM2 prediction for sequence variants."""

import re
from typing import Optional

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.annonars_variant import AnnonarsVariantResponse, VariantResult
from src.defs.auto_acmg import BA1BS1BS2PM2
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar


class AutoBA1BS1BS2PM2:
    """Predicts BA1, BS1, BS2, PM2 criteria for sequence variants."""

    def __init__(
        self,
        seqvar: SeqVar,
        genome_release: GenomeRelease,
        variant_info: VariantResult,
        *,
        config: Optional[Config] = None,
    ):
        #: Configuration to use.
        self.config = config or Config()
        self.seqvar = seqvar
        self.genome_release = genome_release
        self.variant_info = variant_info
        self.annonars_client = AnnonarsClient(api_base_url=self.config.api_base_url_annonars)
        self.prediction: BA1BS1BS2PM2 | None = None

    def _get_variant_info(self, seqvar: SeqVar) -> Optional[AnnonarsVariantResponse]:
        """Get variant information from Annonars.

        Returns:
            dict: Annonars response.
        """
        try:
            logger.debug("Getting variant information for {}.", seqvar)
            return self.annonars_client.get_variant_info(seqvar)
        except AutoAcmgBaseException as e:
            logger.error("Failed to get variant information. Error: {}", e)
            return None

    def evaluate_ba1(self, seqvar: SeqVar, variant_data: VariantResult) -> bool:
        if (
            not variant_data.gnomad_exomes
            or not variant_data.gnomad_exomes.alleleCounts
            or not variant_data.gnomad_exomes.alleleCounts[0].raw
            or not variant_data.gnomad_exomes.alleleCounts[0].raw.af
        ):
            raise MissingDataError("No allele counts found in variant data")

        chrom = seqvar.chrom
        allele_freq = variant_data.gnomad_exomes.alleleCounts[0].raw.af

        if chrom.startswith("chrM"):
            return allele_freq > 0.01
        else:
            return allele_freq > 0.05

    def evaluate_bs1(self, seqvar: SeqVar, variant_data: VariantResult, max_credible_freq) -> bool:
        if (
            not variant_data.gnomad_exomes
            or not variant_data.gnomad_exomes.alleleCounts
            or not variant_data.gnomad_exomes.alleleCounts[0].byAncestryGroup
            or not variant_data.gnomad_exomes.alleleCounts[0].raw
            or not variant_data.gnomad_exomes.alleleCounts[0].raw.af
        ):
            raise MissingDataError("No allele counts found in variant data")
        chrom = seqvar.chrom
        faf = 0.0
        for group in variant_data.gnomad_exomes.alleleCounts[0].byAncestryGroup:
            if group.faf95 and group.faf95 > faf:
                faf = group.faf95

        if chrom.startswith("chrM"):
            return variant_data.gnomad_exomes.alleleCounts[0].raw.af > 0.005
        else:
            return faf > max_credible_freq

    def evaluate_bs2(self, variant_data: VariantResult) -> bool:
        if (
            not variant_data.gnomad_exomes
            or not variant_data.gnomad_exomes.alleleCounts
            or not variant_data.gnomad_exomes.alleleCounts[0].raw
            or not variant_data.gnomad_exomes.alleleCounts[0].raw.ac
            or not variant_data.gnomad_exomes.alleleCounts[0].raw.nhomalt
        ):
            raise MissingDataError("No allele counts found in variant data")

        gene_mode = "inplace"  # Field to identify recessive/dominant/X-linked
        allele_count = variant_data.gnomad_exomes.alleleCounts[0].raw.ac

        if gene_mode in ["recessive", "X-linked"]:
            return allele_count > 2
        elif gene_mode == "dominant":
            return allele_count > 5
        return False

    def evaluate_pm2(self, seqvar: SeqVar, variant_data: VariantResult) -> bool:
        if (
            not variant_data.gnomad_exomes
            or not variant_data.gnomad_exomes.alleleCounts
            or not variant_data.gnomad_exomes.alleleCounts[0].raw
            or not variant_data.gnomad_exomes.alleleCounts[0].raw.af
            or not variant_data.gnomad_exomes.alleleCounts[0].raw.ac
            or not variant_data.gnomad_exomes.alleleCounts[0].raw.nhomalt
        ):
            raise MissingDataError("No allele counts found in variant data")

        chrom = seqvar.chrom
        allele_freq = variant_data.gnomad_exomes.alleleCounts[0].raw.af
        gene_mode = "inplace"  # Field to identify recessive/dominant/X-linked

        if chrom.startswith("chrM"):
            return allele_freq < 0.00002
        else:
            if gene_mode in ["recessive", "X-linked"]:
                allele_count = variant_data.gnomad_exomes.alleleCounts[0].raw.ac
                if allele_count <= 4:
                    return True
            elif gene_mode == "dominant":
                homozygous_count = variant_data.gnomad_exomes.alleleCounts[0].raw.nhomalt
                if homozygous_count <= 1:
                    return True
            return allele_freq < 0.0001

    def predict(self) -> Optional[BA1BS1BS2PM2]:
        """
        Predicts the BA1, BS1, BS2, PM2 criteria for the sequence variant.


        Note:
            Rules:
            BA1: Allele frequency is >5% in Exome Sequencing Project, 1000 Genomes Project, or Exome
            Aggregation Consortium.

            BS1: Allele frequency greater than expected for disorder.

            BS2: Observed in a healthy adult individual for a recessive (homozygous), dominant
            (heterozygous), or X-linked (hemizygous) disorder, with full penetrance expected at an
            early age.


            PM2: Absent from controls (or at extremely low frequency if recessive) in Exome
            Sequencing Project, 1000 Genomes or ExAC.

        Returns:
            BA1BS1BS2PM2: The prediction result.
        """
        try:
            # Instantiate the prediction result
            self.prediction = BA1BS1BS2PM2()

            # Evaluate each criterion
            self.prediction.BA1 = self.evaluate_ba1(self.seqvar, self.variant_info)
            self.prediction.BS1 = self.evaluate_bs1(
                self.seqvar, self.variant_info, max_credible_freq=0.01
            )  # Adjust max_credible_freq as needed
            self.prediction.BS2 = self.evaluate_bs2(self.variant_info)
            self.prediction.PM2 = self.evaluate_pm2(self.seqvar, self.variant_info)
        except AutoAcmgBaseException as e:
            logger.error("Error occurred during BA1, BS1, BS2, PM2 prediction. Error: {}", e)
            self.prediction = None

        # Return the prediction result
        return self.prediction
