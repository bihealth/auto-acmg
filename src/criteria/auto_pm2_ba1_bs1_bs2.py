"""Implementation of BA1, BS1, BS2, PM2 prediction for sequence variants."""

import re
from typing import List, Optional, Tuple

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.annonars_variant import (
    AlleleCount,
    AnnonarsVariantResponse,
    GnomadExomes,
    GnomadMtDna,
    VariantResult,
)
from src.defs.auto_acmg import (
    PM2BA1BS1BS2,
    AlleleCondition,
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
    ClingenDosageMap,
)
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.seqvar import SeqVar
from src.utils import AutoACMGHelper, SeqVarTranscriptsHelper


class AutoPM2BA1BS1BS2(AutoACMGHelper):
    """Predicts PM2, BA1, BS1, BS2 criteria for sequence variants."""

    def __init__(self):
        super().__init__()
        #: Prediction result.
        self.prediction_pm2ba1bs1bs2: Optional[PM2BA1BS1BS2] = None
        #: comment_pm2ba1bs1bs2 to store the prediction explanation.
        self.comment_pm2ba1bs1bs2: str = ""

    def _get_control_af(self, gnomad_exomes: Optional[GnomadExomes]) -> Optional[AlleleCount]:
        """
        Get the allele frequency information for the control population.

        Args:
            variant_data: The variant data.

        Returns:
            The allele frequency for the control population. None if no data found.
        """
        if not gnomad_exomes or not gnomad_exomes.alleleCounts:
            return None
            # raise MissingDataError("No allele counts found in variant data")
        for af in gnomad_exomes.alleleCounts:
            if af.cohort == "controls":
                return af
        return None

    def _get_any_af(self, gnomad_exomes: Optional[GnomadExomes]) -> Optional[AlleleCount]:
        """
        Get the allele frequency information for any population.

        Args:
            variant_data: The variant data.

        Returns:
            The allele frequency for any population. None if no data found.
        """
        if not gnomad_exomes or not gnomad_exomes.alleleCounts:
            self.comment_pm2ba1bs1bs2 += "Missing gnomad data."
            return None
            # raise MissingDataError("No allele counts found in variant data")
        best_af = None
        for af in gnomad_exomes.alleleCounts:
            if not best_af:
                best_af = af
            elif af.cohort == "controls":
                return af
            else:
                if best_af.afGrpmax < af.afGrpmax:
                    best_af = af
        return best_af

    def _get_af(
        self,
        seqvar: SeqVar,
        gnomad_mtdna: Optional[GnomadMtDna],
        gnomad_exomes: Optional[GnomadExomes],
    ) -> Optional[float]:
        """
        Get the allele frequency for the sequence variant.

        Args:
            seqvar: The sequence variant.
            variant_data: The variant data.

        Returns:
            The allele frequency. None if no controls data
        """
        if seqvar.chrom.startswith("M"):
            if not gnomad_mtdna or not gnomad_mtdna.afHet:
                self.comment_pm2ba1bs1bs2 += "No gnomad data found for mitochondrial variant.\n"
                return None
                # raise MissingDataError("No gnomad data found for mitochondrial variant.")
            else:
                self.comment_pm2ba1bs1bs2 += (
                    "Mitochondrial variant with allele frequency: " f"{gnomad_mtdna.afHet}.\n"
                )
                return gnomad_mtdna.afHet
        else:
            controls_af = self._get_control_af(gnomad_exomes=gnomad_exomes)
            any_af = self._get_any_af(gnomad_exomes=gnomad_exomes)
            af = controls_af or any_af
            if not af or not af.afGrpmax:
                self.comment_pm2ba1bs1bs2 += "No allele frequency data found."
                return None
                # raise MissingDataError("No allele frequency found in data.")
            self.comment_pm2ba1bs1bs2 += f"Allele frequency: {af.afGrpmax}.\n"
            return af.afGrpmax

    def _get_allele_cond(self, seqvar: SeqVar) -> AlleleCondition:
        """
        Get the allele condition for the sequence variant.

        Get the Clingen dosage for the gene from the gene transcript data (mehari).
        If the Clingen dosage is unknown, try Decipher and Domino data.
        We use the following thresholds:
        - Decipher: pHi >= 0.9 for dominant
        - Domino: score > 0.5934 for dominant, score < 0.3422 for recessive

        Args:
            seqvar: The sequence variant.

        Returns:
            The allele condition.
        """
        # Fetch transcript data
        seqvar_transcript_helper = SeqVarTranscriptsHelper(seqvar, config=self.config)
        seqvar_transcript_helper.initialize()
        (
            _,
            gene_transcript,
            _,
            _,
            _,
        ) = seqvar_transcript_helper.get_ts_info()
        if not gene_transcript or not gene_transcript.geneId:
            raise AlgorithmError("No gene found for the transcript.")
        gene_info = self.annonars_client.get_gene_info(gene_transcript.geneId)
        if not gene_info or not gene_info.genes or not gene_info.genes.root:
            raise AlgorithmError("No clingen gene information found for the gene.")
        # Get the Clingen dosage
        clingen_dosage = AlleleCondition.Unknown
        if (
            not (gene_data := gene_info.genes.root.get(gene_transcript.geneId))
            or not (clingen_data := gene_data.clingen)
            or not clingen_data.haploinsufficiencyScore
        ):
            self.comment_pm2ba1bs1bs2 += "No Clingen dosage information found for the gene.\n"
        else:
            clingen_dosage = ClingenDosageMap[clingen_data.haploinsufficiencyScore]
        # Try Decipher and Domino if Clingen dosage is unknown
        if clingen_dosage == AlleleCondition.Unknown:
            # Decipher
            if (
                not (gene_data := gene_info.genes.root.get(gene_transcript.geneId))
                or not (decipher_data := gene_data.decipherHi)
                or not decipher_data.pHi
            ):
                self.comment_pm2ba1bs1bs2 += "No Decipher information found for the gene.\n"
            else:
                if decipher_data.pHi >= 0.9:
                    clingen_dosage = AlleleCondition.Dominant
            # Domino
            if (
                not (gene_data := gene_info.genes.root.get(gene_transcript.geneId))
                or not (domino_data := gene_data.domino)
                or not domino_data.score
            ):
                self.comment_pm2ba1bs1bs2 += "No Domino information found for the gene.\n"
            else:
                if domino_data.score > 0.5934:
                    clingen_dosage = AlleleCondition.Dominant
                elif domino_data.score < 0.3422:
                    clingen_dosage = AlleleCondition.Recessive
        return clingen_dosage

    def _check_zyg(self, seqvar: SeqVar, gnomad_exomes: Optional[GnomadExomes]) -> bool:
        """
        Check the zygosity of the sequence variant.

        If the variant is mitochondrial, it is not considered for BS2 criteria.
        Otherwise, parse the allele condition and check the zygosity:
        If the variant is on X chromosome, check the allele count for XX and XY as follows:
        - for dominant: XX allele count - 2 * XX nhomalt + XY allele count > 2
        - for recessive: XX nhomalt + XY nhomalt > 2
        - for dominant/recessive: XX allele count - 2 * XX nhomalt + XY allele count > 2 and
        XX nhomalt + XY nhomalt > 2
        If the variant is on autosomal chromosomes, check the allele count as follows:
        - for dominant: allele count - 2 * nhomalt > 5
        - for recessive: nhomalt > 5
        - for dominant/recessive: allele count - 2 * nhomalt > 5 and nhomalt > 5
        Return True if the variant is in a recessive (homozygous), dominant (heterozygous), or
        X-linked (hemizygous) disorder (condition is met).

        Args:
            variant_data: The variant data.

        Returns:
            True if the variant is recessive (homozygous), dominant (heterozygous), or X-linked
            (hemizygous) disorder.
        """
        if seqvar.chrom.startswith("M"):
            self.comment_pm2ba1bs1bs2 += (
                "Mitochondrial variants are not considered for BS2 criteria."
            )
            return False

        allele_condition = self._get_allele_cond(seqvar)
        self.comment_pm2ba1bs1bs2 += f"Allele condition: {allele_condition.name}.\n"
        controls_af = self._get_control_af(gnomad_exomes)
        any_af = self._get_any_af(gnomad_exomes)
        af = controls_af or any_af
        if not af or not af.bySex:
            self.comment_pm2ba1bs1bs2 += "No controls allele data found in control data.\n"
            raise MissingDataError("No raw data found in control data.")

        # X-linked disorders
        if seqvar.chrom == "X":
            if not af.bySex.xx or not af.bySex.xy:
                self.comment_pm2ba1bs1bs2 += "No allele data found for XX or XY in control data.\n"
                raise MissingDataError("No allele data found for XX or XY in control data.")
            xx_ac = af.bySex.xx.ac if af.bySex.xx.ac else 0
            xy_ac = af.bySex.xy.ac if af.bySex.xy.ac else 0
            xx_nhomalt = af.bySex.xx.nhomalt if af.bySex.xx.nhomalt else 0
            xy_nhomalt = af.bySex.xy.nhomalt if af.bySex.xy.nhomalt else 0
            self.comment_pm2ba1bs1bs2 += (
                f"Allele count for XX: {xx_ac}, XY: {xy_ac}, "
                f"Nhomalt for XX: {xx_nhomalt}, XY: {xy_nhomalt}.\n"
            )
            if allele_condition == AlleleCondition.Dominant:
                if xx_ac - 2 * xx_nhomalt + xy_ac > 2:
                    self.comment_pm2ba1bs1bs2 += (
                        "XX allele count - 2 * XX nhomalt + XY allele count "
                        f"({xx_ac - 2 * xx_nhomalt + xy_ac}) > 2.\n"
                        "The variant is in a dominant X-linked disorder."
                    )
                    return True
            elif allele_condition == AlleleCondition.Recessive:
                if xx_nhomalt + xy_nhomalt > 2:
                    self.comment_pm2ba1bs1bs2 += (
                        "XX nhomalt + XY nhomalt "
                        f"({xx_nhomalt + xy_nhomalt}) > 2.\n"
                        "The variant is in a recessive X-linked disorder."
                    )
                    return True
            else:
                if xx_ac - 2 * xx_nhomalt + xy_ac > 2 and xx_ac + xy_ac > 2:
                    self.comment_pm2ba1bs1bs2 += (
                        "XX allele count - 2 * XX nhomalt + XY allele count "
                        f"({xx_ac - 2 * xx_nhomalt + xy_ac}) > 2.\n"
                        "XX nhomalt + XY nhomalt "
                        f"({xx_nhomalt + xy_nhomalt}) > 2.\n"
                        "The variant is in a dominant/recessive X-linked disorder."
                    )
                    return True
        # Autosomal disorders
        else:
            if not af.bySex.overall:
                self.comment_pm2ba1bs1bs2 += "No allele data found for overall in control data.\n"
                raise MissingDataError("No allele data found for overall in control data.")
            ac = af.bySex.overall.ac if af.bySex.overall.ac else 0
            nhomalt = af.bySex.overall.nhomalt if af.bySex.overall.nhomalt else 0
            self.comment_pm2ba1bs1bs2 += f"Allele count: {ac}, Nhomalt: {nhomalt}.\n"
            if allele_condition == AlleleCondition.Dominant:
                if ac - 2 * nhomalt > 5:
                    self.comment_pm2ba1bs1bs2 += (
                        "Allele count - 2 * Nhomalt "
                        f"({ac - 2 * nhomalt}) > 5.\n"
                        "The variant is in a dominant (heterozygous) disorder."
                    )
                    return True
            elif allele_condition == AlleleCondition.Recessive:
                if nhomalt > 5:
                    self.comment_pm2ba1bs1bs2 += (
                        f"Nhomalt {nhomalt} > 0.\n"
                        "The variant is in a recessive (homozygous) disorder."
                    )
                    return True
            else:
                if ac - 2 * nhomalt > 5 and nhomalt > 5:
                    self.comment_pm2ba1bs1bs2 += (
                        "Allele count - 2 * Nhomalt "
                        f"({ac - 2 * nhomalt}) > 5.\n"
                        f"Nhomalt {nhomalt} > 5.\n"
                        "The variant is in a dominant/recessive disorder."
                    )
                    return True
        return False

    def _ba1_exception(self, seqvar: SeqVar) -> bool:
        """
        Check the exception for BA1 criteria.

        If the variant in exception list, return True.

        Args:
            seqvar: The sequence variant.

        Returns:
            True if the variant is in exception list.
        """
        exception_list = [
            SeqVar(
                genome_release=GenomeRelease.GRCh37,
                chrom="chr3",
                pos=128598490,
                delete="C",
                insert="CTAAG",
                user_repr="NM_014049.4:c.-44_-41dupTAAG",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh38,
                chrom="chr3",
                pos=128879647,
                delete="C",
                insert="CTAAG",
                user_repr="NM_014049.4:c.-44_-41dupTAAG",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh37,
                chrom="chr13",
                pos=20763612,
                delete="C",
                insert="T",
                user_repr="NM_004004.5:c.109G>A",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh38,
                chrom="chr13",
                pos=20189473,
                delete="C",
                insert="T",
                user_repr="NM_004004.5:c.109G>A",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh37,
                chrom="chr6",
                pos=26091179,
                delete="C",
                insert="G",
                user_repr="NM_000410.3:c.187C>G",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh38,
                chrom="chr6",
                pos=26090951,
                delete="C",
                insert="G",
                user_repr="NM_000410.3:c.187C>G",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh37,
                chrom="chr6",
                pos=26093141,
                delete="G",
                insert="A",
                user_repr="NM_000410.3:c.845G>A",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh38,
                chrom="chr6",
                pos=26092913,
                delete="G",
                insert="A",
                user_repr="NM_000410.3:c.845G>A",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh37,
                chrom="chr16",
                pos=3299586,
                delete="G",
                insert="A",
                user_repr="NM_000243.2:c.1105C>T",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh38,
                chrom="chr16",
                pos=3249586,
                delete="G",
                insert="A",
                user_repr="NM_000243.2:c.1105C>T",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh37,
                chrom="chr16",
                pos=3299468,
                delete="C",
                insert="T",
                user_repr="NM_000243.2:c.1223G>A",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh38,
                chrom="chr16",
                pos=3249468,
                delete="C",
                insert="T",
                user_repr="NM_000243.2:c.1223G>A",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh37,
                chrom="chr13",
                pos=73409497,
                delete="G",
                insert="A",
                user_repr="NM_006346.2:c.1214G>A",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh38,
                chrom="chr13",
                pos=72835359,
                delete="G",
                insert="A",
                user_repr="NM_006346.2:c.1214G>A",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh37,
                chrom="chr12",
                pos=121175678,
                delete="C",
                insert="T",
                user_repr="NM_000017.3:c.511C>T",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh38,
                chrom="chr12",
                pos=120737875,
                delete="C",
                insert="T",
                user_repr="NM_000017.3:c.511C>T",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh37,
                chrom="chr3",
                pos=15686693,
                delete="G",
                insert="C",
                user_repr="NM_000060.4:c.1330G>C",
            ),
            SeqVar(
                genome_release=GenomeRelease.GRCh38,
                chrom="chr3",
                pos=15645186,
                delete="G",
                insert="C",
                user_repr="NM_000060.4:c.1330G>C",
            ),
        ]
        if seqvar in exception_list:
            self.comment_pm2ba1bs1bs2 += "The variant is in the exception list for BA1 criteria."
            return True
        return False

    def verify_pm2ba1bs1bs2(
        self,
        seqvar: SeqVar,
        gnomad_exones: Optional[GnomadExomes],
        gnomad_mtdna: Optional[GnomadMtDna],
    ) -> Tuple[Optional[PM2BA1BS1BS2], str]:
        """
        Predicts the PM2, BA1, BS1, BS2 criteria for the sequence variant.

        Note:
            Rules:
            PM2: Absent from controls allele frequency data.

            BA1: Allele frequency is >5%.

            BS1: Allele frequency is between 1% and 5%.

            BS2: Observed in a healthy adult individual for a recessive (homozygous), dominant
            (heterozygous), or X-linked (hemizygous) disorder, with full penetrance expected at an
            early age.

        Returns:
            BA1BS1BS2PM2: The prediction result.
        """
        self.prediction_pm2ba1bs1bs2 = PM2BA1BS1BS2()
        self.comment_pm2ba1bs1bs2 = ""
        try:
            self.comment_pm2ba1bs1bs2 += "Check allele frequency for the control population.\n"
            af = self._get_af(seqvar, gnomad_mtdna, gnomad_exones)
            if not af:
                self.comment_pm2ba1bs1bs2 += "No allele frequency data found.\n"
            elif af >= 0.05 and not self._ba1_exception(seqvar):
                self.comment_pm2ba1bs1bs2 += "Allele frequency > 5%: BA1 is met."
                self.prediction_pm2ba1bs1bs2.BA1 = True
            elif af >= 0.00015 and not self._ba1_exception(seqvar):
                self.comment_pm2ba1bs1bs2 += "Allele frequency > 1%: BS1 is met."
                self.prediction_pm2ba1bs1bs2.BS1 = True
            elif af <= 0.0001:
                self.comment_pm2ba1bs1bs2 += "Allele frequency <= 1%: PM2 is met."
                self.prediction_pm2ba1bs1bs2.PM2 = True

            self.comment_pm2ba1bs1bs2 += "Check zygosity.\n"
            if (
                not self.prediction_pm2ba1bs1bs2.BA1
                and af
                and self._check_zyg(seqvar, gnomad_exones)
            ):
                self.prediction_pm2ba1bs1bs2.BS2 = True

        except AutoAcmgBaseException as e:
            logger.error("Error occurred during PM2, BA1, BS1, BS2 prediction. Error: {}", e)
            self.comment_pm2ba1bs1bs2 = (
                f"An error occurred while predicting PM2, BA1, BS1, BS2 criteria: {e}"
            )
            self.prediction_pm2ba1bs1bs2 = None

        return self.prediction_pm2ba1bs1bs2, self.comment_pm2ba1bs1bs2

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """
        Predicts the PM2, BA1, BS1, BS2 criteria for the sequence variant.

        Args:
            seqvar: The sequence variant.
            var_data: The variant data.

        Returns:
            The prediction result.
        """
        pred, comment = self.verify_pm2ba1bs1bs2(
            seqvar, var_data.gnomad_exomes, var_data.gnomad_mtdna
        )
        if pred:
            pm2_pred = (
                AutoACMGPrediction.Met
                if pred.PM2
                else (AutoACMGPrediction.NotMet if pred.PM2 is False else AutoACMGPrediction.Failed)
            )
            ba1_pred = (
                AutoACMGPrediction.Met
                if pred.BA1
                else (AutoACMGPrediction.NotMet if pred.BA1 is False else AutoACMGPrediction.Failed)
            )
            bs1_pred = (
                AutoACMGPrediction.Met
                if pred.BS1
                else (AutoACMGPrediction.NotMet if pred.BS1 is False else AutoACMGPrediction.Failed)
            )
            bs2_pred = (
                AutoACMGPrediction.Met
                if pred.BS2
                else (AutoACMGPrediction.NotMet if pred.BS2 is False else AutoACMGPrediction.Failed)
            )
        else:
            pm2_pred = AutoACMGPrediction.Failed
            ba1_pred = AutoACMGPrediction.Failed
            bs1_pred = AutoACMGPrediction.Failed
            bs2_pred = AutoACMGPrediction.Failed
        return (
            AutoACMGCriteria(
                name="PM2",
                prediction=pm2_pred,
                strength=pred.PM2_strength if pred else AutoACMGStrength.PathogenicModerate,
                summary=comment,
            ),
            AutoACMGCriteria(
                name="BA1",
                prediction=ba1_pred,
                strength=pred.BA1_strength if pred else AutoACMGStrength.BenignStandAlone,
                summary=comment,
            ),
            AutoACMGCriteria(
                name="BS1",
                prediction=bs1_pred,
                strength=pred.BS1_strength if pred else AutoACMGStrength.BenignStrong,
                summary=comment,
            ),
            AutoACMGCriteria(
                name="BS2",
                prediction=bs2_pred,
                strength=pred.BS2_strength if pred else AutoACMGStrength.BenignStrong,
                summary=comment,
            ),
        )
