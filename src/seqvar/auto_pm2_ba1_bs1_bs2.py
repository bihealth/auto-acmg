"""Implementation of BA1, BS1, BS2, PM2 prediction for sequence variants."""

from typing import Optional, Tuple

from loguru import logger

from src.defs.annonars_variant import AlleleCount
from src.defs.auto_acmg import (
    BA1_ESCEPTION_LIST,
    PM2BA1BS1BS2,
    AlleleCondition,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    ClingenDosageMap,
)
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, MissingDataError
from src.defs.seqvar import SeqVar
from src.utils import AutoACMGHelper, SeqVarTranscriptsHelper


class AutoPM2BA1BS1BS2(AutoACMGHelper):
    """Class for PM2, BA1, BS1, BS2 prediction."""

    def __init__(self):
        super().__init__()
        #: Prediction result.
        self.prediction_pm2ba1bs1bs2: Optional[PM2BA1BS1BS2] = None
        #: comment_pm2ba1bs1bs2 to store the prediction explanation.
        self.comment_pm2ba1bs1bs2: str = ""

    def _get_control_af(self, var_data: AutoACMGSeqVarData) -> Optional[AlleleCount]:
        """
        Get the allele frequency information for the control population.

        Args:
            var_data: The variant data.

        Returns:
            Optional[AlleleCount]: The allele frequency for the control population. None if no data
            found.
        """
        if not var_data.gnomad_exomes or not var_data.gnomad_exomes.alleleCounts:
            return None
        for af in var_data.gnomad_exomes.alleleCounts:
            if af.cohort == "controls":
                return af
        return None

    def _get_any_af(self, var_data: AutoACMGSeqVarData) -> Optional[AlleleCount]:
        """
        Get the highest allele frequency information for any population. The control group has
        priority.

        Args:
            var_data: The variant data.

        Returns:
            Optional[AlleleCount]: The highest allele frequency for any population. None if no data
            found.
        """
        if not var_data.gnomad_exomes or not var_data.gnomad_exomes.alleleCounts:
            return None
        max_af = None
        for af in var_data.gnomad_exomes.alleleCounts:
            if not max_af and af.anGrpmax and af.anGrpmax > var_data.thresholds.an_min:
                max_af = af
            elif af.cohort == "controls":
                return af
            else:
                if (
                    max_af
                    and max_af.afGrpmax
                    and af.afGrpmax
                    and af.anGrpmax
                    and max_af.afGrpmax < af.afGrpmax
                    and af.anGrpmax > var_data.thresholds.an_min
                ):
                    max_af = af
        return max_af

    def _get_af(
        self,
        seqvar: SeqVar,
        var_data: AutoACMGSeqVarData,
    ) -> Optional[float]:
        """
        Get the allele frequency for the sequence variant.

        Args:
            seqvar: The sequence variant.
            variant_data: The variant data.

        Returns:
            Optional[float]: The allele frequency. None if no controls data
        """
        controls_af = self._get_control_af(var_data)
        any_af = self._get_any_af(var_data)
        af = controls_af or any_af
        if not af or not af.afGrpmax:
            raise MissingDataError("No allele frequency found in gnomad data.")
        return af.afGrpmax

    def _get_m_af(
        self,
        var_data: AutoACMGSeqVarData,
    ) -> Optional[float]:
        """
        Get the allele frequency for the mitochondrial sequence variant.

        Args:
            variant_data: The variant data.

        Returns:
            Optional[float]: The allele frequency. None if no controls data
        """
        if not var_data.gnomad_mtdna or not var_data.gnomad_mtdna.afHet:
            raise MissingDataError("No allele frequency found in mitochondrial gnomad data.")
        else:
            return var_data.gnomad_mtdna.afHet

    def _get_allele_cond(self, seqvar: SeqVar) -> AlleleCondition:
        """
        Get the allele condition for the sequence variant.

        Get the Clingen dosage for the gene from the gene transcript data (mehari).
        If the Clingen dosage is unknown, try Decipher and Domino scores. Compare the scores to
        specific thresholds to determine the allele condition.

        Args:
            seqvar: The sequence variant.

        Returns:
            AlleleCOndition: The allele condition.
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

    def _check_zyg(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
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
            bool: True if the variant is recessive (homozygous), dominant (heterozygous), or
            X-linked (hemizygous) disorder.
        """
        if seqvar.chrom.startswith("M"):
            self.comment_pm2ba1bs1bs2 += (
                "Mitochondrial variants are not considered for BS2 criteria."
            )
            return False

        allele_condition = self._get_allele_cond(seqvar)
        self.comment_pm2ba1bs1bs2 += f"Allele condition: {allele_condition.name}.\n"
        controls_af = self._get_control_af(var_data)
        any_af = self._get_any_af(var_data)
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
        Check the exception for BA1 criteria, specified by VCEP modification. If the variant in the
        exception list, return True.

        Args:
            seqvar: The sequence variant.

        Returns:
            bool: True if the variant is in exception list.
        """
        if seqvar in BA1_ESCEPTION_LIST:
            return True
        return False

    def _bs2_not_applicable(self, var_data: AutoACMGSeqVarData) -> bool:
        """
        Check if the BS2 criteria is not applicable.

        Per default, the BS2 criteria is applicable. Only some specific VCEP modifications can
        exclude the BS2 criteria.

        Args:
            seqvar: The sequence variant.

        Returns:
            bool: True if the BS2 criteria is not applicable.
        """
        return False

    def verify_pm2ba1bs1bs2(
        self,
        seqvar: SeqVar,
        var_data: AutoACMGSeqVarData,
    ) -> Tuple[Optional[PM2BA1BS1BS2], str]:
        """
        Predicts the PM2, BA1, BS1, BS2 criteria for the sequence variant.

        Predict criteria by checking the allele frequency data and comparing it to the thresholds.
        Assign PM2, BA1 and BS1 criteria based on the allele frequency data. For BS2 criteria,
        check zygosity of the variant.

        Note:
            Rules:
            PM2: Absent from controls allele frequency data.

            BA1: Allele frequency is >5%.

            BS1: Allele frequency is between 1% and 5%.

            BS2: Observed in a healthy adult individual for a recessive (homozygous), dominant
            (heterozygous), or X-linked (hemizygous) disorder, with full penetrance expected at an
            early age.

        Args:
            seqvar: The sequence variant.
            var_data: The variant data.

        Returns:
            Tuple[Optional[PM2BA1BS1BS2], str]: The prediction result and the explanation.
        """
        self.prediction_pm2ba1bs1bs2 = PM2BA1BS1BS2()
        self.comment_pm2ba1bs1bs2 = ""
        try:
            if seqvar.chrom == "MT":
                af = self._get_m_af(var_data)
                if not af:
                    self.comment_pm2ba1bs1bs2 = "No allele frequency data found. "
                elif af <= 0.00002:
                    self.comment_pm2ba1bs1bs2 = "Allele frequency <= 0.002%: PM2 is met. "
                    self.prediction_pm2ba1bs1bs2.PM2 = True
                elif af > 0.01:
                    self.comment_pm2ba1bs1bs2 = "Allele frequency > 1%: BA1 is met. "
                    self.prediction_pm2ba1bs1bs2.BA1 = True
                elif af > 0.005:
                    self.comment_pm2ba1bs1bs2 = "Allele frequency > 0.5%: BS1 is met. "
                    self.prediction_pm2ba1bs1bs2.BS1 = True
                return self.prediction_pm2ba1bs1bs2, self.comment_pm2ba1bs1bs2

            af = self._get_af(seqvar, var_data)
            if not af:
                self.comment_pm2ba1bs1bs2 = "No allele frequency data found. "
            elif self._ba1_exception(seqvar):
                self.comment_pm2ba1bs1bs2 = "The variant is in the exception list for BA1 criteria."
                self.prediction_pm2ba1bs1bs2.BA1 = False
                self.prediction_pm2ba1bs1bs2.BS1 = False
            elif af >= var_data.thresholds.ba1_benign:
                self.comment_pm2ba1bs1bs2 = (
                    f"Allele frequency > {var_data.thresholds.ba1_benign}: BA1 is met. "
                )
                self.prediction_pm2ba1bs1bs2.BA1 = True
            elif af >= var_data.thresholds.bs1_benign:
                self.comment_pm2ba1bs1bs2 = (
                    f"Allele frequency > {var_data.thresholds.bs1_benign}: BS1 is met. "
                )
                self.prediction_pm2ba1bs1bs2.BS1 = True
            elif af <= var_data.thresholds.pm2_pathogenic:
                self.comment_pm2ba1bs1bs2 = (
                    f"Allele frequency <= {var_data.thresholds.pm2_pathogenic}: PM2 is met. "
                )
                self.prediction_pm2ba1bs1bs2.PM2 = True

            if not self._bs2_not_applicable(var_data) and self._check_zyg(seqvar, var_data):
                self.comment_pm2ba1bs1bs2 += (
                    "The variant is in a recessive, dominant, or X-linked disorder: BS2 is met."
                )
                self.prediction_pm2ba1bs1bs2.BS2 = True

        except AutoAcmgBaseException as e:
            self.comment_pm2ba1bs1bs2 = (
                f"An error occurred while predicting PM2, BA1, BS1, BS2 criteria: {e}"
            )
            self.prediction_pm2ba1bs1bs2 = None

        return self.prediction_pm2ba1bs1bs2, self.comment_pm2ba1bs1bs2

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Predict PM2, BA1, BS1, BS2 criteria."""
        logger.info("Predict PM2, BA1, BS1, BS2")
        pred, comment = self.verify_pm2ba1bs1bs2(seqvar, var_data)
        if pred:
            pm2_pred = (
                AutoACMGPrediction.Applicable
                if pred.PM2
                else (
                    AutoACMGPrediction.NotApplicable
                    if pred.PM2 is False
                    else AutoACMGPrediction.Failed
                )
            )
            ba1_pred = (
                AutoACMGPrediction.Applicable
                if pred.BA1
                else (
                    AutoACMGPrediction.NotApplicable
                    if pred.BA1 is False
                    else AutoACMGPrediction.Failed
                )
            )
            bs1_pred = (
                AutoACMGPrediction.Applicable
                if pred.BS1
                else (
                    AutoACMGPrediction.NotApplicable
                    if pred.BS1 is False
                    else AutoACMGPrediction.Failed
                )
            )
            bs2_pred = (
                AutoACMGPrediction.Applicable
                if pred.BS2
                else (
                    AutoACMGPrediction.NotApplicable
                    if pred.BS2 is False
                    else AutoACMGPrediction.Failed
                )
            )
            pm2_strength = pred.PM2_strength
            ba1_strength = pred.BA1_strength
            bs1_strength = pred.BS1_strength
            bs2_strength = pred.BS2_strength
        else:
            pm2_pred = AutoACMGPrediction.Failed
            ba1_pred = AutoACMGPrediction.Failed
            bs1_pred = AutoACMGPrediction.Failed
            bs2_pred = AutoACMGPrediction.Failed
            pm2_strength = AutoACMGStrength.PathogenicModerate
            ba1_strength = AutoACMGStrength.BenignStandAlone
            bs1_strength = AutoACMGStrength.BenignStrong
            bs2_strength = AutoACMGStrength.BenignStrong
        return (
            AutoACMGCriteria(
                name="PM2",
                prediction=pm2_pred,
                strength=pm2_strength,
                summary=comment,
            ),
            AutoACMGCriteria(
                name="BA1",
                prediction=ba1_pred,
                strength=ba1_strength,
                summary=comment,
            ),
            AutoACMGCriteria(
                name="BS1",
                prediction=bs1_pred,
                strength=bs1_strength,
                summary=comment,
            ),
            AutoACMGCriteria(
                name="BS2",
                prediction=bs2_pred,
                strength=bs2_strength,
                summary=comment,
            ),
        )
