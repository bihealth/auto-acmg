"""
Predictor for Coagulation Factor Deficiency VCEP.
Included genes:
F8 (HGNC:3546),
F9 (HGNC:3551).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN071
https://cspec.genome.network/cspec/ui/svi/doc/GN080
"""

from typing import Dict, List, Optional, Tuple

from loguru import logger

from src.defs.annonars_variant import GnomadExomes
from src.defs.auto_acmg import (
    PM2BA1BS1BS2,
    PS1PM5,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, MissingDataError
from src.defs.seqvar import SeqVar
from src.seqvar.auto_ps1_pm5 import DNA_BASES
from src.seqvar.default_predictor import DefaultSeqVarPredictor

#: VCEP specifications for Coagulation Factor Deficiency.
SPECs: List[VcepSpec] = [
    VcepSpec(
        identifier="GN071",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN080",
        version="1.0.0",
    ),
]

PM1_CLUSTER: Dict[str, Dict[str, Dict[str, List]]] = {
    "HGNC:3546": {  # F8
        "strong": {
            "residues": [
                (391, 392),
                (759, 760),
                (1701, 1705),
                (1708, 1709),
                (1683, 1683),
                (1689, 1689),
                (737, 737),
                (742, 742),
            ]
        },
        "moderate": {
            "residues": [(1667, 1667), (1332, 1332)],  # Residues affecting secretion
            "domains": [(2267, 2304)],  # FXa-binding residues, excluding Ser2283
            "excluded_residues": [(2283, 2283)],  # Excluded residue in FXa-binding region
        },
    },
    "HGNC:3551": {  # F9
        "strong": {
            "residues": [
                (267, 267),
                (315, 315),
                (411, 411),  # Catalytic residues
                (191, 192),
                (226, 227),  # Activation residues
            ]
        },
        "moderate": {"exons": [3, 4, 5]},  # Moderate level for variants within exons 3, 4, or 5
    },
}


class CoagulationFactorDeficiencyPredictor(DefaultSeqVarPredictor):

    def verify_ps1pm5(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PS1PM5], str]:
        """Override PS1/PM5 for Coagulation Factor Deficiency."""
        self.prediction_ps1pm5, self.comment_ps1pm5 = super().verify_ps1pm5(seqvar, var_data)
        if self.prediction_ps1pm5 and self._affect_splicing(var_data):
            self.prediction_ps1pm5.PS1 = False
            self.comment_ps1pm5 = "Variant affects splicing. PS1 is not applicable."
            self.prediction_ps1pm5.PM5 = False
            self.comment_ps1pm5 = "Variant affects splicing. PM5 is not applicable."
        return self.prediction_ps1pm5, self.comment_ps1pm5

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Specify PM1 domains and residues for Coagulation Factor Deficiency."""
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            raise AlgorithmError(
                f"PM1 criteria of CoagulationFactorDeficiencyPredictor VCEP "
                f"is not defined for {var_data.hgnc_id}."
            )

        # Check strong level criteria
        for start, end in gene_cluster.get("strong", {}).get("residues", []):
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant affects a critical residue for {var_data.hgnc_id} "
                    f"between positions {start}-{end}. PM1 is met at the Strong level."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicStrong,
                    summary=comment,
                )

        # Check moderate level criteria for residues
        for start, end in gene_cluster.get("moderate", {}).get("residues", []):
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant affects a moderate-level residue for {var_data.hgnc_id} "
                    f"between positions {start}-{end}. PM1 is met at the Moderate level."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicModerate,
                    summary=comment,
                )

        # Check moderate level criteria for domains, excluding specific residues
        for start, end in gene_cluster.get("moderate", {}).get("domains", []):
            if start <= var_data.prot_pos <= end:
                excluded_residues = gene_cluster.get("moderate", {}).get("excluded_residues", [])
                if not any(start <= var_data.prot_pos <= end for start, end in excluded_residues):
                    comment = (
                        f"Variant falls within a moderate-level region for {var_data.hgnc_id} "
                        f"between positions {start}-{end}. PM1 is met at the Moderate level."
                    )
                    return AutoACMGCriteria(
                        name="PM1",
                        prediction=AutoACMGPrediction.Met,
                        strength=AutoACMGStrength.PathogenicModerate,
                        summary=comment,
                    )

        # Check moderate level criteria for exons in F9
        affected_exon = self._get_affected_exon(var_data, seqvar)
        if var_data.hgnc_id == "HGNC:3551" and affected_exon in [3, 4, 5]:
            comment = (
                f"Variant falls within exon {affected_exon} for {var_data.hgnc_id}. "
                f"PM1 is met at the Moderate level."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=comment,
            )

        # If no criteria match
        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary=f"Variant does not meet the PM1 criteria for {var_data.hgnc_id}.",
        )

    def _absent_in_males(self, var_data: AutoACMGSeqVarData) -> bool:
        """Check if the allele frequency is absent in males."""
        if not var_data.gnomad_exomes:
            raise MissingDataError("No gnomAD exomes data found.")
        ac = self._get_any_af(var_data)
        if ac and ac.bySex and ac.bySex.xx and ac.bySex.xx.ac and ac.bySex.xx.ac == 0:
            return True
        return False

    def verify_pm2ba1bs1bs2(
        self,
        seqvar: SeqVar,
        var_data: AutoACMGSeqVarData,
    ) -> Tuple[Optional[PM2BA1BS1BS2], str]:
        """
        Predicts the PM2, BA1, BS1, BS2 criteria for the sequence variant.

        Note:
            Rules:
            PM2: Absent from males in gnomAD exomes.

            BA1: Allele frequency is greater than 0.033% for F8 and 0.00556% for F9.

            BS1: Allele frequency is greater than 0.0033% for F8 and 0.000556% for F9.

            BS2: No change from the reference implementation.

        Returns:
            BA1BS1BS2PM2: The prediction result.
        """
        self.prediction_pm2ba1bs1bs2 = PM2BA1BS1BS2()
        self.comment_pm2ba1bs1bs2 = ""
        if var_data.hgnc_id == "HGNC:3546":
            var_data.thresholds.ba1_benign = 0.00033
            var_data.thresholds.bs1_benign = 0.000033
        elif var_data.hgnc_id == "HGNC:3551":
            var_data.thresholds.ba1_benign = 0.0000556
            var_data.thresholds.bs1_benign = 0.00000556
        try:
            af = self._get_af(seqvar, var_data)
            if not af:
                self.comment_pm2ba1bs1bs2 = "No allele frequency data found. "
            elif self._ba1_exception(seqvar):
                self.comment_pm2ba1bs1bs2 = "The variant is in the exception list for BA1 criteria."
                self.prediction_pm2ba1bs1bs2.BA1 = False
                self.prediction_pm2ba1bs1bs2.BS1 = False
            elif af >= var_data.thresholds.ba1_benign:
                self.comment_pm2ba1bs1bs2 = "Allele frequency > 5%: BA1 is met. "
                self.prediction_pm2ba1bs1bs2.BA1 = True
            elif af >= var_data.thresholds.bs1_benign:
                self.comment_pm2ba1bs1bs2 = "Allele frequency > 1%: BS1 is met. "
                self.prediction_pm2ba1bs1bs2.BS1 = True
            elif self._absent_in_males(var_data):
                self.comment_pm2ba1bs1bs2 = "The variant is absent in males: PM2 is met. "
                self.prediction_pm2ba1bs1bs2.PM2 = True

            if not self.prediction_pm2ba1bs1bs2.BA1 and af and self._check_zyg(seqvar, var_data):
                self.comment_pm2ba1bs1bs2 += (
                    "The variant is in a recessive, dominant, or X-linked disorder: BS2 is met."
                )
                self.prediction_pm2ba1bs1bs2.BS2 = True

        except AutoAcmgBaseException as e:
            logger.error("Error occurred during PM2, BA1, BS1, BS2 prediction. Error: {}", e)
            self.comment_pm2ba1bs1bs2 = (
                f"An error occurred while predicting PM2, BA1, BS1, BS2 criteria: {e}"
            )
            self.prediction_pm2ba1bs1bs2 = None

        return self.prediction_pm2ba1bs1bs2, self.comment_pm2ba1bs1bs2

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 is not applicable for Coagulation Factor Deficiency."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override PP2, BP1 for Coagulation Factor Deficiency to return not applicable status."""
        return (
            AutoACMGCriteria(
                name="PP2",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary="PP2 is not applicable for the gene.",
            ),
            AutoACMGCriteria(
                name="BP1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary="BP1 is not applicable for the gene.",
            ),
        )

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Change the spliceAI and phyloP threshold for BP7."""
        if var_data.hgnc_id == "HGNC:3546":
            var_data.thresholds.spliceAI_acceptor_gain = 0.05
            var_data.thresholds.spliceAI_acceptor_loss = 0.05
            var_data.thresholds.spliceAI_donor_gain = 0.05
            var_data.thresholds.spliceAI_donor_loss = 0.05
        elif var_data.hgnc_id == "HGNC:3551":
            var_data.thresholds.spliceAI_acceptor_gain = 0.01
            var_data.thresholds.spliceAI_acceptor_loss = 0.01
            var_data.thresholds.spliceAI_donor_gain = 0.01
            var_data.thresholds.spliceAI_donor_loss = 0.01
        var_data.thresholds.phyloP100 = 0.1
        return super().predict_bp7(seqvar, var_data)
