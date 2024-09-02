"""
Predictor for Rett and Angelman-like Disorders VCEP.
Included genes:
TCF4 (HGNC:11634),
SLC9A6 (HGNC:11079),
CDKL5 (HGNC:11411),
FOXG1 (HGNC:3811),
MECP2 (HGNC:6990),
UBE3A (HGNC:12496).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN032
https://cspec.genome.network/cspec/ui/svi/doc/GN033
https://cspec.genome.network/cspec/ui/svi/doc/GN034
https://cspec.genome.network/cspec/ui/svi/doc/GN035
https://cspec.genome.network/cspec/ui/svi/doc/GN036
https://cspec.genome.network/cspec/ui/svi/doc/GN037
"""

from typing import Dict, List, Optional, Tuple, Union

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import (
    PM2BA1BS1BS2,
    PM4BP3,
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.exceptions import AutoAcmgBaseException
from src.defs.seqvar import SeqVar

#: VCEP specifications for Rett and Angelman-like Disorders.
SPECs: List[VcepSpec] = [
    VcepSpec(
        identifier="GN032",
        version="3.0.0",
    ),
    VcepSpec(
        identifier="GN033",
        version="3.0.0",
    ),
    VcepSpec(
        identifier="GN034",
        version="3.0.0",
    ),
    VcepSpec(
        identifier="GN035",
        version="3.0.0",
    ),
    VcepSpec(
        identifier="GN036",
        version="3.0.0",
    ),
    VcepSpec(
        identifier="GN037",
        version="4.0.0",
    ),
]

# fmt: off
PM1_CLUSTER: Dict[str, Dict[str, List[Tuple[int, int]]]] = {
    "HGNC:11634": {  # TCF4
        "domains": [(564, 617)],  # Basic Helix-Loop-Helix (bHLH) domain
    },
    "HGNC:11411": {  # CDKL5
        "domains": [(19, 43), (169, 171)],  # ATP binding region, TEY phosphorylation site
    },
    "HGNC:3811": {  # FOXG1
        "domains": [(181, 275)],  # Forkhead domain
    },
    "HGNC:6990": {  # MECP2
        "domains": [(90, 162), (302, 306)],  # Methyl-DNA binding (MBD), Transcriptional repression domain (TRD)
    },
    "HGNC:12496": {  # UBE3A
        "domains": [(820, 820)],  # 3’ cysteine binding site
    },
}

#: Domains to exclude from PM4
PM4_EXCLUDE: Dict[str, List[Tuple[int, int]]] = {
    "HGNC:11411": [(904, int(1e6))],  # Exclude C-terminal region
    "HGNC:3811": [
        (35, 57),  # Histine-rich region
        (58, 86),  # Prolinne- and Glutamate-rich region
        (105, 112)  # Proline-rich region
    ],
    "HGNC:6990": [(381, 405)],   # Proline-rich region
}

#: FOXG1 BP3 region
FOXG1_BP3_REGION: List[Tuple[int, int]] = [
    (47, 57),  # Poly-His region
    (70, 73),  # Poly-Glutamin region
    (58, 61),  # Poly-Proline region
    (65, 69),  # Poly-Proline region
    (74, 80),  # Poly-Proline region
]


class RettAngelmanPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override predict_pm1 to specify critical domains for Rett and Angelman-like Disorders."""
        logger.info("Predict PM1")

        # PM1 is not applicable for SLC9A6
        if var_data.hgnc_id == "HGNC:11079":
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=f"PM1 is not applicable for SLC9A6.",
            )

        gene_info = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_info:
            return super().predict_pm1(seqvar, var_data)

        # Check if the variant falls within a critical region
        for start, end in gene_info.get("domains", []):
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant affects a critical residue in {var_data.hgnc_id} "
                    f"within positions {start}-{end}. PM1 is met."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicModerate,
                    summary=comment,
                )

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary=f"Variant does not meet the PM1 criteria for {var_data.hgnc_id}.",
        )

    def verify_pm2ba1bs1bs2(
        self,
        seqvar: SeqVar,
        var_data: AutoACMGData,
    ) -> Tuple[Optional[PM2BA1BS1BS2], str]:
        """
        Predicts the PM2, BA1, BS1, BS2 criteria for the sequence variant.

        Note:
            Rules:
            PM2: Absent from controls allele frequency data.

            BA1: Allele frequency is greater than 0.05%.

            BS1: Allele frequency is between 0.025% and 0.05%.

            BS2: No change from the reference implementation.

        Returns:
            BA1BS1BS2PM2: The prediction result.
        """
        self.prediction_pm2ba1bs1bs2 = PM2BA1BS1BS2()
        self.comment_pm2ba1bs1bs2 = ""
        var_data.thresholds.ba1_benign = 0.0003
        var_data.thresholds.bs1_benign = 0.00008
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
            elif not self._get_control_af(var_data):
                self.comment_pm2ba1bs1bs2 = "No control allele frequency data found. "
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

    def _exclude_pm4(self, seqvar: SeqVar, var_data: AutoACMGData) -> bool:
        """Check if the variant should be excluded from PM4."""
        if var_data.hgnc_id in PM4_EXCLUDE:
            for start, end in PM4_EXCLUDE[var_data.hgnc_id]:
                if start <= var_data.prot_pos <= end:
                    return True
        return False

    def _in_foxg1_bp3_region(self, var_data: AutoACMGData) -> bool:
        """Check if the variant is in the BP3 region for FOXG1."""
        if var_data.hgnc_id != "HGNC:3811":
            return False
        for start, end in FOXG1_BP3_REGION:
            if start <= var_data.prot_pos <= end:
                return True
        return False

    def verify_pm4bp3(self, seqvar: SeqVar, var_data: AutoACMGData) -> Tuple[Optional[PM4BP3], str]:
        """
        Override PM4 and BP3 for Rett and Angelman-like Disorders. PM4 is met for stop-loss variants
        and in-frame deletions/insertions that are not in repeat regions or conserved domains.
        BP3 is met for in-frame deletions/insertions in the BP3 region for FOXG1.
        """
        self.prediction_pm4bp3 = PM4BP3()
        self.comment_pm4bp3 = ""
        try:
            # Stop-loss variants are considered as PM4
            if self._is_stop_loss(var_data):
                self.comment_pm4bp3 = "Variant consequence is stop-loss. PM4 is met."
                self.prediction_pm4bp3.PM4 = True
                self.prediction_pm4bp3.BP3 = False
            # In-frame deletions/insertions
            elif self.is_inframe_delins(var_data):
                self.comment_pm4bp3 = f"Variant consequence is in-frame deletion/insertion. "
                if not self._in_repeat_region(seqvar) and not self._exclude_pm4(seqvar, var_data):
                    self.comment_pm4bp3 += (
                        "Variant is not in a repeat region or a conserved domain. PM4 is met."
                    )
                    self.prediction_pm4bp3.PM4 = True
                    self.prediction_pm4bp3.BP3 = False
                else:
                    self.comment_pm4bp3 += (
                        "Variant is in a repeat region or not in a conserved domain or in excluded "
                        "region."
                    )
                    self.prediction_pm4bp3.PM4 = False
                    self.prediction_pm4bp3.BP3 = False
                # BP3 for FOXG1
                if self._in_foxg1_bp3_region(var_data):
                    self.comment_pm4bp3 += " Variant is in the BP3 region for FOXG1."
                    self.prediction_pm4bp3.BP3 = True
            else:
                self.comment_pm4bp3 = (
                    "Variant consequence is not indel or stop-loss. PM4 and BP3 are not met."
                )
                self.prediction_pm4bp3.PM4 = False
                self.prediction_pm4bp3.BP3 = False

        except AutoAcmgBaseException as e:
            logger.error("Failed to predict PM4 and BP3 criteria. Error: {}", e)
            self.comment_pm4bp3 = f"An error occured while predicting PM4 and BP3 criteria: {e}"
            self.prediction_pm4bp3 = None
        return self.prediction_pm4bp3, self.comment_pm4bp3

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """
        Override PP2 and BP1 to return not applicable status for Rett and Angelman-like Disorders.
        """
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

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Change BP7 thresholds for Rett and Angelman-like Disorders VCEP."""
        var_data.thresholds.phyloP100 = 0.1
        return super().predict_bp7(seqvar, var_data)
