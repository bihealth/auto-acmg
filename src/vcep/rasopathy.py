"""
Predictor for RASopathy VCEP.
Included genes:
SHOC2 (HGNC:15454),
NRAS (HGNC:7989),
RAF1 (HGNC:9829),
SOS1 (HGNC:11187),
SOS2 (HGNC:11188),
PTPN11 (HGNC:9644),
KRAS (HGNC:6407),
MAP2K1 (HGNC:6840),
HRAS (HGNC:5173),
RIT1 (HGNC:10023),
MAP2K2 (HGNC:6842),
BRAF (HGNC:1097),
MRAS (HGNC:7227),
LZTR1 (HGNC:6742),
RRAS2 (HGNC:17271),
PPP1CB (HGNC:9282).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN038
https://cspec.genome.network/cspec/ui/svi/doc/GN039
https://cspec.genome.network/cspec/ui/svi/doc/GN040
https://cspec.genome.network/cspec/ui/svi/doc/GN041
https://cspec.genome.network/cspec/ui/svi/doc/GN042
https://cspec.genome.network/cspec/ui/svi/doc/GN043
https://cspec.genome.network/cspec/ui/svi/doc/GN044
https://cspec.genome.network/cspec/ui/svi/doc/GN045
https://cspec.genome.network/cspec/ui/svi/doc/GN046
https://cspec.genome.network/cspec/ui/svi/doc/GN047
https://cspec.genome.network/cspec/ui/svi/doc/GN048
https://cspec.genome.network/cspec/ui/svi/doc/GN049
https://cspec.genome.network/cspec/ui/svi/doc/GN087
https://cspec.genome.network/cspec/ui/svi/doc/GN094
https://cspec.genome.network/cspec/ui/svi/doc/GN127
https://cspec.genome.network/cspec/ui/svi/doc/GN128
"""

from typing import Dict, List, Optional, Tuple, Union

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import (
    PM2BA1BS1BS2,
    AutoACMGCriteria,
    AutoACMGData,
    AutoACMGPrediction,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.exceptions import AutoAcmgBaseException
from src.defs.seqvar import SeqVar

#: VCEP specifications for RASopathy.
SPECs: List[VcepSpec] = [
    VcepSpec(
        identifier="GN038",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN039",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN040",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN041",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN042",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN043",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN044",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN045",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN046",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN047",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN048",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN049",
        version="2.0.0",
    ),
    VcepSpec(
        identifier="GN087",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN094",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN127",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN128",
        version="1.0.0",
    ),
]

# fmt: off
PM1_CLUSTER_RASOPATHY: Dict[str, Dict[str, List[Union[int, Tuple[int, int]]]]] = {
    "HGNC:7989": {  # NRAS
        "domains": [(10, 17), (25, 40), (57, 64), (145, 156)],  # P-loop, SW1, SW2, SAK
    },
    "HGNC:9829": {  # RAF1
        "domains": [(251, 266)],  # CR2 domain in exon 7
        "exons": [14, 17],  # Critical exons
    },
    "HGNC:11187": {  # SOS1
        "domains": [(420, 500)],  # PH domain
    },
    "HGNC:11188": {  # SOS2
        "domains": [(418, 498)],  # PH domain
    },
    "HGNC:9644": {  # PTPN11
        "domains": [(4, 4), (7, 9), (58, 63), (69, 77), (247, 247), (251, 251),
                    (255, 256), (258, 258), (261, 261), (265, 265),
                    (278, 281), (284, 284)],  # Critical N-SH2 and PTPN domain residues
    },
    "HGNC:6407": {  # KRAS
        "domains": [(10, 17), (25, 40), (57, 64), (145, 156)],  # P-loop, SW1, SW2, SAK
    },
    "HGNC:6840": {  # MAP2K1
        "domains": [(43, 61), (124, 134)],  # Specific domains
    },
    "HGNC:5173": {  # HRAS
        "domains": [(10, 17), (25, 40), (57, 64), (145, 156)],  # P-loop, SW1, SW2, SAK
    },
    "HGNC:10023": {  # RIT1
        "domains": [(28, 35), (43, 58), (75, 82)],  # P-loop, SW1, SW2
    },
    "HGNC:6842": {  # MAP2K2
        "domains": [(47, 65), (128, 138)],  # Specific domains
    },
    "HGNC:1097": {  # BRAF
        "domains": [(459, 474), (594, 627)],  # P-loop, CR3 activation segment
        "exons": [6, 11],  # Critical exons
    },
    "HGNC:7227": {  # MRAS
        "domains": [(20, 27), (35, 50), (67, 74)],  # P-loop, SW1, SW2
    },
    "HGNC:17271": {  # RRAS2
        "domains": [(21, 28), (36, 51), (68, 75)],  # P-loop, SW1, SW2
    },
}


class RASopathyPredictor(DefaultPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGData) -> AutoACMGCriteria:
        """Override predict_pm1 to specify critical domains for RASopathy."""
        logger.info("Predict PM1")

        if var_data.hgnc_id in [
            "HGNC:15454",  # SHOC2
            "HGNC:6742",  # LZTR1
            "HGNC:9282",  # PPP1CB
        ]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=f"PM1 is not applicable for {var_data.hgnc_id}.",
            )

        gene_info = PM1_CLUSTER_RASOPATHY.get(var_data.hgnc_id, None)
        if not gene_info:
            return super().predict_pm1(seqvar, var_data)

        # Check if the variant falls within a critical region
        for start, end in gene_info.get("domains", []):   # type: ignore
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

        # Check if the variant is in a critical exon
        affected_exon = self._get_affected_exon(var_data, seqvar)
        if affected_exon in gene_info.get("exons", []):
            comment = (
                f"Variant affects a critical exon {affected_exon} in {var_data.hgnc_id}. PM1 is met."
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
        var_data.thresholds.ba1_benign = 0.0005
        var_data.thresholds.bs1_benign = 0.00025
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

    def predict_pp2bp1(self, seqvar: SeqVar, var_data: AutoACMGData) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """
        Override PP2 and BP1 for RASopathy. PP2 is met for missense changes. BP1 is not applicable.
        """
        pp2 = False
        comment = "Not applicable for the gene."
        if var_data.hgnc_id in ["HGNC:9644", "HGNC:6840", "HGNC:1097", "HGNC:9282"]:
            if self._is_missense(var_data):
                pp2 = True
                comment = f"PP2 is met for {var_data.hgnc_id} as the variant is a missense change."
            else:
                pp2 = False
                comment = (
                    f"PP2 is not met for {var_data.hgnc_id} as the variant is not a missense "
                    "change."
                )
        return (
            AutoACMGCriteria(
                name="PP2",
                prediction=AutoACMGPrediction.Met if pp2 else AutoACMGPrediction.NotMet,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=comment,
            ),
            AutoACMGCriteria(
                name="BP1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.BenignSupporting,
                summary="BP1 is not applicable for the gene.",
            ),
        )


    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGData) -> bool:
        """BP3 is not applicable for RASopathy."""
        return True
