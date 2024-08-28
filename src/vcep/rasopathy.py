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

from typing import Dict, List, Tuple, Union

from loguru import logger

from src.criteria.default_predictor import DefaultPredictor
from src.defs.auto_acmg import AutoACMGCriteria, AutoACMGData, AutoACMGPrediction, AutoACMGStrength
from src.defs.seqvar import SeqVar

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
        """
        Override predict_pm1 to include VCEP-specific logic for RASopathy.
        """
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

    def predict_pp2bp1(self, seqvar: SeqVar, var_data: AutoACMGData) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override PP2 and BP1 for RASopathy."""
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
        """Override BP3 for RASopathy."""
        return True
