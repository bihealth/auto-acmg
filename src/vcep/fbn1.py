"""
Predictor for FBN1 VCEP.
Included gene: FBN1 (HGNC:3603).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN022
"""

from typing import Optional, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    PS1PM5,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

#: VCEP specification for FBN1.
SPEC: VcepSpec = VcepSpec(
    identifier="GN022",
    version="1.0.0",
)

# Define the critical residues for PM1 strong consideration
# fmt: off
PM1_CLUSTER = {
    "HGNC:3603": {
        "strong_residues": [
            250, 257, 262, 271, 273, 286, 292, 299, 304, 313, 315, 328, 494, 499, 504, 513, 515, 528,
            534, 541, 546, 555, 557, 570, 576, 582, 587, 596, 598, 611, 617, 623, 628, 637, 639, 652,
            727, 734, 739, 748, 750, 763, 769, 776, 781, 790, 792, 805, 811, 816, 821, 830, 832, 845,
            914, 921, 926, 935, 937, 950, 1032, 1039, 1044, 1053, 1055, 1068, 1074, 1081, 1086, 1095,
            1097, 1111, 1117, 1124, 1129, 1138, 1140, 1153, 1159, 1166, 1171, 1180, 1182, 1195, 1201,
            1208, 1212, 1221, 1223, 1236, 1242, 1249, 1254, 1263, 1265, 1278, 1284, 1291, 1296, 1305,
            1307, 1320, 1326, 1333, 1339, 1348, 1350, 1361, 1367, 1374, 1380, 1389, 1391, 1402, 1408,
            1415, 1420, 1429, 1431, 1444, 1450, 1456, 1461, 1470, 1472, 1485, 1491, 1497, 1502, 1511,
            1513, 1526, 1610, 1617, 1622, 1631, 1633, 1646, 1652, 1658, 1663, 1672, 1674, 1687, 1770,
            1777, 1782, 1791, 1793, 1806, 1812, 1818, 1824, 1833, 1835, 1847, 1853, 1860, 1865, 1874,
            1876, 1889, 1895, 1900, 1905, 1914, 1916, 1928, 1934, 1942, 1947, 1956, 1958, 1971, 1977,
            1984, 1989, 1998, 2000, 2011, 2017, 2024, 2029, 2038, 2040, 2053, 2131, 2137, 2142, 2151,
            2153, 2164, 2170, 2176, 2181, 2190, 2192, 2204, 2210, 2217, 2221, 2230, 2232, 2245, 2251,
            2258, 2265, 2274, 2276, 2289, 2295, 2302, 2307, 2316, 2318, 2331, 2406, 2413, 2418, 2427,
            2429, 2442, 2448, 2455, 2459, 2468, 2470, 2483, 2489, 2496, 2500, 2509, 2511, 2522, 2528,
            2535, 2541, 2550, 2552, 2565, 2571, 2577, 2581, 2590, 2592, 2605, 2611, 2617, 2622, 2631,
            2633, 2646, 2652, 2659, 2663, 2672, 2674, 2686,
        ],
        "moderate_residues": [
            85, 89, 94, 100, 102, 111, 119, 123, 129, 134, 136, 145, 150, 154, 160, 166, 168, 177,
            186, 195, 204, 209, 210, 221, 224, 231, 244, 453, 460, 465, 474, 476, 488, 336, 345, 358,
            359, 360, 365, 377, 389, 661, 670, 683, 684, 685, 696, 699, 711, 853, 862, 875, 876, 887,
            890, 896, 958, 967, 980, 981, 982, 993, 996, 1008, 1534, 1549, 1562, 1563, 1564, 1574,
            1577, 1589, 1695, 1706, 1719, 1720, 1721, 1733, 1736, 1748, 2061, 2070, 2083, 2084, 2085,
            2096, 2099, 2111, 2339, 2348, 2363, 2364, 2365, 2375, 2378, 2390,
        ]
    }
}


class FBN1Predictor(DefaultSeqVarPredictor):

    def verify_ps1pm5(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PS1PM5], str]:
        """Override PS1/PM5 for FBN1."""
        self.prediction_ps1pm5, self.comment_ps1pm5 = super().verify_ps1pm5(seqvar, var_data)
        if self.prediction_ps1pm5 and self._is_missense(var_data) and self._affect_splicing(var_data):
            self.prediction_ps1pm5.PS1 = False
            self.comment_ps1pm5 = "Variant affects splicing. PS1 is not applicable."
            self.prediction_ps1pm5.PM5 = False
            self.comment_ps1pm5 = "Variant affects splicing. PM5 is not applicable."
        return self.prediction_ps1pm5, self.comment_ps1pm5

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override predict_pm1 to specify critical residues for FBN1."""
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check strong level criteria
        if var_data.prot_pos in gene_cluster["strong_residues"]:
            comment = (
                f"Variant affects a critical cysteine residue in FBN1 at position "
                f"{var_data.prot_pos}, leading to disulfide bond folding defects. PM1 is met at the Strong level."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicStrong,
                summary=comment,
            )

        # Check moderate level criteria
        if var_data.prot_pos in gene_cluster["moderate_residues"]:
            comment = (
                f"Variant affects a residue in FBN1 at position {var_data.prot_pos}, which can "
                f"lead to disulfide bond folding defects. PM1 is met at the Moderate level."
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
            summary=f"Variant does not meet the PM1 criteria for FBN1.",
        )

    def _bs2_not_applicable(self, var_data: AutoACMGSeqVarData) -> bool:
        """BS2 is not applicable for FBN1."""
        return True

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.000005
        var_data.thresholds.ba1_benign = 0.001
        var_data.thresholds.bs1_benign = 0.00005
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 is not applicable for FBN1."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """
        Override predict_pp2bp1 to include VCEP-specific logic for FBN1. PP2 is not changed, but
        BP1 is not applicable.
        """
        logger.info("Predict PP2 and BP1")
        pred, comment = self.verify_pp2bp1(seqvar, var_data)
        if pred:
            pp2_pred = (
                AutoACMGPrediction.Met
                if pred.PP2
                else (AutoACMGPrediction.NotMet if pred.PP2 is False else AutoACMGPrediction.Failed)
            )
            pp2_strength = pred.PP2_strength
        else:
            pp2_pred = AutoACMGPrediction.Failed
            pp2_strength = AutoACMGStrength.PathogenicSupporting
        return (
            AutoACMGCriteria(
                name="PP2",
                prediction=pp2_pred,
                strength=pp2_strength,
                summary=comment,
            ),
            AutoACMGCriteria(
                name="BP1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary="BP1 is not applicable for the gene.",
            ),
        )
