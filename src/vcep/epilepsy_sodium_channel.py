"""
Predictor for Epilepsy Sodium Channel VCEP.
Included genes:
SCN1A (HGNC:10585),
SCN2A (HGNC:10588),
SCN3A (HGNC:10590),
SCN8A (HGNC:10596),
SCN1B (HGNC:10586).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN067
https://cspec.genome.network/cspec/ui/svi/doc/GN068
https://cspec.genome.network/cspec/ui/svi/doc/GN069
https://cspec.genome.network/cspec/ui/svi/doc/GN070
https://cspec.genome.network/cspec/ui/svi/doc/GN076
"""

from typing import List, Tuple

from loguru import logger

from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

#: VCEP specifications for Epilepsy Sodium Channel.
SPECs: List[VcepSpec] = [
    VcepSpec(
        identifier="GN067",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN068",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN069",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN070",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN076",
        version="1.0.0",
    ),
]

PM1_CLUSTER = {
    "HGNC:10585": {  # SCN1A
        "domains": [
            (226, 246),
            (261, 273),
            (451, 464),
            (921, 929),
            (941, 949),
            (951, 964),
            (966, 974),
            (1001, 1009),
            (1051, 1069),
            (1466, 1509),
            (1616, 1624),
            (1626, 1639),
            (1641, 1658),
            (1781, 1799),
            (1806, 1824),
            (1921, 1934),
        ],
    },
    "HGNC:10588": {  # SCN2A
        "domains": [
            (227, 247),
            (262, 272),
            (453, 465),
            (912, 920),
            (932, 940),
            (942, 955),
            (957, 965),
            (992, 1000),
            (1042, 1068),
            (1456, 1489),
            (1606, 1623),
            (1616, 1629),
            (1631, 1649),
            (1771, 1789),
            (1796, 1814),
            (1911, 1924),
        ],
    },
    "HGNC:10590": {  # SCN3A
        "domains": [
            (226, 248),
            (261, 273),
            (452, 465),
            (913, 921),
            (933, 941),
            (943, 956),
            (958, 966),
            (993, 1001),
            (1042, 1052),
            (1451, 1494),
            (1601, 1619),
            (1611, 1624),
            (1626, 1635),
            (1766, 1784),
            (1791, 1809),
            (1906, 1919),
        ],
    },
    "HGNC:10596": {  # SCN8A
        "domains": [
            (230, 244),
            (265, 273),
            (439, 452),
            (906, 914),
            (926, 934),
            (936, 949),
            (951, 959),
            (986, 994),
            (1034, 1053),
            (1447, 1498),
            (1597, 1618),
            (1607, 1620),
            (1607, 1628),
            (1761, 1788),
            (1786, 1804),
            (1901, 1914),
        ],
    },
}


class EpilepsySodiumChannelPredictor(DefaultSeqVarPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override predict_pm1 to iclude domains for Epilepsy Sodium Channel genes."""
        logger.info("Predict PM1")

        # Check if SCN1B, where PM1 is not applicable
        if var_data.hgnc_id == "HGNC:10586":  # SCN1B
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicModerate,
                summary="PM1 is not applicable for SCN1B.",
            )

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return super().predict_pm1(seqvar, var_data)

        # Check if the variant falls within any of the critical residues
        for start_aa, end_aa in gene_cluster.get("domains", []):
            if start_aa <= var_data.prot_pos <= end_aa:
                comment = (
                    f"Variant falls within a critical residue region for {var_data.hgnc_id} "
                    f"between positions {start_aa}-{end_aa}. PM1 is met."
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

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.000001  # Very rare
        if var_data.hgnc_id == "HGNC:10585":  # SCN1A
            var_data.thresholds.ba1_benign = 0.0002
            var_data.thresholds.bs1_benign = 0.000004
        elif var_data.hgnc_id in ["HGNC:10588", "HGNC:10590", "HGNC:10596"]:  # SCN2A, SCN3A, SCN8A
            var_data.thresholds.ba1_benign = 0.0001
            var_data.thresholds.bs1_benign = 0.000002
        elif var_data.hgnc_id == "HGNC:10586":  # SCN1B
            var_data.thresholds.ba1_benign = 0.003
            var_data.thresholds.bs1_benign = 0.0001
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override PP2 and BP1 for Epilepsy Sodium Channel to return not applicable status."""
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

    def predict_pp3bp4(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Use REVEL for PP3 and BP4 for Epilepsy Sodium Channel."""
        var_data.thresholds.pp3bp4_strategy = "revel"
        return super().predict_pp3bp4(seqvar, var_data)
