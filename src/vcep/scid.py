"""
Predictor for Severe Combined Immunodeficiency Disease VCEP.
Included genes:
FOXN1 (HGNC:12765),
ADA (HGNC:186),
DCLRE1C (HGNC:17642),
IL7R (HGNC:6024),
JAK3 (HGNC:6193),
RAG1 (HGNC:9831),
RAG2 (HGNC:9832),
IL2RG (HGNC:6010).
Links:
https://cspec.genome.network/cspec/ui/svi/doc/GN113
https://cspec.genome.network/cspec/ui/svi/doc/GN114
https://cspec.genome.network/cspec/ui/svi/doc/GN116
https://cspec.genome.network/cspec/ui/svi/doc/GN119
https://cspec.genome.network/cspec/ui/svi/doc/GN121
https://cspec.genome.network/cspec/ui/svi/doc/GN123
https://cspec.genome.network/cspec/ui/svi/doc/GN124
https://cspec.genome.network/cspec/ui/svi/doc/GN129
"""

from typing import Dict, List, Optional, Tuple, Union

from loguru import logger

from src.defs.auto_acmg import (
    PP3BP4,
    AlleleCondition,
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGSeqVarData,
    AutoACMGStrength,
    VcepSpec,
)
from src.defs.exceptions import AutoAcmgBaseException, MissingDataError
from src.defs.seqvar import SeqVar
from src.seqvar.default_predictor import DefaultSeqVarPredictor

#: VCEP specifications for Severe Combined Immunodeficiency Disease.
SPECs: List[VcepSpec] = [
    VcepSpec(
        identifier="GN113",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN114",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN116",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN119",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN121",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN123",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN124",
        version="1.0.0",
    ),
    VcepSpec(
        identifier="GN129",
        version="1.0.0",
    ),
]

# fmt: off
PM1_CLUSTER_SCID: Dict[str, Dict[str, List[Union[int, Tuple[int, int]]]]] = {
    "HGNC:12765": {  # FOXN1
        "domains": [(270, 367)],  # DNA binding forkhead domain
    },
    "HGNC:6193": {  # JAK3
        "domains": [(651, 759)],  # JH2 domain residues R651W and C759R
    },
    "HGNC:9831": {  # RAG1
        "domains": [(394, 460), (461, 517)],  # NBD domain, DDBD domain
        "supporting_domains": [(387, 1011)],  # Core domain (supporting strength)
    },
    "HGNC:9832": {  # RAG2
        "domains": [(414, 487)],  # PHD domain
        "supporting_domains": [(1, 383)],  # Core domain (supporting strength)
    },
    "HGNC:6010": {  # IL2RG
        "strong_residues": [
            62, 72, 102, 115,  # Conserved cysteine residues
            224, 226, 691, 285,  # CpG dinucleotides
            237, 238, 239, 240, 241,  # WSxWS motif
        ],
        "strong_domains": [(263, 283)],  # Transmembrane domain
    },
}


class SCIDPredictor(DefaultSeqVarPredictor):

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """
        Override predict_pm1 to specify critical domains for Severe Combined Immunodeficiency
        Disease.
        """
        logger.info("Predict PM1")

        if var_data.hgnc_id in [
            "HGNC:186",  # ADA
            "HGNC:17642",  # DCLRE1C
            "HGNC:6024",  # IL7R
        ]:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotApplicable,
                strength=AutoACMGStrength.PathogenicSupporting,
                summary=f"PM1 is not applicable for {var_data.hgnc_id}.",
            )

        gene_info = PM1_CLUSTER_SCID.get(var_data.hgnc_id, None)
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

        # Check if the variant falls within a supporting region
        for start, end in gene_info.get("supporting_domains", []):   # type: ignore
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant affects a residue in the supporting region of {var_data.hgnc_id} "
                    f"within positions {start}-{end}. PM1 is met at a supporting level."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicSupporting,
                    summary=comment,
                )

        # Check if the variant affects a strong residue
        if var_data.prot_pos in gene_info.get("strong_residues", []):
            comment = (
                f"Variant affects a strong residue in {var_data.hgnc_id} "
                f"at position {var_data.prot_pos}. PM1 is met."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicStrong,
                summary=comment,
            )

        # Check if the variant falls within a strong domain
        for start, end in gene_info.get("strong_domains", []):   # type: ignore
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant affects a residue in the strong domain of {var_data.hgnc_id} "
                    f"within positions {start}-{end}. PM1 is met at a strong level."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicStrong,
                    summary=comment,
                )

        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary=f"Variant does not meet the PM1 criteria for {var_data.hgnc_id}.",
        )

    def _check_zyg(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """
        Check the zygosity of the sequence variant.

        BS2 only to be used when variant is observed in the homozygous state in a healthy adult.
        For FOXN1, the variant is not in a recessive disorder.

        Args:
            variant_data: The variant data.

        Returns:
            True if the variant is recessive (homozygous), dominant (heterozygous), or X-linked
            (hemizygous) disorder.
        """
        if var_data.hgnc_id == "HGNC:12765":  # FOXN1
            return False
        allele_condition = self._get_allele_cond(seqvar)
        self.comment_pm2ba1bs1bs2 += f"Allele condition: {allele_condition.name}.\n"
        controls_af = self._get_control_af(var_data)
        any_af = self._get_any_af(var_data)
        af = controls_af or any_af
        if not af or not af.bySex:
            self.comment_pm2ba1bs1bs2 += "No controls allele data found in control data.\n"
            raise MissingDataError("No raw data found in control data.")

        if not af.bySex.overall:
            self.comment_pm2ba1bs1bs2 += "No allele data found for overall in control data.\n"
            raise MissingDataError("No allele data found for overall in control data.")
        ac = af.bySex.overall.ac if af.bySex.overall.ac else 0
        nhomalt = af.bySex.overall.nhomalt if af.bySex.overall.nhomalt else 0
        self.comment_pm2ba1bs1bs2 += f"Allele count: {ac}, Nhomalt: {nhomalt}.\n"
        if allele_condition == AlleleCondition.Recessive:
            if nhomalt > 5:
                self.comment_pm2ba1bs1bs2 += (
                    f"Nhomalt {nhomalt} > 5.\n"
                    "The variant is in a recessive (homozygous) disorder."
                )
                return True
        return False

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        if var_data.hgnc_id == "HGNC:12765":
            var_data.thresholds.pm2_pathogenic = 0.00002412
            var_data.thresholds.ba1_benign = 0.00447
            var_data.thresholds.bs1_benign = 0.00141
        elif var_data.hgnc_id == "HGNC:186":
            var_data.thresholds.pm2_pathogenic = 0.0001742
            var_data.thresholds.ba1_benign = 0.00721
            var_data.thresholds.bs1_benign = 0.00161
        elif var_data.hgnc_id == "HGNC:17642":
            var_data.thresholds.pm2_pathogenic = 0.00003266
            var_data.thresholds.ba1_benign = 0.00346
            var_data.thresholds.bs1_benign = 0.00078
        elif var_data.hgnc_id == "HGNC:6024":
            var_data.thresholds.pm2_pathogenic = 0.00004129
            var_data.thresholds.ba1_benign = 0.00566
            var_data.thresholds.bs1_benign = 0.00126
        elif var_data.hgnc_id == "HGNC:6193":
            var_data.thresholds.pm2_pathogenic = 0.000115
            var_data.thresholds.ba1_benign = 0.00447
            var_data.thresholds.bs1_benign = 0.001
        elif var_data.hgnc_id == "HGNC:9831":
            var_data.thresholds.pm2_pathogenic = 0.000102
            var_data.thresholds.ba1_benign = 0.00872
            var_data.thresholds.bs1_benign = 0.00195
        elif var_data.hgnc_id == "HGNC:9832":
            var_data.thresholds.pm2_pathogenic = 0.0000588
            var_data.thresholds.ba1_benign = 0.00872
            var_data.thresholds.bs1_benign = 0.00195
        elif var_data.hgnc_id == "HGNC:6010":
            var_data.thresholds.pm2_pathogenic = 0.000124
            var_data.thresholds.ba1_benign = 0.01110
            var_data.thresholds.bs1_benign = 0.00249
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _is_conserved(self, var_data: AutoACMGSeqVarData) -> bool:
        """
        Override the default _is_conserved method to ignore this check for SCID genes.
        """
        if var_data.hgnc_id == "HGNC:12765":
            return super()._is_conserved(var_data)
        return False  # Ignore conservation check for SCID genes

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """
        Override PP2 and BP1 for Severe Combined Immunodeficiency Disease to return not applicable
        status.
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

    def verify_pp3bp4(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PP3BP4], str]:
        """Predict PP3 and BP4 criteria."""
        self.prediction_pp3bp4 = PP3BP4()
        self.comment_pp3bp4 = ""
        try:
            if var_data.hgnc_id == "HGNC:12765":
                score = "revel"
                var_data.thresholds.revel_pathogenic = 0.644
                var_data.thresholds.revel_benign = 0.29
                self.prediction_pp3bp4.PP3 = self._is_pathogenic_score(
                    var_data,
                    (score, getattr(var_data.thresholds, f"{score}_pathogenic")),
                )
                self.prediction_pp3bp4.BP4 = self._is_benign_score(
                    var_data,
                    (score, getattr(var_data.thresholds, f"{score}_benign")),
                )

                var_data.thresholds.spliceAI_acceptor_gain = 0.2
                var_data.thresholds.spliceAI_acceptor_loss = 0.2
                var_data.thresholds.spliceAI_donor_gain = 0.2
                var_data.thresholds.spliceAI_donor_loss = 0.2
                self.prediction_pp3bp4.PP3 = self.prediction_pp3bp4.PP3 or self._affect_spliceAI(
                    var_data
                )
                var_data.thresholds.spliceAI_acceptor_gain = 0.1
                var_data.thresholds.spliceAI_acceptor_loss = 0.1
                var_data.thresholds.spliceAI_donor_gain = 0.1
                var_data.thresholds.spliceAI_donor_loss = 0.1
                self.prediction_pp3bp4.BP4 = self.prediction_pp3bp4.BP4 and not self._affect_spliceAI(
                    var_data
                )
            else:
                var_data.thresholds.spliceAI_acceptor_gain = 0.2
                var_data.thresholds.spliceAI_acceptor_loss = 0.2
                var_data.thresholds.spliceAI_donor_gain = 0.2
                var_data.thresholds.spliceAI_donor_loss = 0.2
                self.prediction_pp3bp4.PP3 = self.prediction_pp3bp4.PP3 or self._affect_spliceAI(
                    var_data
                )
                # Not applicable BP4
                self.prediction_pp3bp4.BP4 = False

        except AutoAcmgBaseException as e:
            self.comment_pp3bp4 = f"An error occurred during prediction. Error: {e}"
            self.prediction_pp3bp4 = None
        return self.prediction_pp3bp4, self.comment_pp3bp4

    def predict_bp7(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """
        Override donor and acceptor positions for Severe Combined Immunodeficiency Disease genes
        VCEP.
        """
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
