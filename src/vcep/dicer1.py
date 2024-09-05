"""
Predictor for DICER1 and miRNA-Processing Gene VCEP.
Included gene: DICER1 (HGNC:17098).
Link: https://cspec.genome.network/cspec/ui/svi/doc/GN024
"""

from typing import Dict, List, Optional, Tuple

from loguru import logger

from src.defs.annonars_variant import VariantResult
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

#: VCEP specification for DICER1.
SPEC: VcepSpec = VcepSpec(
    identifier="GN024",
    version="1.3.0",
)


PM1_CLUSTER: Dict[str, Dict[str, Dict[str, List]]] = {
    # DICER1
    "HGNC:17098": {
        "moderate": {
            "residues": [1344, 1705, 1709, 1713, 1809, 1810, 1813],  # Metal ion-binding residues
        },
        "supporting": {
            "domains": [(1682, 1846)],  # RNase IIIb domain excluding the metal ion-binding residues
        },
    },
}


class DICER1Predictor(DefaultSeqVarPredictor):

    def _is_pathogenic(self, variant_info: VariantResult) -> bool:
        """Override _is_pathogenic for DICER1."""
        if variant_info.clinvar and variant_info.clinvar.records:
            for rec in variant_info.clinvar.records:
                if (
                    (r := rec.classifications)
                    and (g := r.germlineClassification)
                    and g.description in ["Pathogenic"]
                ):
                    return True
        return False

    def verify_ps1pm5(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[Optional[PS1PM5], str]:
        """Override PS1/PM5 for DICER1."""
        self.prediction_ps1pm5, self.comment_ps1pm5 = super().verify_ps1pm5(seqvar, var_data)
        if (
            self.prediction_ps1pm5
            and self._is_missense(var_data)
            and self._affect_splicing(var_data)
        ):
            self.prediction_ps1pm5.PS1 = False
            self.comment_ps1pm5 = "Variant affects splicing. PS1 is not applicable."
            self.prediction_ps1pm5.PM5 = False
            self.comment_ps1pm5 = "Variant affects splicing. PM5 is not applicable."
        return self.prediction_ps1pm5, self.comment_ps1pm5

    def predict_pm1(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> AutoACMGCriteria:
        """Override predict_pm1 to specify critical domains for DICER1."""
        logger.info("Predict PM1")

        gene_cluster = PM1_CLUSTER.get(var_data.hgnc_id, None)
        if not gene_cluster:
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.NotMet,
                strength=AutoACMGStrength.PathogenicModerate,
                summary="Variant does not affect a critical domain for DICER1.",
            )

        # Check moderate level criteria for metal ion-binding residues
        if var_data.prot_pos in gene_cluster.get("moderate", {}).get("residues", []):
            comment = (
                f"Variant affects a metal ion-binding residue in DICER1 "
                f"(position {var_data.prot_pos}). PM1 is met at the Moderate level."
            )
            return AutoACMGCriteria(
                name="PM1",
                prediction=AutoACMGPrediction.Met,
                strength=AutoACMGStrength.PathogenicModerate,
                summary=comment,
            )

        # Check supporting level criteria for RNase IIIb domain excluding metal ion-binding residues
        for start, end in gene_cluster.get("supporting", {}).get("domains", []):
            if start <= var_data.prot_pos <= end:
                comment = (
                    f"Variant affects a residue in the RNase IIIb domain of DICER1 "
                    f"(position {var_data.prot_pos}) outside the key metal ion-binding residues. "
                    f"PM1 is met at the Supporting level."
                )
                return AutoACMGCriteria(
                    name="PM1",
                    prediction=AutoACMGPrediction.Met,
                    strength=AutoACMGStrength.PathogenicSupporting,
                    summary=comment,
                )

        # If no criteria match
        return AutoACMGCriteria(
            name="PM1",
            prediction=AutoACMGPrediction.NotMet,
            strength=AutoACMGStrength.PathogenicModerate,
            summary="Variant does not affect a critical domain for DICER1.",
        )

    def predict_pm2ba1bs1bs2(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria, AutoACMGCriteria]:
        """Change the thresholds for PM2, BA1 and BS1."""
        var_data.thresholds.pm2_pathogenic = 0.000005
        var_data.thresholds.ba1_benign = 0.003
        var_data.thresholds.bs1_benign = 0.0003
        return super().predict_pm2ba1bs1bs2(seqvar, var_data)

    def _bp3_not_applicable(self, seqvar: SeqVar, var_data: AutoACMGSeqVarData) -> bool:
        """BP3 is not applicable for DICER1."""
        return True

    def predict_pp2bp1(
        self, seqvar: SeqVar, var_data: AutoACMGSeqVarData
    ) -> Tuple[AutoACMGCriteria, AutoACMGCriteria]:
        """Override PP2, BP1 for DICER1 to return not applicable status."""
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
        """Override donor and acceptor positions for DICER1 VCEP."""
        var_data.thresholds.bp7_donor = 7
        var_data.thresholds.bp7_acceptor = 21
        return super().predict_bp7(seqvar, var_data)
