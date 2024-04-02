"""PVS1 criteria for Sequence Variants (SeqVar)."""

import re
from typing import List

import typer

from src.api.annonars import AnnonarsClient
from src.core.exceptions import InvalidAPIResposeError
from src.seqvar import SeqVar
from src.types.enums import PVS1Prediction, PVS1PredictionSeqVarPath, SeqVarConsequence
from src.types.mehari import CdsPos, Exon, TranscriptGene, TranscriptSeqvar


class SetupSeqVarData:
    """Mixin class for setting up the SeqVar data."""

    def __init__(
        self,
        seqvar: SeqVar,
        seqvar_transcript: TranscriptSeqvar,
        gene_transcript: TranscriptGene,
        consequence: SeqVarConsequence = SeqVarConsequence.NonsenseFrameshift,
    ):
        self.seqvar = seqvar
        self.consequence = consequence
        self._seqvar_transcript = seqvar_transcript
        self._gene_transcript = gene_transcript
        self.HGVS: str = ""
        self.pHGVS: str = ""
        self.tHGVS: str = ""
        self.gene_hgnc_id: str = ""
        self.transcript_tags: List[str] = []
        self.exons: List[Exon] = []
        self.cds_pos: int | None = None

        self._initialize()

    def _initialize(self):
        """Setup the PVS1 class."""
        self.HGVS = self._gene_transcript.id
        self.pHGVS = self.HGVS + ":" + (self._seqvar_transcript.hgvs_p or "")
        self.tHGVS = self.HGVS + ":" + (self._seqvar_transcript.hgvs_t or "")
        self.gene_hgnc_id = self._seqvar_transcript.gene_id
        self.transcript_tags = self._seqvar_transcript.feature_tag
        self.exons = self._gene_transcript.genomeAlignments[0].exons
        self.cds_pos = (
            self._seqvar_transcript.cds_pos.ord
            if isinstance(self._seqvar_transcript.cds_pos, CdsPos)
            else None
        )


class SeqVarPVS1(SetupSeqVarData):
    """PVS1 criteria for transcript."""

    def __init__(self):
        self.prediction: PVS1Prediction = PVS1Prediction.NotPVS1
        self.prediction_path: PVS1PredictionSeqVarPath = PVS1PredictionSeqVarPath.NotSet

    def verify_PVS1(self):
        """Make the PVS1 prediction."""
        if self.consequence == SeqVarConsequence.NonsenseFrameshift:
            if self.gene_hgnc_id == "HGNC:9588":  # Follow guidelines for PTEN
                if self._get_pHGVS_termination(self.pHGVS) < 374:
                    self.prediction = PVS1Prediction.PVS1
                    return

            if self._undergo_nmd(self.exons, self.pHGVS, self.gene_hgnc_id):
                if self._in_biologically_relevant_transcript(self.transcript_tags):
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.NF1
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.NF2
            else:
                if self._critical4protein_function(self.seqvar, self.cds_pos, self.exons):
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionSeqVarPath.NF3
                else:
                    if self._lof_is_frequent_in_population(
                        self.seqvar,
                        self._gene_transcript.genomeAlignments[0].cdsStart,
                        self._gene_transcript.genomeAlignments[0].cdsEnd,
                    ) or not self._in_biologically_relevant_transcript(self.transcript_tags):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionSeqVarPath.NF4
                    else:
                        if self._lof_removes_more_then_10_percent_of_protein(
                            self.pHGVS, self.exons
                        ):
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionSeqVarPath.NF5
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionSeqVarPath.NF6

        elif self.consequence == SeqVarConsequence.SpliceSites:
            if self._exon_skipping_or_cryptic_ss_disruption() and self._undergo_nmd(
                self.exons, self.pHGVS, self.gene_hgnc_id
            ):
                if self._in_biologically_relevant_transcript(self.transcript_tags):
                    self.prediction = PVS1Prediction.PVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.SS1
                else:
                    self.prediction = PVS1Prediction.NotPVS1
                    self.prediction_path = PVS1PredictionSeqVarPath.SS2
            elif self._exon_skipping_or_cryptic_ss_disruption() and not self._undergo_nmd(
                self.exons, self.pHGVS, self.gene_hgnc_id
            ):
                if self._critical4protein_function(self.seqvar, self.cds_pos, self.exons):
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionSeqVarPath.SS3
                else:
                    if self._lof_is_frequent_in_population(
                        self.seqvar,
                        self._gene_transcript.genomeAlignments[0].cdsStart,
                        self._gene_transcript.genomeAlignments[0].cdsEnd,
                    ) or not self._in_biologically_relevant_transcript(self.transcript_tags):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionSeqVarPath.SS4
                    else:
                        if self._lof_removes_more_then_10_percent_of_protein(
                            self.pHGVS, self.exons
                        ):
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionSeqVarPath.SS5
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionSeqVarPath.SS6
            else:
                if self._critical4protein_function(self.seqvar, self.cds_pos, self.exons):
                    self.prediction = PVS1Prediction.PVS1_Strong
                    self.prediction_path = PVS1PredictionSeqVarPath.SS10
                else:
                    if self._lof_is_frequent_in_population(
                        self.seqvar,
                        self._gene_transcript.genomeAlignments[0].cdsStart,
                        self._gene_transcript.genomeAlignments[0].cdsEnd,
                    ) or not self._in_biologically_relevant_transcript(self.transcript_tags):
                        self.prediction = PVS1Prediction.NotPVS1
                        self.prediction_path = PVS1PredictionSeqVarPath.SS7
                    else:
                        if self._lof_removes_more_then_10_percent_of_protein(
                            self.pHGVS, self.exons
                        ):
                            self.prediction = PVS1Prediction.PVS1_Strong
                            self.prediction_path = PVS1PredictionSeqVarPath.SS8
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
                            self.prediction_path = PVS1PredictionSeqVarPath.SS9

        elif self.consequence == SeqVarConsequence.InitiationCodon:
            if self._alternative_start_codon():
                self.prediction = PVS1Prediction.NotPVS1
                self.prediction_path = PVS1PredictionSeqVarPath.IC3
            else:
                if self._upstream_pathogenic_variant():
                    self.prediction = PVS1Prediction.PVS1_Moderate
                    self.prediction_path = PVS1PredictionSeqVarPath.IC1
                else:
                    self.prediction = PVS1Prediction.PVS1_Supporting
                    self.prediction_path = PVS1PredictionSeqVarPath.IC2
        else:
            self.prediction = PVS1Prediction.NotPVS1
            typer.echo(f"Consequence {self.consequence} is not supported for PVS1 prediction.")

    @staticmethod
    def _get_pHGVS_termination(pHGVS: str) -> int:
        """
        Get termination position from pHGVS.
        **Note:** If the position is not found, return -1.
        Examples:
        - NM_031475.2:p.Gln98*
        - NM_031475.2:p.Ala586Glyfs*73
        - NP_000305.3:p.Arg378SerfsTer5
        - p.Arg97Glyfs*26 (alternatively p.Arg97GlyfsTer26, or short p.Arg97fs)
        """
        if "fs" in pHGVS:  # If frameshift
            pattern1 = re.compile(r"p\.\D+(\d+)\D+fs(\*|X|Ter)(\d+)")
            match1 = pattern1.search(pHGVS)
            pattern2 = re.compile(r"p\.\D+(\d+)fs")
            match2 = pattern2.search(pHGVS)

            if match1:
                # if int(match1.group(1)) / (self.transcript.cds_length/3) > 0.5:
                termination = int(match1.group(1)) + int(match1.group(3))
                # else:
                #    termination = int((self.transcript.cds_length/3)/2)
            elif match2:
                termination = int(match2.group(1))
            else:
                termination = -1

        elif [char in pHGVS for char in ["*", "X", "Ter"]]:  # If premature termination codon
            pattern = re.compile(r"p\.\D+(\d+)(\*|X|Ter)")
            match = pattern.search(pHGVS)
            termination = int(match.group(1)) if match else -1
        else:
            termination = -1
        return termination

    def _undergo_nmd(self, exons: List[Exon], pHGVS: str, hgnc_id: str) -> bool:
        """
        Nonsense-mediated decay (NMD) classification. Return if the variant
        undergoes NMD prediction.
        **Rule:** If the variant is located in the last exon or in the last 50 nucleotides
        of the penultimate exon, it is NOT predicted to undergo NMD.
        """
        new_stop_codon = self._get_pHGVS_termination(pHGVS)
        cds_sizes = [exon.altEndI - exon.altStartI for exon in exons]
        if hgnc_id == "HGNC:4284":  # Hearing Loss Guidelines GJB2
            return True
        elif len(cds_sizes) <= 1:
            return False
        else:
            nmd_cutoff = sum(cds_sizes[:-1]) - min(50, cds_sizes[-2])
            return new_stop_codon * 3 <= nmd_cutoff

    @staticmethod
    def _in_biologically_relevant_transcript(transcript_tags: List[str]) -> bool:
        """
        Check if the exon with SeqVar is in a biologically relevant transcript.
        **Rule:** If the variant is located in a transcript with a MANE Select tag, it is
        predicted to be in a biologically relevant transcript.
        """
        # Ensure that necessary data is available
        return "ManeSelect" in transcript_tags

    @staticmethod
    def _critical4protein_function(seqvar: SeqVar, cds_pos: int | None, exons: List[Exon]) -> bool:
        """
        Check if the truncated/altered region is critical for the protein function.
        **Rule:** Affect of the variant on the protein function is indicated
        by experimental or clinical evidence:
        - Presence of Pathogenic variants downstream of the new stop codon
        """
        if not cds_pos:
            return False
        start_pos = 0
        for exon in exons:
            if exon.altCdsStartI < cds_pos and exon.altCdsEndI > cds_pos:
                start_pos = exon.altStartI + (cds_pos - exon.altCdsStartI)
                break
        end_pos = exons[-1].altEndI

        try:
            annonars_client = AnnonarsClient()
            response = annonars_client.get_variant_from_range(seqvar, start_pos, end_pos)
            if (
                response
                and response.result.gnomad_genomes
                and response.result.gnomad_genomes[0].vep
            ):
                pathogenic_variants = 0
                for variant in response.result.gnomad_genomes[0].vep:
                    # TODO: count only pathogenic variants
                    if variant.consequence in ["pathogenic"]:
                        pathogenic_variants += 1
                # TODO: Proove that this is the correct threshold
                if (
                    pathogenic_variants > 2
                    and pathogenic_variants / len(response.result.gnomad_genomes[0].vep) > 0.005
                ):
                    return True
                return False
            else:
                raise InvalidAPIResposeError("Failed to get variant from range.")
        except Exception as e:
            typer.secho(
                f"Failed to get variant from range for variant {seqvar.user_representation}.",
                err=True,
                fg=typer.colors.RED,
            )
            typer.echo(e, err=True)
            return False

    @staticmethod
    def _lof_is_frequent_in_population(seqvar: SeqVar, start: int, end: int) -> bool:
        """
        Check if the LoF variants in the exon are frequent in the general population.
        **Rule:** If the LoF variants in the exon > 0.1% in the general population, it is
        predicted to be frequent in the general population.
        """
        try:
            annonars_client = AnnonarsClient()
            response = annonars_client.get_variant_from_range(seqvar, start, end)
            if (
                response
                and response.result.gnomad_genomes
                and response.result.gnomad_genomes[0].vep
            ):
                lof_variants = 0
                all_variants = len(response.result.gnomad_genomes[0].vep)
                for variant in response.result.gnomad_genomes[0].vep:
                    if variant.consequence in ["frameshift_variant", "stop_gained"]:
                        lof_variants += 1
                # TODO: Proove that this is the correct threshold
                # The guideline does not specify the exact threshold. 0.1% was taken from AutoPVS1
                if lof_variants / all_variants > 0.001:
                    return True
                else:
                    return False
            else:
                raise InvalidAPIResposeError("Failed to get variant from range.")
        except Exception as e:
            typer.secho(
                f"Failed to get variant from range for variant {seqvar.user_representation}.",
                err=True,
                fg=typer.colors.RED,
            )
            typer.echo(e, err=True)
            return False

    @staticmethod
    def _lof_removes_more_then_10_percent_of_protein(pHGVS: str, exons: List[Exon]) -> bool:
        """Check if the LoF variant removes more than 10% of the protein."""
        cds_length = sum([exon.altEndI - exon.altStartI for exon in exons])
        pattern = re.compile(r"p\.\D+(\d+)(\D+fs)?(\*|X|Ter)(\d+)?")
        match = pattern.search(pHGVS)
        codon_offset = int(match.group(1)) if match else -1
        codon_length = cds_length / 3
        if codon_offset > 0 and (codon_length - codon_offset) / codon_length > 0.1:
            return True
        else:
            return False

    def _exon_skipping_or_cryptic_ss_disruption(self) -> bool:
        """Check if the variant causes exon skipping or cryptic splice site disruption."""
        # TODO: Implement this method
        return False

    def _alternative_start_codon(self) -> bool:
        """Check if the variant introduces an alternative start codon."""
        # TODO: Implement this method
        return False

    def _upstream_pathogenic_variant(self) -> bool:
        """Check if the variant is an upstream pathogenic variant."""
        # TODO: Implement this method
        return False
