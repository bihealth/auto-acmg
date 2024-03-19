"""PVS1 criteria for Sequence Variants (SeqVar)."""

import logging
import re
from enum import Enum, auto
from typing import Dict, List

from src.api.annonars import AnnonarsClient
from src.core.config import settings
from src.core.exceptions import InvalidAPIResposeError
from src.seqvar import SeqVar
from src.types import PVS1Prediction, SeqVarConsequence

# Setup logging
logging_level = logging.DEBUG if settings.DEBUG else logging.INFO
logging.basicConfig(level=logging_level)
logger = logging.getLogger(__name__)


class SeqVarPVS1:
    """PVS1 criteria for transcript."""

    def __init__(
        self,
        seqvar: SeqVar,
        seqvar_transcript: Dict,
        gene_transcript: Dict,
        consequence: SeqVarConsequence = SeqVarConsequence.NonsenseFrameshift,
    ):
        self.seqvar = seqvar
        self.consequence = consequence
        self.seqvar_transcripts = seqvar_transcript
        self.gene_transcripts = gene_transcript
        self.HGVS: str = ""
        self.pHGVS: str = ""
        self.tHGVS: str = ""
        self.gene_hgnc_id: str = ""
        self.prediction: PVS1Prediction = PVS1Prediction.NotPVS1
        self._initialize()

        # Frameshift/Nonsense attributes
        self.undergo_nmd: bool = False
        self.is_biologically_relevant: bool = False
        self.is_critical_for_protein_function: bool = False
        self.LoF_exon_is_frequent_in_pop: bool = False
        self.LoF_removes_10_percent_of_protein: bool = False

    def _initialize(self):
        """Setup the PVS1 class."""
        self.HGVS = self.gene_transcripts["id"]
        self.pHGVS = self.HGVS + ":" + self.seqvar_transcripts["hgvs_p"]
        self.tHGVS = self.HGVS + ":" + self.seqvar_transcripts["hgvs_t"]
        self.gene_hgnc_id = self.seqvar_transcripts["gene_id"]
        self.transcript_tags = self.seqvar_transcripts["feature_tag"]
        self.cds_sizes = [
            exon["altEndI"] - exon["altStartI"]
            for exon in self.gene_transcripts["genomeAlignments"][0]["exons"]
        ]
        self.cds_length = sum(self.cds_sizes)

    async def verify_PVS1(self):
        """Make the PVS1 prediction."""
        if self.consequence == SeqVarConsequence.NonsenseFrameshift:
            if self.gene_hgnc_id == "HGNC:9588":  # Follow guidelines for PTEN
                if self._get_pHGVS_termination(self.pHGVS) < 374:
                    self.prediction = PVS1Prediction.PVS1
                    return

            if self._undergo_nmd():
                if self._in_biologically_relevant_transcript(self.transcript_tags):
                    self.prediction = PVS1Prediction.PVS1
                else:
                    self.prediction = PVS1Prediction.NotPVS1
            else:
                if self._critical4protein_function():
                    self.prediction = PVS1Prediction.PVS1_Strong
                else:
                    # Guidelines are unclear on "and/or" statement here.
                    if await self._lof_is_frequent_in_population(
                        self.seqvar
                    ) or not self._in_biologically_relevant_transcript(self.transcript_tags):
                        self.prediction = PVS1Prediction.NotPVS1
                    else:
                        if self._lof_removes_more_then_10_percent_of_protein(
                            self.pHGVS, self.cds_length
                        ):
                            self.prediction = PVS1Prediction.PVS1_Strong
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate
        elif self.consequence == SeqVarConsequence.SpliceSites:
            # TODO: Implement SpliceSites PVS1 prediction
            pass
        elif self.consequence == SeqVarConsequence.InitiationCodon:
            # TODO: Implement InitiationCodon PVS1 prediction
            pass
        else:
            self.prediction = PVS1Prediction.NotPVS1
            logger.info(f"Consequence {self.consequence} is not supported for PVS1 prediction.")

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

    def _undergo_nmd(self) -> bool:
        """
        Nonsense-mediated decay (NMD) classification. Return if the variant
        undergoes NMD prediction.
        **Rule:** If the variant is located in the last exon or in the last 50 nucleotides
        of the penultimate exon, it is NOT predicted to undergo NMD.
        """
        # Ensure that necessary data is available
        assert self.cds_sizes is not None

        new_stop_codon = self._get_pHGVS_termination(self.pHGVS)
        if self.gene_hgnc_id == "HGNC:4284":  # Hearing Loss Guidelines GJB2
            return True
        elif len(self.cds_sizes) <= 1:
            return False
        else:
            nmd_cutoff = sum(self.cds_sizes[:-1]) - min(50, self.cds_sizes[-2])
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

    def _critical4protein_function(self) -> bool:
        """Check if the truncated/altered region is critical for the protein function."""
        # TODO: Implement this method
        return False

    @staticmethod
    async def _lof_is_frequent_in_population(seqvar) -> bool:
        """
        Check if the LoF variants in the exon are frequent in the general population.
        **Rule:** If the LoF variants in the exon > 0.1% in the general population, it is
        predicted to be frequent in the general population.
        """
        try:
            annonars_client = AnnonarsClient()
            # TODO: Use correct start and stop positions
            response = await annonars_client.get_variant_from_range(
                seqvar, seqvar.pos - 20, seqvar.pos + 20
            )
            if response:
                lof_variants = 0
                all_variants = len(response["result"]["gnomad_genomes"]["vep"])
                for variant in response["result"]["gnomad_genomes"]["vep"]:
                    if variant["consequence"] in ["frameshift_variant", "stop_gained"]:
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
            logger.error(
                f"Failed to get variant from range for variant {seqvar.user_representation}."
            )
            logger.error(e)
            return False

    @staticmethod
    def _lof_removes_more_then_10_percent_of_protein(pHGVS: str, cds_length: int) -> bool:
        """Check if the LoF variant removes more than 10% of the protein."""
        pattern = re.compile(r"p\.\D+(\d+)(\D+fs)?(\*|X|Ter)(\d+)?")
        match = pattern.search(pHGVS)
        codon_offset = int(match.group(1)) if match else -1
        codon_length = cds_length / 3
        if codon_offset > 0 and (codon_length - codon_offset) / codon_length > 0.1:
            return True
        else:
            return False
