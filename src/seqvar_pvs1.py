"""PVS1 criteria for Sequence Variants (SeqVar)."""

import re
from enum import Enum, auto
from typing import Dict

from src.seqvar import SeqVar
from src.types import PVS1Prediction, SeqVarConsequence


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
        self.in_relevant_transcript: bool = False
        self.critical4protein: bool = False
        self.is_frequent_in_population: bool = False
        self.remove_much_of_protein: bool = False

        # self.id = self.id_generator()
        # self.vep_input = f"/tmp/vep.{self.id}.vcf"
        # self.vep_output = f"/tmp/vep.{self.id}.tab"

        # self.vep_variation = "na"
        # self.vep_symbol = "na"
        # self.vep_trans = "na"
        # self.vep_canonical = "na"
        # self.vep_pick = "na"
        # self.vep_consequence = "na"
        # self.hgvs_c = "na"
        # self.hgvs_p = "na"
        # self.hgvs_g = "na"
        # self.vep_exon = "na"
        # self.vep_intron = "na"

    def _initialize(self):
        """Setup the PVS1 class."""
        self.HGVS = self.gene_transcripts["id"]
        self.pHGVS = self.HGVS + ":" + self.seqvar_transcripts["hgvs_p"]
        self.tHGVS = self.HGVS + ":" + self.seqvar_transcripts["hgvs_t"]
        self.gene_hgnc_id = self.seqvar_transcripts["gene_id"]

    def verify_PVS1(self):
        """Make the PVS1 prediction."""
        if self.consequence == SeqVarConsequence.NonsenseFrameshift:
            self.undergo_nmd = self._undergo_nmd()
            if self.undergo_nmd:
                self.in_relevant_transcript = self._in_biologically_relevant_transcript()
                if self.in_relevant_transcript:
                    self.prediction = PVS1Prediction.PVS1
                else:
                    self.prediction = PVS1Prediction.NotPVS1
            else:
                self.critical4protein = self._critical4protein_function()
                if self.critical4protein:
                    self.prediction = PVS1Prediction.PVS1_Strong
                else:
                    self.is_frequent_in_population = self._lof_is_frequent_in_population()
                    if self.is_frequent_in_population:
                        self.prediction = PVS1Prediction.NotPVS1
                    else:
                        self.remove_much_of_protein = self._remove_more_then_10_percent_of_protein()
                        if self.remove_much_of_protein:
                            self.prediction = PVS1Prediction.PVS1_Strong
                        else:
                            self.prediction = PVS1Prediction.PVS1_Moderate

    # @staticmethod
    # def id_generator(size: int = 6, chars: str = string.ascii_uppercase + string.digits) -> str:
    #     """Generates a unique identifier for the VEP."""
    #     return "".join(random.choice(chars) for _ in range(size))

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

        :param pHGVS: Protein HGVS
        :type pHGVS: str
        :return: Termination position
        :rtype: int
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
        See more at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6185798/#:~:text=Generally%2C%20NMD%20is%20not%20predicted%20to%20occur%20if%20the%20premature%20termination%20codon%20occurs%20in%20the%203%E2%80%99%20most%20exon%20or%20within%20the%203%E2%80%99%2Dmost%2050%20nucleotides%20of%20the%20penultimate%20exon

        :return: Variant undergoes NMD
        :rtype: bool
        """
        assert self.gene_transcripts is not None
        new_stop_codon = self._get_pHGVS_termination(self.pHGVS)
        cds_sizes = [
            exon["altEndI"] - exon["altStartI"]
            for exon in self.gene_transcripts["genomeAlignments"][0]["exons"]
        ]
        if self.gene_hgnc_id == "HGNC:4284":  # Hearing Loss Guidelines GJB2
            return True
        elif len(cds_sizes) <= 1:
            return False
        else:
            nmd_cutoff = sum(cds_sizes[:-1]) - min(50, cds_sizes[-2])
            return new_stop_codon * 3 <= nmd_cutoff

    def _in_biologically_relevant_transcript(self) -> bool:
        """
        Check if the exon is in a biologically relevant transcript.

        :return: Variant's exon in biologically relevant transcript(s)
        :rtype: bool
        """
        return "ManeSelect" in self.seqvar_transcripts["feature_tag"]

    def _critical4protein_function(self) -> bool:
        """
        Check if the truncated/altered region is critical for the protein function.

        :return: Variant is critical for protein function
        :rtype: bool
        """
        # TODO: Implement this method
        return False

    def _lof_is_frequent_in_population(self) -> bool:
        """
        Check if the LoF variants in the exon are frequent in the general population.

        :return: LoF variants are frequent
        :rtype: bool
        """
        return False

    def _remove_more_then_10_percent_of_protein(self) -> bool:
        """
        Check if the protein is 10% of the normal length.

        :return: Protein is 10% of the normal length
        :rtype: bool
        """
        return False
