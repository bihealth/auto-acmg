"""Implementations of the PVS1 algorithm."""

import logging
import re
from enum import Enum, auto
from typing import Dict, List

from src.api.mehari import MehariClient
from src.core.config import settings
from src.genome_builds import GenomeRelease
from src.seqvar import SeqVar, SeqVarResolver

logging_level = logging.DEBUG if settings.DEBUG else logging.INFO
logging.basicConfig(level=logging_level)
logger = logging.getLogger(__name__)


class SeqVarConsequence(Enum):
    """Consequence of a sequence variant."""

    NonsenseFrameshift = auto()
    SpliceSites = auto()
    InitiationCodon = auto()


class PVS1Prediction(Enum):
    """PVS1 prediction."""

    PVS1 = auto()
    PVS1_Strong = auto()
    PVS1_Moderate = auto()
    PVS1_Supporting = auto()
    NotPVS1 = auto()


class AutoPVS1:
    """AutoPVS1 algorithm for PVS1 criteria prediction."""

    def __init__(self, variant_name: str, genome_release: GenomeRelease = GenomeRelease.GRCh38):
        self.variant_name = variant_name
        self.genome_release = genome_release

    async def _get_transcripts_seqvar(self):
        """Get all transcripts for the given sequence variant."""
        assert self.seqvar is not None
        try:
            mehari_client = MehariClient()
            response = await mehari_client.get_seqvar_transcripts(self.seqvar)
            if response["result"]:
                self.seqvar_transcripts = response["result"]
            else:
                self.seqvar_transcripts = None

            if self.seqvar_transcripts and len(self.seqvar_transcripts) > 0:
                self.gene_hgnc_id = self.seqvar_transcripts[0]["gene_id"]

                for transcript in self.seqvar_transcripts:
                    self.HGVSs.append(transcript["feature_id"])

                result = await mehari_client.get_gene_transcripts(
                    self.gene_hgnc_id, self.seqvar.genome_release
                )

                if result["transcripts"]:
                    self.gene_transcripts = result["transcripts"]
                else:
                    self.gene_transcripts = None

        except Exception as e:
            logger.error(
                f"Failed to get transcripts for variant {self.seqvar.user_representation}."
            )
            logger.error(e)

    async def resolve_variant(self) -> SeqVar | None:
        """
        Resolve the variant.

        :return: Sequence variant
        :rtype: SeqVar | None
        """
        # TODO: Add resolve for Structure variants
        try:
            # Try to resolve as Sequence variant
            seqvar_resolver = SeqVarResolver()
            seqvar: SeqVar = await seqvar_resolver.resolve_seqvar(
                self.variant_name, self.genome_release
            )
            logger.debug(f"Resolved variant: {seqvar}.")
            return seqvar
        except Exception as e:
            logger.error(e)
            return None

    async def predict(self):
        """Run the AutoPVS1 algorithm."""
        logger.info(f"Running AutoPVS1 for variant {self.variant_name}.")
        variant = await self.resolve_variant()
        if isinstance(variant, SeqVar):
            self.seqvar = variant
            self.HGVSs: List[str] = []
            self.seqvar_transcripts = None
            self.gene_transcripts = None
            self.gene_hgnc_id = None
            self.predictions: Dict[str, PVS1Prediction] = {}

            logger.debug(f"Retrieving transcripts.")
            await self._get_transcripts_seqvar()
            if self.HGVSs:
                for hgvs in self.HGVSs:
                    logger.info(f"Analyzing transcript {hgvs}.")
                    # Choose transcripts
                    seqvar_t = None
                    gene_t = None
                    for transcript_info in self.seqvar_transcripts:
                        if transcript_info["feature_id"] == hgvs:
                            seqvar_t = transcript_info
                            break
                    for transcript_info in self.gene_transcripts:
                        if transcript_info["id"] == hgvs:
                            gene_t = transcript_info
                            break

                    if seqvar_t and gene_t:
                        # TODO: Add support for SpliceSites and InitiationCodon consequences
                        pvs1 = SeqVarPVS1(self.seqvar, seqvar_t, gene_t)
                        pvs1.verify_PVS1()
                        self.predictions[hgvs] = pvs1.prediction
                        logger.debug(f"PVS1 prediction for {hgvs}: {pvs1.prediction.name}")

            else:
                logger.debug("No transcripts found for the variant.")
        elif isinstance(variant, str):
            # TODO: Add Structure variants
            pass
        else:
            logger.error(f"Failed to resolve variant {self.variant_name}.")
            return


class SeqVarPVS1:
    """PVS1 criteria for transcript."""

    def __init__(
        self,
        seqvar: SeqVar,
        seqvar_transcripts: Dict,
        gene_transcripts: Dict,
        consequence: SeqVarConsequence = SeqVarConsequence.NonsenseFrameshift,
    ):
        self.seqvar = seqvar
        self.consequence = consequence
        self.seqvar_transcripts = seqvar_transcripts
        self.gene_transcripts = gene_transcripts
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
