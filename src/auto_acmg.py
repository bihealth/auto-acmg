"""Implementations of the PVS1 algorithm."""

from typing import Dict, Optional, Type, Union

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.annonars_variant import VariantResult
from src.defs.auto_acmg import AutoACMGSeqVarResult, AutoACMGStrucVarResult, CdsInfo, GenomicStrand
from src.defs.exceptions import AlgorithmError, AutoAcmgBaseException, ParseError
from src.defs.genome_builds import GenomeRelease
from src.defs.mehari import CdsPos, ProteinPos, TxPos
from src.defs.seqvar import SeqVar, SeqVarResolver
from src.defs.strucvar import StrucVar, StrucVarResolver
from src.seqvar.auto_bp7 import AutoBP7
from src.seqvar.auto_pm1 import AutoPM1
from src.seqvar.auto_pm2_ba1_bs1_bs2 import AutoPM2BA1BS1BS2
from src.seqvar.auto_pm4_bp3 import AutoPM4BP3
from src.seqvar.auto_pp2_bp1 import AutoPP2BP1
from src.seqvar.auto_pp3_bp4 import AutoPP3BP4
from src.seqvar.auto_ps1_pm5 import AutoPS1PM5
from src.seqvar.auto_pvs1 import AutoPVS1
from src.seqvar.default_predictor import DefaultSeqVarPredictor
from src.strucvar.default_predictor import DefaultStrucVarPredictor
from src.utils import SeqVarTranscriptsHelper, StrucVarTranscriptsHelper
from src.vcep import (
    ACADVLPredictor,
    BrainMalformationsPredictor,
    CardiomyopathyPredictor,
    CDH1Predictor,
    CerebralCreatineDeficiencySyndromesPredictor,
    CoagulationFactorDeficiencyPredictor,
    CongenitalMyopathiesPredictor,
    DICER1Predictor,
    ENIGMAPredictor,
    EpilepsySodiumChannelPredictor,
    FamilialHypercholesterolemiaPredictor,
    FBN1Predictor,
    GlaucomaPredictor,
    HBOPCPredictor,
    HearingLossPredictor,
    HHTPredictor,
    InsightColorectalCancerPredictor,
    LeberCongenitalAmaurosisPredictor,
    LysosomalDiseasesPredictor,
    MalignantHyperthermiaPredictor,
    MitochondrialDiseasesPredictor,
    MonogenicDiabetesPredictor,
    MyeloidMalignancyPredictor,
    PKUPredictor,
    PlateletDisordersPredictor,
    PTENPredictor,
    PulmonaryHypertensionPredictor,
    RASopathyPredictor,
    RettAngelmanPredictor,
    SCIDPredictor,
    ThrombosisPredictor,
    TP53Predictor,
    VHLPredictor,
    VonWillebrandDiseasePredictor,
)

#: Mapping of HGNC gene identifiers to predictor classes.
VCEP_MAPPING = {
    "HGNC:92": ACADVLPredictor,  # ACADVL
    "HGNC:393": BrainMalformationsPredictor,  # AKT3
    "HGNC:3942": BrainMalformationsPredictor,  # MTOR
    "HGNC:8975": BrainMalformationsPredictor,  # PIK3CA
    "HGNC:8980": BrainMalformationsPredictor,  # PIK3R2
    "HGNC:7577": CardiomyopathyPredictor,  # MYH7
    "HGNC:7551": CardiomyopathyPredictor,  # MYBPC3
    "HGNC:11947": CardiomyopathyPredictor,  # TNNI3
    "HGNC:11949": CardiomyopathyPredictor,  # TNNT2
    "HGNC:12010": CardiomyopathyPredictor,  # TPM1
    "HGNC:143": CardiomyopathyPredictor,  # ACTC1
    "HGNC:7583": CardiomyopathyPredictor,  # MYL2
    "HGNC:7584": CardiomyopathyPredictor,  # MYL3
    "HGNC:1748": CDH1Predictor,  # CDH1
    "HGNC:4175": CerebralCreatineDeficiencySyndromesPredictor,  # GATM
    "HGNC:4136": CerebralCreatineDeficiencySyndromesPredictor,  # GAMT
    "HGNC:11055": CerebralCreatineDeficiencySyndromesPredictor,  # SLC6A8
    "HGNC:3546": CoagulationFactorDeficiencyPredictor,  # F8
    "HGNC:3551": CoagulationFactorDeficiencyPredictor,  # F9
    "HGNC:7720": CongenitalMyopathiesPredictor,  # NEB
    "HGNC:129": CongenitalMyopathiesPredictor,  # ACTA1
    "HGNC:2974": CongenitalMyopathiesPredictor,  # DNM2
    "HGNC:7448": CongenitalMyopathiesPredictor,  # MTM1
    "HGNC:10483": CongenitalMyopathiesPredictor,  # RYR1
    "HGNC:17098": DICER1Predictor,  # DICER1
    "HGNC:1100": ENIGMAPredictor,  # BRCA1
    "HGNC:1101": ENIGMAPredictor,  # BRCA2
    "HGNC:10585": EpilepsySodiumChannelPredictor,  # SCN1A
    "HGNC:10588": EpilepsySodiumChannelPredictor,  # SCN2A
    "HGNC:10590": EpilepsySodiumChannelPredictor,  # SCN3A
    "HGNC:10596": EpilepsySodiumChannelPredictor,  # SCN8A
    "HGNC:10586": EpilepsySodiumChannelPredictor,  # SCN1B
    "HGNC:6547": FamilialHypercholesterolemiaPredictor,  # LDLR
    "HGNC:3603": FBN1Predictor,  # FBN1
    "HGNC:7610": GlaucomaPredictor,  # MYOC
    "HGNC:13733": HearingLossPredictor,  # CDH23
    "HGNC:2180": HearingLossPredictor,  # COCH
    "HGNC:4284": HearingLossPredictor,  # GJB2
    "HGNC:6298": HearingLossPredictor,  # KCNQ4
    "HGNC:7605": HearingLossPredictor,  # MYO6
    "HGNC:7606": HearingLossPredictor,  # MYO7A
    "HGNC:8818": HearingLossPredictor,  # SLC26A4
    "HGNC:11720": HearingLossPredictor,  # TECTA
    "HGNC:12601": HearingLossPredictor,  # USH2A
    "HGNC:7594": HearingLossPredictor,  # MYO15A
    "HGNC:8515": HearingLossPredictor,  # OTOF
    "HGNC:795": HBOPCPredictor,  # ATM
    "HGNC:26144": HBOPCPredictor,  # PALB2
    "HGNC:175": HHTPredictor,  # ACVRL1
    "HGNC:3349": HHTPredictor,  # ENG
    "HGNC:583": InsightColorectalCancerPredictor,  # APC
    "HGNC:7127": InsightColorectalCancerPredictor,  # MLH1
    "HGNC:7325": InsightColorectalCancerPredictor,  # MSH2
    "HGNC:7329": InsightColorectalCancerPredictor,  # MSH6
    "HGNC:9122": InsightColorectalCancerPredictor,  # PMS2
    "HGNC:10294": LeberCongenitalAmaurosisPredictor,  # RPE65
    "HGNC:4065": LysosomalDiseasesPredictor,  # GAA
    "HGNC:10483": MalignantHyperthermiaPredictor,  # GBA
    "HGNC:23287": MitochondrialDiseasesPredictor,  # ETHE1
    "HGNC:8806": MitochondrialDiseasesPredictor,  # PDHA1
    "HGNC:9179": MitochondrialDiseasesPredictor,  # POLG
    "HGNC:16266": MitochondrialDiseasesPredictor,  # SLC19A3
    "HGNC:11621": MonogenicDiabetesPredictor,  # HNF1A
    "HGNC:5024": MonogenicDiabetesPredictor,  # HNF4A
    "HGNC:4195": MonogenicDiabetesPredictor,  # GCK
    "HGNC:10471": MyeloidMalignancyPredictor,  # RUNX1
    "HGNC:8582": PKUPredictor,  # PAH
    "HGNC:6138": PlateletDisordersPredictor,  # ITGA2B
    "HGNC:6156": PlateletDisordersPredictor,  # ITGB3
    "HGNC:9588": PTENPredictor,  # PTEN
    "HGNC:1078": PulmonaryHypertensionPredictor,  # BMPR2
    "HGNC:12726": VonWillebrandDiseasePredictor,  # VWF
    "HGNC:775": ThrombosisPredictor,  # SERPINC1
    "HGNC:11998": TP53Predictor,  # TP53
    "HGNC:12687": VHLPredictor,  # VHL
    "HGNC:15454": RASopathyPredictor,  # SHOC2
    "HGNC:7989": RASopathyPredictor,  # NRAS
    "HGNC:9829": RASopathyPredictor,  # RAF1
    "HGNC:11187": RASopathyPredictor,  # SOS1
    "HGNC:11188": RASopathyPredictor,  # SOS2
    "HGNC:9644": RASopathyPredictor,  # PTPN11
    "HGNC:6407": RASopathyPredictor,  # KRAS
    "HGNC:6840": RASopathyPredictor,  # MAP2K1
    "HGNC:5173": RASopathyPredictor,  # HRAS
    "HGNC:10023": RASopathyPredictor,  # RIT1
    "HGNC:6842": RASopathyPredictor,  # MAP2K2
    "HGNC:1097": RASopathyPredictor,  # BRAF
    "HGNC:7227": RASopathyPredictor,  # MRAS
    "HGNC:6742": RASopathyPredictor,  # LZTR1
    "HGNC:17271": RASopathyPredictor,  # RRAS2
    "HGNC:9282": RASopathyPredictor,  # PPP1CB
    "HGNC:11634": RettAngelmanPredictor,  # TCF4
    "HGNC:11079": RettAngelmanPredictor,  # SLC9A6
    "HGNC:11411": RettAngelmanPredictor,  # CDKL5
    "HGNC:3811": RettAngelmanPredictor,  # FOXG1
    "HGNC:6990": RettAngelmanPredictor,  # MECP2
    "HGNC:12496": RettAngelmanPredictor,  # UBE3A
    "HGNC:12765": SCIDPredictor,  # FOXN1
    "HGNC:186": SCIDPredictor,  # ADA
    "HGNC:17642": SCIDPredictor,  # DCLRE1C
    "HGNC:6024": SCIDPredictor,  # IL7R
    "HGNC:6193": SCIDPredictor,  # JAK3
    "HGNC:9831": SCIDPredictor,  # RAG1
    "HGNC:9832": SCIDPredictor,  # RAG2
    "HGNC:6010": SCIDPredictor,  # IL2RG
}


class AutoACMG:
    """Class for predicting ACMG criteria.

    This class handles both sequence variants and structural variants to determine their potential
    impact under the various criteria of the ACMG guidelines for variant classification. Currently
    it only implements the PVS1 criterion (not finished yet).

    Attributes:
        variant_name (str): The name or identifier of the variant being analyzed.
        genome_release (GenomeRelease): The genome release version, defaults to GRCh38.
    """

    def __init__(
        self,
        variant_name: str,
        genome_release: GenomeRelease = GenomeRelease.GRCh38,
        *,
        config: Optional[Config] = None,
    ):
        """Initializes the AutoACMG with the specified variant and genome release.

        Args:
            variant_name: The name or identifier of the variant.
            genome_release (Optional): The genome release version, such as GRCh38 or GRCh37.
        """
        #: Configuration to use.
        self.config = config or Config()
        #: Annonars client.
        self.annonars_client: AnnonarsClient = AnnonarsClient(
            api_base_url=self.config.api_base_url_annonars
        )
        #: The name or identifier of the variant.
        self.variant_name = variant_name
        #: The genome release version.
        self.genome_release = genome_release
        #: The resolved sequence variant.
        self.seqvar: Optional[SeqVar] = None
        #: The resolved structural variant.
        self.strucvar: Optional[StrucVar] = None
        #: Prediction result for sequence variants.
        self.seqvar_result: AutoACMGSeqVarResult = AutoACMGSeqVarResult()
        #: Prediction result for structural variants.
        self.strucvar_result: AutoACMGStrucVarResult = AutoACMGStrucVarResult()

    def _get_variant_info(self, seqvar: SeqVar) -> Optional[VariantResult]:
        """
        Get variant information from Annonars.

        Args:
            seqvar: The variant to get information for.

        Returns:
            Optional[VariantResult]: The variant information. None if the information is not
            available.
        """
        try:
            logger.debug("Getting variant information for {}.", seqvar)
            return self.annonars_client.get_variant_info(seqvar).result
        except AutoAcmgBaseException as e:
            logger.error("Failed to get variant information. Error: {}", e)
            return None

    def _convert_score_val(self, score_value: Optional[Union[str, float, int]]) -> Optional[float]:
        """
        Convert score value to float.

        Since the score values can be represented as strings (with ";" as separator), we pick the
        maximum value that is not empty ("."). If the value is already numeric, we return it as is.

        Args:
            score_value (Optional[Union[str, float, int]]): Score value to convert.

        Returns:
            Optional[float]: Converted score value.

        Raises:
            AlgorithmError: If the score value cannot be converted to float.
        """
        if score_value is None:
            return None
        if isinstance(score_value, (float, int)):
            return float(score_value)
        try:
            scores = [float(score) for score in score_value.split(";") if score not in [".", ""]]
            return max(scores) if scores else None
        except ValueError as e:
            raise AlgorithmError("Failed to convert score value to float.") from e

    def resolve_variant(self) -> Union[SeqVar, StrucVar, None]:
        """Attempts to resolve the specified variant as either a sequence or structural variant.

        This method first tries to resolve the variant as a sequence variant. If it fails, it then
        attempts to resolve it as a structural variant.

        Returns:
            SeqVar or None: The resolved variant object or None if resolution fails.

        Raises:
            Exception: Specific exceptions are caught and logged, but generic exceptions may be
            raised if both resolutions fail.
        """
        logger.debug("Resolving variant: {}", self.variant_name)
        try:
            seqvar_resolver = SeqVarResolver(config=self.config)
            seqvar: SeqVar = seqvar_resolver.resolve_seqvar(self.variant_name, self.genome_release)
            logger.debug("Resolved sequence variant: {}", seqvar)
            return seqvar
        except ParseError:
            logger.debug("Failed to resolve sequence variant.")
            try:
                logger.debug("Trying to resolve structural variant.")
                strucvar_resolver = StrucVarResolver(config=self.config)
                strucvar: StrucVar = strucvar_resolver.resolve_strucvar(
                    self.variant_name, self.genome_release
                )
                return strucvar
            except ParseError as e:
                logger.debug("Failed to resolve the structure variant")
                return None
        except AutoAcmgBaseException as e:
            logger.error("An unexpected error occurred: {}", e)
            return None

    def parse_seqvar_data(self, seqvar: SeqVar) -> AutoACMGSeqVarResult:
        """Parses the data for the prediction."""
        # Mehari data
        ts_helper = SeqVarTranscriptsHelper(seqvar)
        ts_helper.initialize()
        seqvar_transcript, gene_transcript, _, all_genes_tx, _ = ts_helper.get_ts_info()
        if not seqvar_transcript or not gene_transcript:
            raise AutoAcmgBaseException("Transcript information is missing.")

        self.seqvar_result.data.consequence.mehari = seqvar_transcript.consequences
        self.seqvar_result.data.gene_symbol = seqvar_transcript.gene_symbol
        self.seqvar_result.data.hgnc_id = seqvar_transcript.gene_id
        self.seqvar_result.data.transcript_id = seqvar_transcript.feature_id
        self.seqvar_result.data.transcript_tags = seqvar_transcript.feature_tag
        self.seqvar_result.data.tx_pos_utr = (
            seqvar_transcript.tx_pos.ord if isinstance(seqvar_transcript.tx_pos, TxPos) else -1
        )
        self.seqvar_result.data.cds_pos = (
            seqvar_transcript.cds_pos.ord if isinstance(seqvar_transcript.cds_pos, CdsPos) else 0
        )
        self.seqvar_result.data.prot_pos = (
            seqvar_transcript.protein_pos.ord
            if isinstance(seqvar_transcript.protein_pos, ProteinPos)
            else -1
        )
        self.seqvar_result.data.prot_length = (
            seqvar_transcript.protein_pos.total
            if isinstance(seqvar_transcript.protein_pos, ProteinPos)
            else -1
        )
        self.seqvar_result.data.cds_info = {
            ts.id: CdsInfo(
                start_codon=ts.startCodon,
                stop_codon=ts.stopCodon,
                cds_start=ts.genomeAlignments[0].cdsStart,
                cds_end=ts.genomeAlignments[0].cdsEnd,
                cds_strand=GenomicStrand.from_string(ts.genomeAlignments[0].strand),
                exons=ts.genomeAlignments[0].exons,
            )
            for ts in all_genes_tx
        }
        self.seqvar_result.data.cds_start = gene_transcript.genomeAlignments[0].cdsStart
        self.seqvar_result.data.cds_end = gene_transcript.genomeAlignments[0].cdsEnd
        self.seqvar_result.data.strand = GenomicStrand.from_string(
            gene_transcript.genomeAlignments[0].strand
        )
        self.seqvar_result.data.exons = gene_transcript.genomeAlignments[0].exons

        # Annonars data
        variant_info = self._get_variant_info(seqvar)
        if not variant_info:
            raise AutoAcmgBaseException("Failed to get variant information.")

        if cadd := variant_info.cadd:
            self.seqvar_result.data.consequence.cadd = cadd.ConsDetail or ""
            self.seqvar_result.data.consequence.cadd_consequence = cadd.Consequence or ""
            self.seqvar_result.data.scores.cadd.phyloP100 = cadd.verPhyloP
            self.seqvar_result.data.scores.cadd.gerp = cadd.GerpS
            self.seqvar_result.data.scores.cadd.spliceAI_acceptor_gain = cadd.SpliceAI_acc_gain
            self.seqvar_result.data.scores.cadd.spliceAI_acceptor_loss = cadd.SpliceAI_acc_loss
            self.seqvar_result.data.scores.cadd.spliceAI_donor_gain = cadd.SpliceAI_don_gain
            self.seqvar_result.data.scores.cadd.spliceAI_donor_loss = cadd.SpliceAI_don_loss
            self.seqvar_result.data.scores.cadd.ada = cadd.dbscSNV_ada
            self.seqvar_result.data.scores.cadd.rf = cadd.dbscSNV_rf
        if dbsnfp := variant_info.dbnsfp:
            self.seqvar_result.data.pHGVS = dbsnfp.HGVSp_VEP or ""
            self.seqvar_result.data.scores.dbnsfp.alpha_missense = self._convert_score_val(
                dbsnfp.AlphaMissense_score
            )
            self.seqvar_result.data.scores.dbnsfp.metaRNN = self._convert_score_val(
                dbsnfp.MetaRNN_score
            )
            self.seqvar_result.data.scores.dbnsfp.bayesDel_noAF = self._convert_score_val(
                dbsnfp.BayesDel_noAF_score
            )
            self.seqvar_result.data.scores.dbnsfp.revel = self._convert_score_val(
                dbsnfp.REVEL_score
            )
            self.seqvar_result.data.scores.dbnsfp.phyloP100 = self._convert_score_val(
                dbsnfp.phyloP100way_vertebrate
            )
            self.seqvar_result.data.scores.dbnsfp.sift = self._convert_score_val(dbsnfp.SIFT_score)
            self.seqvar_result.data.scores.dbnsfp.polyphen2 = self._convert_score_val(
                dbsnfp.Polyphen2_HVAR_score
            )
            self.seqvar_result.data.scores.dbnsfp.mutationTaster = self._convert_score_val(
                dbsnfp.MutationTaster_score
            )
            self.seqvar_result.data.scores.dbnsfp.fathmm = self._convert_score_val(
                dbsnfp.FATHMM_score
            )
            self.seqvar_result.data.scores.dbnsfp.provean = self._convert_score_val(
                dbsnfp.PROVEAN_score
            )
            self.seqvar_result.data.scores.dbnsfp.vest4 = self._convert_score_val(
                dbsnfp.VEST4_score
            )
            self.seqvar_result.data.scores.dbnsfp.mutpred = self._convert_score_val(
                dbsnfp.MutPred_score
            )
            self.seqvar_result.data.scores.dbnsfp.primateAI = self._convert_score_val(
                dbsnfp.PrimateAI_score
            )
        if dbscsnv := variant_info.dbscsnv:
            self.seqvar_result.data.scores.dbscsnv.ada = dbscsnv.ada_score
            self.seqvar_result.data.scores.dbscsnv.rf = dbscsnv.rf_score
        self.seqvar_result.data.gnomad_exomes = variant_info.gnomad_exomes
        self.seqvar_result.data.gnomad_mtdna = variant_info.gnomad_mtdna

        # Gene info from Annonars
        gene_info = self.annonars_client.get_gene_info(self.seqvar_result.data.hgnc_id)
        if not gene_info:
            raise AutoAcmgBaseException("Failed to get gene information.")
        if gnomad_constraints := gene_info.genes.root[
            self.seqvar_result.data.hgnc_id
        ].gnomadConstraints:
            self.seqvar_result.data.scores.misZ = gnomad_constraints.misZ
        return self.seqvar_result

    def parse_strucvar_data(self, strucvar: StrucVar) -> AutoACMGStrucVarResult:
        """Parse the data for the prediction."""
        # Mehari data
        ts_helper = StrucVarTranscriptsHelper(strucvar)
        ts_helper.initialize()
        strucvar_transcript, gene_transcript, all_strucvar_tx, all_genes_tx = (
            ts_helper.get_ts_info()
        )
        if not strucvar_transcript or not gene_transcript:
            raise AutoAcmgBaseException("Transcript information is missing.")

        self.strucvar_result.data.hgnc_id = strucvar_transcript.hgnc_id
        self.strucvar_result.data.gene_symbol = gene_transcript.geneSymbol
        self.strucvar_result.data.transcript_id = gene_transcript.id
        self.strucvar_result.data.transcript_tags = gene_transcript.tags or []
        self.strucvar_result.data.strand = GenomicStrand.from_string(
            gene_transcript.genomeAlignments[0].strand
        )
        self.strucvar_result.data.exons = gene_transcript.genomeAlignments[0].exons
        self.strucvar_result.data.start_cdn = gene_transcript.startCodon
        self.strucvar_result.data.stop_cdn = gene_transcript.stopCodon

        return self.strucvar_result

    def select_predictor(self, hgnc_id: str) -> Type[DefaultSeqVarPredictor]:
        """Selects the predictor for the specified gene.

        Args:
            hgnc_id: The HGNC gene identifier.

        Returns:
            Type[Predictor]: The predictor class for the specified gene.
        """
        if hgnc_id in VCEP_MAPPING:
            return VCEP_MAPPING[hgnc_id]
        return DefaultSeqVarPredictor

    def predict(self) -> Union[AutoACMGSeqVarResult, AutoACMGStrucVarResult, None]:
        """
        Predict ACMG criteria for the specified variant.

        This method first resolves the variant, then predicts the PVS1 criterion for sequence
        variants and other ACMG criteria.

        Note:
            The method can resolve both sequence and structural variants, but currently only
            sequence variants are supported for ACMG criteria prediction.

        Raises:
            Exception: Specific exceptions are caught and logged, but generic exceptions may be
            raised if the prediction fails.
        """
        logger.info("Predicting ACMG criteria for variant: {}", self.variant_name)
        variant = self.resolve_variant()
        if not variant:
            logger.error("Failed to resolve the variant.")
            return None
        elif isinstance(variant, SeqVar):
            self.seqvar = variant
            self.seqvar_result.seqvar = self.seqvar
        elif isinstance(variant, StrucVar):
            self.strucvar = variant
            self.strucvar_result.strucvar = self.strucvar

        if isinstance(self.seqvar, SeqVar):
            if not self.seqvar:
                logger.error("Failed to resolve the sequence variant.")
                return None
            # ====== Setup the data ======
            self.parse_seqvar_data(self.seqvar)

            # ====== Predict ======
            predictor_class = self.select_predictor(self.seqvar_result.data.hgnc_id)
            predictor = predictor_class(self.seqvar, self.seqvar_result, self.config)
            seqvar_prediction = predictor.predict()
            # logger.info("Prediction: {}", seqvar_prediction)
            return seqvar_prediction

        elif isinstance(self.strucvar, StrucVar):
            if not self.strucvar:
                logger.error("Failed to resolve the structural variant.")
                return None
            # ====== Setup the data ======
            self.parse_strucvar_data(self.strucvar)

            # ====== Predict ======
            sp = DefaultStrucVarPredictor(self.strucvar, self.strucvar_result, self.config)
            strucvar_prediction = sp.predict()
            logger.info("Prediction: {}", strucvar_prediction)
            return strucvar_prediction

        else:
            logger.info("Structural variants are not supported for ACMG criteria prediction yet.")
            return None
