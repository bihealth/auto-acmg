from typing import List, Optional

from loguru import logger

from src.api.annonars import AnnonarsClient
from src.core.config import Config
from src.defs.auto_acmg import AutoACMGSeqVarResult
from src.defs.auto_pvs1 import PVS1Prediction
from src.defs.seqvar import SeqVar
from src.seqvar.auto_bp7 import AutoBP7
from src.seqvar.auto_pm1 import AutoPM1
from src.seqvar.auto_pm2_ba1_bs1_bs2 import AutoPM2BA1BS1BS2
from src.seqvar.auto_pm4_bp3 import AutoPM4BP3
from src.seqvar.auto_pp2_bp1 import AutoPP2BP1
from src.seqvar.auto_pp3_bp4 import AutoPP3BP4
from src.seqvar.auto_ps1_pm5 import AutoPS1PM5
from src.seqvar.auto_pvs1 import AutoPVS1

#: Pathogenic PVS1 predictions of sequence variants.
PVS1_POSITIVE_SEQVAR_PREDICTIONS = [
    PVS1Prediction.PVS1,
    PVS1Prediction.PVS1_Strong,
    PVS1Prediction.PVS1_Moderate,
    PVS1Prediction.PVS1_Supporting,
]

#: Criteria, which have no automatic prediction yet.
NOT_IMPLEMENTED_CRITERIA: List[str] = [
    "PS2",  # De novo (both parents not tested)
    "PS3",  # Well-established in vitro or in vivo functional studies
    "PS4",  # Increased prevalence in affected individuals vs. controls
    "PM3",  # Recessive disorder (zygosity unknown)
    "PM6",  # Assumed de novo (both parents not tested)
    "PP1",  # Cosegregation with disease in multiple affected family members
    "PP4",  # Patient's phenotype or family history specificity
    "BS3",  # Well-established in vitro or in vivo functional studies
    "BS4",  # Lack of segregation in affected family members
    "BP2",  # Observed in healthy individuals (zygosity unknown)
    "BP5",  # Case with an alternate molecular basis for disease
]


class DefaultSeqVarPredictor(
    AutoPVS1,
    AutoPS1PM5,
    AutoPM1,
    AutoPM2BA1BS1BS2,
    AutoPM4BP3,
    AutoPP2BP1,
    AutoPP3BP4,
    AutoBP7,
):
    def __init__(self, seqvar: SeqVar, result: AutoACMGSeqVarResult, config: Config):
        #: Configuration to use.
        self.config = config or Config()
        #: Sequence variant to predict.
        self.seqvar = seqvar
        #: Annonars client.
        self.annonars_client: AnnonarsClient = AnnonarsClient(
            api_base_url=self.config.api_base_url_annonars
        )
        #: Prediction result.
        self.result = result

    def predict(self) -> Optional[AutoACMGSeqVarResult]:
        # PP5 and BP6 criteria are depricated
        logger.warning("Note, that PP5 and BP6 criteria are depricated and not predicted.")
        # Not implemented criteria
        logger.warning(
            "Some criteria are not implemented yet: {}",
            NOT_IMPLEMENTED_CRITERIA,
        )

        # PVS1
        self.result.criteria.pvs1 = self.predict_pvs1(self.seqvar, self.result.data)

        # PS1 and PM5
        self.result.criteria.ps1, self.result.criteria.pm5 = self.predict_ps1pm5(
            self.seqvar, self.result.data
        )

        # PM1
        self.result.criteria.pm1 = self.predict_pm1(self.seqvar, self.result.data)

        # PM2, BA1, BS1 and BS2
        (
            self.result.criteria.pm2,
            self.result.criteria.ba1,
            self.result.criteria.bs1,
            self.result.criteria.bs2,
        ) = self.predict_pm2ba1bs1bs2(self.seqvar, self.result.data)

        # PM4 and BP3
        self.result.criteria.pm4, self.result.criteria.bp3 = self.predict_pm4bp3(
            self.seqvar, self.result.data
        )

        # PP2 and BP1
        self.result.criteria.pp2, self.result.criteria.bp1 = self.predict_pp2bp1(
            self.seqvar, self.result.data
        )

        # PP3 and BP4
        self.result.criteria.pp3, self.result.criteria.bp4 = self.predict_pp3bp4(
            self.seqvar, self.result.data
        )

        # BP7
        self.result.criteria.bp7 = self.predict_bp7(self.seqvar, self.result.data)

        logger.info("ACMG criteria prediction completed.")
        return self.result
