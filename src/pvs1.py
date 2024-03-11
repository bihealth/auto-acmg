"""Implementations of the PVS1 algorithm."""

import random
import string

from src.seqvar import SeqVar


class PVS1:
    """PVS1 algorithm implementation."""

    def __init__(self, seqvar: SeqVar):
        self.seqvar = seqvar
        self.id = self.id_generator()
        self.vep_input = f"/tmp/vep.{self.id}.vcf"
        self.vep_output = f"/tmp/vep.{self.id}.tab"

        self.vep_variation = "na"
        self.vep_symbol = "na"
        self.vep_trans = "na"
        self.vep_canonical = "na"
        self.vep_pick = "na"
        self.vep_consequence = "na"
        self.hgvs_c = "na"
        self.hgvs_p = "na"
        self.hgvs_g = "na"
        self.vep_exon = "na"
        self.vep_intron = "na"

    @staticmethod
    def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
        """Generates a unique identifier for the VEP."""
        return "".join(random.choice(chars) for _ in range(size))
