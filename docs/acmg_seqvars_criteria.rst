.. _acmg_seqvars_criteria:

==============================
ACMG Sequence Variant Criteria
==============================

This section describes the criteria as implemented for the automated classification of ACMG sequence variants.

A large focus is given towards the transparency to the user.
The user is always given a report of the evidence for and against a given criterion.
Even if not further described, the reason for skipping a criterion is always reported (e.g., BA1 disables all other benign criteria, or if a criterion is only valid for missense variants).

.. _acmg_seqvars_criteria-inheritance:

-------------------
Mode of Inheritance
-------------------

The mode of inheritance for a is derived from the following list of sources.
The sources are iterated in the order given below and the first one with a match is used for deriving mode of inheritance.

1. **NHGRI CGD**
   The `National Human Genome Research Institute Clinical Genomic Database <https://research.nhgri.nih.gov/CGD/>`__
2. **ClinGen Disease Validity** `Website <https://clinicalgenome.org/curation-activities/gene-disease-validity/>`__
3. **EBI gene2phenotype** `Website <https://www.ebi.ac.uk/gene2phenotype>`__
4. **ClinGen GenCC** `Website <https://thegencc.org/>`__
5. **Genomics England PanelApp** `Web-App <https://panelapp.genomicsengland.co.uk/>`__
6. **Domino** is a machine learning method for prediction of mode of inheritance and described in `PMID:28985496 <https://pubmed.ncbi.nlm.nih.gov/28985496/>`__.
   We use the thresholds from `PMID:30376034 <https://pubmed.ncbi.nlm.nih.gov/30376034/>`__.

   - dominant: if score >= 0.5934
   - recessive: if score <= 0.3422
   - dominant/recessive: if score is between 0.3422 and 0.5934

.. _acmg_seqvars_criteria-inheritance-literature:

Literature
==========

Refer to :ref:`ACMG details<acmg_seqvars_details-references>` for the literature.


.. _acmg_seqvars_criteria-frequency:

----------------
Allele Frequency
----------------

The databases gnomAD exomes and genomes will be used to derived allele frequencies.
We trust the gnomAD quality control pipeline and filter solely above the gnomAD variant annotation.
Frequencies are considered valid if:

- allele number is >=2000
- gnomAD quality filter is PASS

The criteria BA1 and BS1 will ignore the following "bottlenecked" populations:
Amish (ami), Ashkenazi Jewish (asj), European Finnish (fin), Middle Eastern (mid), and "Remaining Individuals (rmi) groups".

.. _acmg_seqvars_criteria-calibration:

-----------
Calibration
-----------

We currently do not have our own calibration yet.
This is a big TODO once we have all the data.

--------------------
Transcript Selection
--------------------

We use the MANE Plus Clinical set of transcripts for the genes that have such information.
Otherwise, we use the longest transcript as the MANE transcript and consider all alternate RefSeq transcripts.

MANE annotation is not available for GRCh37 and will not be provided by ENSEMBL/NCBI.
We thus map the MANE (and MANE Plus Clinical) transcripts from GRCh38 by using the latest version of the transcript with the same identifier.
If mapping either MANE or any of the MANE Plus Clinical transcripts fails then we fall back to the longest transcript rule.

.. _acmg_seqvars_criteria-patho-predictions:

-------------------------
Pathogenicity Predictions
-------------------------

We currently use the thresholds from :footcite:t:`pejaver:2022` only.

Once we have our own calibration, we can extend our predictions to novel tools such as AlphaMissense.

.. _acmg_seqvars_mods:

------------------------------
Code Modification Nomenclature
------------------------------

In accordance with ClinGen Sequence Variant Interpretation Work Group (2017), modification codes are used.
That is, for a criterion ``${CRIT}``, the modification codes ``${CRIT}_Supporting``, ``${CRIT}_Moderate``, ``${CRIT}_Strong``, ``${CRIT}_VeryStrong``.

.. _acmg_seqvars_criteria-rules:

--------------------
Rules / Point System
--------------------

We consider three rule systems:

- The original ACMG 2015 rules
- The ACGS 2020 rules
- The 2020 Point system described by Tavtigian et al. (2020)

ACMG 2015 Rules
===============

The following rules have been defined in Richards et al. (2015).

Pathogenic
----------

If one of the following criteria 1-3 are fulfilled.

1. 1 very strong (PVS1) AND one of the following
    a. >=1 strong (PS1-PS4)
    b. >=2 moderate (PM1-PM6)
    c. >=1 moderate (PM1-PM6) AND >=1 supporting (PP1-PP5)
    d. >=2 supporting (PP1-PP5)
2. >=2 strong (PS1-PS4)
3. 1 strong (PS1-PS4) AND
    a. >=3 moderate (PM1-PM6)
    b. 2 moderate (PM1-PM6) AND >=2 supporting (PP1-PP5)
    c. 1 moderate (PM1-PM6) AND >=4 supporting (PP1-PP5)

Likely Pathogenic
-----------------

If one of the following criteria 1-7 are fulfilled.

1. 1 very strong (PVS1) AND 1 moderate (PM1-PM6)
2. 1 strong (Ps1-PS4) AND 1-2 moderate (PM1-PM6)
3. 1 strong (PS1-PS4) AND >=2 supporting (PP1-PP5)
4. >=3 moderate (PM1-PM6)
5. 2 moderate (PM1-PM6) AND >=2 supporting (PP1-PP5)
6. 1 moderate (PM1-PM6) AND >=4 supporting (PP1-PP5)

Benign
------

If one of the following criteria 1-2 are fulfilled.

1. 1 standalone (BA1)
2. >=2 strong (BS1-BS4)

Likely Benign
-------------

If one of the following criteria 1-2 are fulfilled.

1. 1 strong (BS1-BS4) AND 1 supporting (BP1-BP7)
2. >=2 supporting (BP1-BP7)

Uncertain Significance
----------------------

If if one of the following criteria 1-2 are fulfilled.

1. Other criteria shown above are not met
2. the criteria for benign and pathogenic are contradictory

ACGS 2020 Rules
===============

The following is a refinement of the rules above set by the Ellard et al. (2020).

Pathogenic
----------

1. 1 very strong (PVS) AND one of the following
    a. >=1 strong
    b. >=1 moderate
    c. >=2 supporting
2. >=3 strong
3. 2 strong AND one of the following
    a. >=1 moderate
    b. >=2 supporting
4. 1 strong AND one of the following
    a. >=3 moderate
    b. >=2 moderate AND >=2 supporting
    c. >=1 moderate AND >=4 supporting

Likely Pathogenic
-----------------


1. >=2 strong
2. 1 strong AND one of teh following
    a. 1-2 moderate OR
    b. >=2 supporting
3. >=3 moderate OR
    a. 2 moderate AND >=2 supporting
    b. 1 modereate AND >=4 supporting

Tavtigian et al. (2020) Rules
=============================

Alternatively, Tavtigian et al. (2020) formulated the rules as an integer point system.

Table 2 from this manuscript gives point values:

.. list-table:: Points per proband

    * - evidence strength
      - points pathogenic
      - points benign
    * - indeterminate
      - 0
      - 0
    * - supporting
      - 1
      - -1
    * - moderate
      - 2
      - -2
    * - strong
      - 4
      - -4
    * - very strong
      - 8
      - -8

The point-based variant classification categories are then given in their Table 3:

.. list-table:: Categories and point ranges

    * - category
      - point ranges
    * - pathogenic
      - >= 10
    * - likely pathogenic
      - 6 to 9
    * - uncertain significance
      - 0 to 5
    * - likely benign
      - -1 to -6
    * - benign
      - <= -7

--------
Criteria
--------

.. _acmg_seqvars_criteria-pvs1:

PVS1 (null variant)
===================

.. note::

    For full information visit the dedicated section of :ref:`PVS1 Algorithm <acmg_pvs1_algorithm>`.

Original Definition
-------------------

    Null variant (nonsense, frameshift, canonical +/-1 or 2 splice sites, initiation codon, single or multi-exon deletion) in a gene where loss of function (LOF) is a known mechanism of disease.

    Caveats:

    - Beware of genes where LOF is not a known disease mechanism (e.g. GFAP, MYH7)
    - Use caution interpreting LOF variants at the extreme 3' end of a gene
    - Use caution with splice variants that are predicted to lead to exon skipping but leave the remainder of the protein intact
    - Use caution in the presence of multiple transcripts

    -- Richards et al. (2015); Table 4

Preconditions / Precomputations
-------------------------------

- Criterion establishes whether LoF is a known mechanism of disease:
    - If at least 2 LoF variants are reported in ClinVar with two or more stars then this criterion is triggered.
    - If the gnomAD LOF Observed/Expected is less than 0.7555 then this criterion is triggered.
- Criterion establishes whether a stop_gain variant introduced nonsense mediated decay (NMD) consistent with Abou Youn et al. (2018) and the VEP NMD plugin.
    - If the variant is on chrMT then it cannot be NMD.
    - If the variant is not_stop gain then then it cannot be NMD, else:
    - If the variant is in the last exon of the transcript then it is predicted to escape NMD.
    - If the variant falls 50bp upstream of the penuultimate (second to the last) exon then it is predicted to escape NMD.
    - If the variant falls int the first 100 coding bases in teh transcript then it is predicted to escape NMD.
    - If the variant is in an intronless transcript, meaning only one exon exists in the transcript, then it is predicted to escape NMD.
    - Else, the variant is predicted to be NMD.
- The MANE Plus Clinical transcripts are used for "biologically relevant transcripts" in this criterion.

Implemented Criterion
---------------------

While the original description is somewhat vague, the specification in Abou Tayoun et al. (2018) is more precise but complex to implement.
We plan to implement it as closely as possible.

Refer to the dedicated section of `PVS1 Algorithm <acmg_pvs1>`__ for the algorithm.

Literature
----------

- Richards et al. (2015) describes the original criterion.
- Abou Tayoun et al. (2018) describe refined criteria for PVS1.
- McCormick et al. (2020) describe the ACMG criteria for chrMT variants.
- The following are from the VEP NMD plugin:
    - Identifying Genes Whose Mutant Transcripts Cause Dominant Disease Traits by Potential Gain-of-Function Alleles (Coban-Akdemir, 2018)
    - The criteria and impact of nonsense-mediated mRNA decay in human cancers (Lindeboom, 2016)

User Report
-----------

The following information is reported to the user:

- The evidence for / against LoF as disease mechanism.
- Whether NMD and NMD escape is predicted for this variant and the reason.
- The use of MANE Plus Clinical or alternate transcripts for locating alternate start codons.
- Further information of interest from the Abou Tayoun et al. (2018) decision tree.

Caveats
-------

- We use the thresholds from `PMID:30376034 <https://pubmed.ncbi.nlm.nih.gov/30376034/>`__ but should reconsider, e.g., switching to LOEUF here with our own thresholds.
- This is currently not implementing the full criteria set from Abou Tayoun et al. (2018).

Notes
-----

- If this criterion is triggered then PP3 and PM4 will be disabled.

.. _acmg_seqvars_criteria-ps1:

PS1 (same amino acid)
=====================

Original Definition
-------------------

    Same amino acid change as a previously established pathogenic variant regardless of nucleotide change.

    Caveat: Beware of changes that impact splicing rather than at the amino acid/protein level.

    -- Richards et al. (2015); Table 4

Preconditions / Precomputations
-------------------------------

- If the variant is not a missense variant then this criterion is skipped.

Implemented Criterion
---------------------

- Consider all equivalent missense variants in ClinVar.
- If at least one of the variant then this criterion is triggered.
    - If the variant has zero stars in ClinVar then we report PS1_Supporting only
    - If the variant has only one star in ClinVar then we report PS1_Moderate only
    - If the variant has two stars in ClinVar then we report PS1
    - If the variant has three stars or above in ClinVar then we report PS1_VeryStrong

User Report
-----------

- The selected variant in ClinVar and with assessment its star status with accession.
- All alternate variants in Clinvar with assessments and star status with accessions.

Literature
----------

N/A

Caveats
-------

- The wording of "established pathogenic" variant is not clear so we use the steps from above.
- Note that this also depends on disease match which the user must confirm manually.

.. _acmg_seqvars_criteria-ps2:

PS2 (confirmed *de novo*)
=========================

No automation has been implemented.

Original Definition
-------------------

    De novo (both maternity and paternity confirmed) in a patient with the disease and no family history

    Note: Confirmation of paternity only is insufficient.
    Egg donation, surrogate motherhood, errors in embryo transfer, etc. can contribute to non-maternity.

    -- Richards et al. (2015); Table 4

.. _acmg_seqvars_criteria-ps3:

PS3 (functional studies)
========================

No automation has been implemented.

Original Definition
-------------------

    Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product.

    Note: Functional studies that have been validated and shown to be reproducible and robust in a clinical diagnostic laboratory setting are considered the most well-established.

    -- Richards et al. (2015); Table 4

.. _acmg_seqvars_criteria-ps4:

PS4 (prevalence)
================

No automation has been implemented.

Original Definition
-------------------

    The prevalence of the variant in affected individuals is significantly increased compared to the prevalence in controls

    Note 1: Relative risk (RR) or odds ratio (OR), as obtained from case-control studies, is >5.0 and the confidence interval around the estimate of RR or OR does not include 1.0. See manuscript for detailed guidance.

	Note 2: In instances of very rare variants where case-control studies may not reach statistical significance, the prior observation of the variant in multiple unrelated patients with the same phenotype, and its absence in controls, may be used as moderate level of evidence.

    -- Richards et al. (2015); Table 4

.. _acmg_seqvars_criteria-pm1:

PM1 (hotspot)
=============

Original Definition
-------------------

    Located in a mutational hot spot and/or critical and well-established functional domain (e.g. active site of an enzyme) without benign variation.

    -- Richards et al. (2015); Table 4

Preconditions / Precomputations
-------------------------------

- If the variant is on chrMT then this criterion is skipped according to McCormick et al. (2020).

Implemented Criterion
---------------------

- If the variant is within a hotspot (at least 4 pathogenic missense/in-frame variants within 25bp radius) then this criterion is triggered.
- If the variant is within an annotated UniProt domain and the domain contains at least 2 pathogenic variants then this criterion is triggered.

User Report
-----------

- The hotspot region definition and the number of pathogenic variants in the region.

Literature
----------

- McCormick et al. (2020) describe the ACMG criteria for chrMT variants.

Caveats
-------

- We currently use the threshold from `PMID:30376034 <https://pubmed.ncbi.nlm.nih.gov/30376034/>`__ and are lacking our own calibration.

.. _acmg_seqvars_criteria-pm3:

PM3 (recessive in *trans*)
==========================

No automation has been implemented.

Original Definition
-------------------

    For recessive disorders, detected in trans with a pathogenic variant.

	Note: This requires testing of parents (or offspring) to determine phase.

    -- Richards et al. (2015); Table 4

.. _acmg_seqvars_criteria-pm4:

PM4 (protein length)
====================

Original Definition
-------------------

    Protein length changes due to in-frame deletions/insertions in a non-repeat region or stop-loss variants.

    -- Richards et al. (2015); Table 4

Preconditions / Precomputations
-------------------------------

- If PVS1 was triggered then this criterion is skipped to avoid double counting.
- If the variant is not an in-frame indel and not a stop-loss variant then this criterion is skipped.

Implemented Criterion
---------------------

- If the variant is an in-frame indel
    - If the variant is inside a repeat masked region then it is skipped
    - If the variant is inside a repeat as annotated by UniProt then it is skipped
    - Otherwise, this criterion is triggered.
- If the variant is a stop-loss variant then this criterion is triggered.

User Report
-----------

- Any reasons for skipping in repeat regions.
- The transcript identifier.

Literature
----------

N/A

Caveats
-------

- Richards et al. (2015) state that the size of the indel and amount of change in amino acids should influence the classification.
  We currently do not have this implemented.

.. _acmg_seqvars_criteria-pm5:

PM5 (overlapping missense)
==========================

Original Definition
-------------------

    Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before.

    Caveat: Beware of changes that impact splicing rather than at the amino acid/protein level.

    -- Richards et al. (2015); Table 4

Preconditions / Precomputations
-------------------------------

- If the variant is on a nuclear chromosome
    - If it is not a missense variant then this criterion is skipped.
- If the variant is on chrMT and not missense and not on a tRNA gene then this criterion is skipped.

Implemented Criterion
---------------------

- If the variant is on a nuclear chromosome:
    - If the variant is at the same position as a pathogenic missense variant then this criterion is triggered.
- If the variant is on chrMT:
    - If the variant is a missense variant and at the same position as a pathogenic one then the criterion is triggered.
    - If the variant is on a tRNA gene and at the same position as a pathogenic one then the criterion is triggered as PM5_Supporting.

User Report
-----------

- The overlapping variant used for criterion.
- Any alternative overlapping variants not chosen.

Literature
----------

- Richards et al. (2018) describes the criterion for nuclear chromosomes.
- McCormick et al. (2020) describes the criterion for chrMT.

Caveats
-------

N/A

.. _acmg_seqvars_criteria-pm6:

PM6 (assumed *de novo*)
=======================

No automation has been implemented.

Original Definition
-------------------

    Assumed de novo, but without confirmation of paternity and maternity.

    -- Richards et al. (2015); Table 4

.. _acmg_seqvars_criteria-pm2:

PM2_Supporting (absent from controls)
=====================================

Original Definition
-------------------

    Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes or ExAC.

    -- Richards et al. (2015); Table 4

Preconditions / Precomputations
-------------------------------

- Determine :ref:`acmg_seqvars_criteria-inheritance` for the gene.
- Determine :ref:`acmg_seqvars_criteria-frequency`.
- If the allele frequency is invalid then this criterion is skipped.

Implemented Criterion
---------------------

- If the variant is on a nuclear chromosome:
    - If the gene is marked as recessive or X-linked:
        - If the variant allele count is <=4 then this criterion is triggered.
    - If the gene is marked as dominant:
        - If the homozygous allele count is <=1 then this criterion is triggered.
        - If the allele frequency is less than 0.0001 then this criterion is triggered.
- If the variant is on chrMT:
    If the variant frequency is below 0.00002=0.002%=1/50,000 then this criterion is triggered.

User Report
-----------

- The values and thresholds used by the criterion even if failed.

Literature
----------

- Richards et al. (2015) describes the original criterion.
- ClinGen Sequence Variant Interpretation Work Group (2020): SVI Recommendation for Absence/Rarity (PM2) - Version 1.0 describes the downgrade to supporting.
- McCormick et al. (2020) describe the ACMG criteria for chrMT variants.

Caveats
-------

- We currently use the threshold from `PMID:30376034 <https://pubmed.ncbi.nlm.nih.gov/30376034/>`__ and are lacking our own calibration.
- This criterion has been downgraded by default to supporting from strong in accordance to ClinGen Sequence Variant Interpretation Work Group (2020): *SVI Recommendation for Absence/Rarity (PM2) - Version 1.0*

.. _acmg_seqvars_criteria-pp1:

PP1 (cosegregation)
===================

No automation has been implemented.

.. _acmg_seqvars_criteria-pp2:

PP2 (missense)
==============

Original Definition
-------------------

    Missense variant in a gene that has a low rate of benign missense variation and where missense variants are a common mechanism of disease.

    -- Richards et al. (2015); Table 4

Preconditions / Precomputations
-------------------------------

- If the variant is on chrMT then this criterion is skipped according to McCormick et al. (2020).
- If the variant is not a missense variant then this criterion is skipped.

Implemented Criterion
---------------------

- If the ratio of pathogenic missense variants over all non-VUS missense variants is greater than 0.808 then this criterion is triggered.

User Report
-----------

- Report the ratio of pathogenic missense variants over all non-VUS missense variants.

Literature
----------

- McCormick et al. (2020) describe the ACMG criteria for chrMT variants.

Caveats
-------

- We currently use the threshold from `PMID:30376034 <https://pubmed.ncbi.nlm.nih.gov/30376034/>`__ and are lacking our own calibration.

Notes
-----

- This criterion is similar to :ref:`acmg_seqvars_criteria-bp1`

.. _acmg_seqvars_criteria-pp3:

PP3 (*in silico* predictions)
=============================

Original Definition
-------------------

    Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact, etc).

    Caveats:

    - As many in silico algorithms use the same or very similar input for their predictions, each algorithm should not be counted as an independent criterion.
    - PP3 can be used only once in any evaluation of a variant.

    -- Richards et al. (2015); Table 4

Preconditions / Precomputations
-------------------------------

- If the criterion PVS1 was triggered then this criterion is skipped.
- If the variant is on chrMT then it is skipped, as we don't have calibration for chrMT yet.
- If the variant is not found in dbNSFP or CADD precomputed scores then it is skipped as we don't have calibration for chrMT yet.

Implemented Criterion
---------------------

An initial prediction is fist done using the general purpose pathogenicity predictors.

- If we have a score from the following, then the prediction is used (in descending order of priority):
    - REVEL, MutPred2, CADD, BayesDel, VEST4, ..., PhyloP
    - we will use the modifiers from :footcite:t:`pejaver:2022`
- If predictions are missing then then PhyloP of the position of the variant is used as a fallback.

Then, for splicing the following is done.

- If a SpliceAI prediction is performed then it is interpreted according to :footcite:t:`walker:2023`.

The highest-scoring variant is used for the final prediction.

User Report
-----------

- The scores and predictions from the predictors.

Caveats
-------

- As described in :ref:`acmg_seqvars_criteria-patho-predictions`, we are currently limited to the precomputed threshold from the literature.
  This hinders us in adopting AlphaMissense effectively, for example.
- We need to compute accuracy to rank the implemented methods.
- We need our own calibration for chrMT.

Notes
-----

- This criterion is similar to :ref:`acmg_seqvars_criteria-bp4`

.. _acmg_seqvars_criteria-pp4:

PP4 (monogenetic)
=================

No automation has been implemented.

.. _acmg_seqvars_criteria-ba1:

BA1 (5% frequency)
==================

Original Definition
-------------------

    Allele frequency is >5% in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium

    -- Richards et al. (2015); Table 4

Preconditions / Precomputations
-------------------------------

- The variant is absent from the exception list from Ghosh et al. (2018).
  If the variant is present on this list, then this criterion is skipped.

Implemented Criterion
---------------------

- If the variant is nuclear (not on chrMT)
    - If the allele frequency is above 0.05 in gnomAD global population then this criterion is triggered.
- else (the variant is on chrMT)
    - If the allele frequency is above 0.01 in gnomAD-mtDNA global population then this criterion is triggered.

User Report
-----------

- The variant frequency.

Literature
----------

- Richards et al. (2015) describes the 5% allele frequency threshold.
- Ghosh et al. (2018) introduce the exception list and ClinGen maintains it.
- McCormick et al. (2020) describe the 1% allele frequency threshold as appropriate for chrMT variants.

Caveats
-------

- The exception *"However, there must be no additional conflicting evidence to support pathogenicity, such as a novel occurrence in a certain haplogroup" from McCormick et al. (2020)* is not implemented yet.

.. _acmg_seqvars_criteria-bs1:

BS1 (expected frequency)
========================

Original Definition
-------------------

    Allele frequency greater than expected for disorder.

    -- Richards et al. (2015); Table 4

Preconditions / Precomputations
-------------------------------

- Determine :ref:`acmg_seqvars_criteria-frequency`.
- If the allele frequency is invalid then this criterion is skipped.

Implemented Criterion
---------------------

- If the variant is on a nuclear chromosome and the user provided a maximal credible population frequency:
    - If the FAF from gnomAD is above the maximal credible population frequency then this criterion is triggered.
- If the variant is on chrMT:
    - If the population frequency is above 0.5% then this criterion is triggered in accordance to McCormick et al. (2020).

User Report
-----------

- The variant frequency and again the user specified maximal credible population frequency for nuclear variants.
- The variant frequency and the 0.5% threshold for chrMT variants.

Literature
----------

- Richards et al. (2015) describes the original criterion without thresholds.
- Gudmundsson et al. (2022) describe the FAF threshold provided by gnomAD.
- McCormick et al. (2020) describe the ACMG criteria for chrMT variants.

.. _acmg_seqvars_criteria-bs2:

BS2 (healthy adult)
===================

Original Definition
-------------------

    Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked (hemizygous) disorder, with full penetrance expected at an early age.

    -- Richards et al. (2015); Table 4

Preconditions / Precomputations
-------------------------------

- If the criterion BA1 triggered then this criterion is skipped.
- Determine :ref:`acmg_seqvars_criteria-inheritance` for the gene.
- Determine :ref:`acmg_seqvars_criteria-frequency`.
- If the allele frequency is invalid then this criterion is skipped.
- If the criterion BA1 was triggered then this criterion is skipped.

Implemented Criterion
---------------------

- If the gene is marked as recessive or X-linked:
    - If the variant allele count is above 2 then this criterion is triggered.
- If the gene is marked as dominant:
    - If the variant allele count is above 5 then this criterion is triggered.

User Report
-----------

- The variant frequency and the threshold used.

Literature
----------

- Chen et al. (2022), Karczewski et al. (2020), etc. describe gnomAD.
- The modes of inheritance for the genes are taken from different sources as described in :ref:`acmg_seqvars_criteria-inheritance`.

Caveats
-------

- The conditions of "full penetrance" and "expected at an early age" need to be checked by the user.

Notes
-----

- Genes can be marked as both recessive and dominant.
- We use the thresholds from `PMID:30376034 <https://pubmed.ncbi.nlm.nih.gov/30376034/>`__.

.. _acmg_seqvars_criteria-bs3:

BS3 (functional studies)
========================

No automation has been implemented.

Original Definition
-------------------

    Well-established in vitro or in vivo functional studies shows no damaging effect on protein function or splicing.

    -- Richards et al. (2015); Table 4

.. _acmg_seqvars_criteria-bs4:

BS4 (lack of segregation)
=========================

No automation has been implemented.

Original Definition
-------------------

    Lack of segregation in affected members of a family

    Caveats:

    - The presence of phenocopies for common phenotypes (i.e. cancer, epilepsy) can mimic lack of segregation among affected individuals.
    - Also, families may have more than one pathogenic variant contributing to an autosomal dominant disorder, further confounding an apparent lack of segregation.

    -- Richards et al. (2015); Table 4

.. _acmg_seqvars_criteria-bp1:

BP1 (missense)
==============

Original Definition
-------------------

    Missense variant in a gene for which primarily truncating variants are known to cause disease

    -- Richards et al. (2015); Table 4

Preconditions / Precomputations
-------------------------------

- If the criterion BA1 triggered then this criterion is skipped.
- If the variant is on chrMT then this criterion is skipped according to McCormick et al. (2020).
- If the variant is not a missense variant then this criterion is skipped.

Implemented Criterion
---------------------

- If the ratio of benign missense variants over all non-VUS missense variants is greater than 0.569 then this criterion is triggered.

User Report
-----------

- Report the ratio of benign missense variants over all non-VUS missense variants together with threshold.

Literature
----------

- McCormick et al. (2020) describe the ACMG criteria for chrMT variants.

Caveats
-------

- We currently use the threshold from `PMID:30376034 <https://pubmed.ncbi.nlm.nih.gov/30376034/>`__ and are lacking our own calibration.

Notes
-----

- This criterion is similar to :ref:`acmg_seqvars_criteria-pp2`

.. _acmg_seqvars_criteria-bp2:

BP2 (recessive in *trans*)
==========================

No automation has been implemented.

Original Definition
-------------------

    Observed in trans with a pathogenic variant for a fully penetrant dominant gene/disorder; or observed in cis with a pathogenic variant in any inheritance pattern

    -- Richards et al. (2015); Table 4

.. _acmg_seqvars_criteria-bp3:

BP3 (in-frame repetitive)
=========================

.. note::

    - We do not have proper Uniprot data yet (domain / repeat)
    - Similar to repeat masker.
    - Probably same for phylop100way?

Original Definition
-------------------

    In-frame deletions/insertions in a repetitive region without a known function.

    -- Richards et al. (2015); Table 4

Preconditions / Precomputations
-------------------------------

- If the criterion BA1 triggered then this criterion is skipped.
- If the variant is on chrMT then this criterion is skipped.

Implemented Criterion
---------------------

- If the variant is in a known functional domain according to UniProt then this criterion is skipped.
- If the variant is in a repeat region according to UniProt repeat annotation genome repeat masker then this criterion is skipped.
- If the variant is in a region of low conservation (PhyloP100Way less than 3.58, same as `PMID:30376034 <https://pubmed.ncbi.nlm.nih.gov/30376034/>`__) then this criterion is skipped.
- If all conditions above fail then this criterion is triggered.

User Report
-----------

- The variant position and the reason for triggering or skipping.

Literature
----------

- McCormick et al. (2020) describe the ACMG criteria for chrMT variants.

Caveats
-------

- We currently use the conservation threshold from `PMID:30376034 <https://pubmed.ncbi.nlm.nih.gov/30376034/>`__ and are lacking our own calibration.
- Different from `PMID:30376034 <https://pubmed.ncbi.nlm.nih.gov/30376034/>`__, we do not check whether there are known pathogenic variants in the region.

.. _acmg_seqvars_criteria-bp4:

BP4 (*in silico* predictions)
=============================

.. note::

    - we have not implemented MitoTip or MitImpact yet
    - we are lacking phylop scores yet
    - we don't have live CADD scores yet

Original Definition
-------------------

    Multiple lines of computational evidence suggest no impact on gene or gene product (conservation, evolutionary, splicing impact, etc).

    Caveat: As many in silico algorithms use the same or very similar input for their predictions, each algorithm cannot be counted as an independent criterion.
    BP4 can be used only once in any evaluation of a variant.

    -- Richards et al. (2015); Table 4

Preconditions / Precomputations
-------------------------------

- If the criterion BA1 triggered then this criterion is skipped.
- If the variant is on chrMT then it is skipped, as we don't have calibration for chrMT yet.
- If the variant is not found in dbNSFP or CADD precomputed scores then it is skipped as we don't have calibration for chrMT yet.

Implemented Criterion
---------------------

See :ref:`acmg_seqvars_criteria-pp3` for details.

User Report
-----------

See :ref:`acmg_seqvars_criteria-pp3` for details.

Literature
----------

See :ref:`acmg_seqvars_criteria-pp3` for details.

Caveats
-------

See :ref:`acmg_seqvars_criteria-pp3` for details.

Notes
-----

- This criterion is similar to :ref:`acmg_seqvars_criteria-pp3`

.. _acmg_seqvars_criteria-bp5:

BP5 (found in solved)
=====================

No automation has been implemented.

Original Definition
-------------------

    Variant found in a case with an alternate molecular basis for disease.

    -- Richards et al. (2015); Table 4

.. _acmg_seqvars_criteria-bp7:

BP7 (synonymous)
================

Original Definition
-------------------

    A synonymous (silent) variant for which splicing prediction algorithms predict no impact to the splice consensus sequence nor the creation of a new splice site AND the nucleotide is not highly conserved.

    -- Richards et al. (2015); Table 4

Preconditions / Precomputations
-------------------------------

- If the variant is on chrMT then this criterion is skipped according to McCormick et al. (2020).

Implemented Criterion
---------------------

- If there is a pathogenic variant +/- 2bp of the position in ClinVar then the criterion is skipped.
- If the variant is closer than 2bp to a splice site then the criterion is skipped.
- If the variant is not predicted to alter the splice site using SpliceAI then the criterion is triggered.

User Report
-----------

- The variant position and the reason for triggering or skipping.

Literature
----------

- McCormick et al. (2020) describe the ACMG criteria for chrMT variants.

Caveats
-------

N/A

Notes
-----

- We use the thresholds from `PMID:30376034 <https://pubmed.ncbi.nlm.nih.gov/30376034/>`__.