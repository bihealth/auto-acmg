.. _auto_acmg_implementation:

=============================
AutoACMG Implementation
=============================


AutoACMG operates through a structured process to determine the most appropriate predictor based on
the type and context of the genetic variant being analyzed. This document describes the general
workflow and the mechanisms behind the decisions that auto-acmg makes to predict the classification
of genetic variants.

Workflow Overview
-----------------

1. **Variant Resolution**:
   The first step involves determining the type of genetic variant. If the variant is recognized as
   a sequence variant (seqvar), auto-acmg employs the SeqVar default predictor. For structural
   variants (strucvar), the StrucVar default predictor is utilized, which primarily evaluates the
   PVS1 criteria.

2. **Gene-Specific Evaluation**:
   For seqvars, auto-acmg checks if the variant is associated with any gene that has a specific
   Variant Curation Expert Panel (VCEP) implementation. If a specific VCEP is applicable, auto-acmg
   will switch to the corresponding VCEP predictor for that gene to leverage expert-specified
   evaluation criteria and thresholds.

   - **SeqVar Implementation**: If no specific VCEP is applicable, the default SeqVar predictor is
     used. The implementation details of the default SeqVar predictor are described in
     :ref:`acmg_seqvar_implementation`.

3. **Structural Variants**:
   For structural variants, the default predictor primarily evaluates the PVS1 criteria. You can
   find the implementation of the default StrucVar predictor in
   :ref:`acmg_strucvar_implementation`.

VCEP-Specific Overrides
-----------------------

Each VCEP-specific predictor may override certain logic or thresholds of the Default Predictor:

- **Threshold Adjustments**: Some VCEP implementations may just adjust thresholds for criteria like
    PP3 or BP4 and don't override the prediction logic.

- **Prediction Logic Changes**: Other VCEPs may override the entire logic for certain criteria. This
    can include changing how predictions are made based on the variant's effects on protein function,
    splicing, or other molecular mechanisms.


.. note::
    For the **RYR1 gene**, there are two VCEPs: *Congenital Myopathies VCEP* and
    *Malignant Hyperthermia Susceptibility*. In AutoACMG, we only consider the
    **Malignant Hyperthermia Susceptibility VCEP** for RYR1 and avoid the Congenital Myopathies
    VCEP. For more details on the VCEP implementations, please refer to the source code in the
    ``src/vcep`` directory of our GitHub repository.

.. note::
    For **von Willebrand Disease**, there are two different rulesets. We implement the combination
    of both rulesets. For more details on the VCEP implementations, please refer to the source code
    in the ``src/vcep`` directory of our GitHub repository.


For a detailed overview of how specific VCEP predictors modify or extend the default behavior, you
can review the source code on the GitHub repository under ``src/vcep`` and for the details refer to
the `ClinGen VCEP Criteria Specification Registry <https://cspec.genome.network/cspec/ui/svi/>`__.

Link to Source Code
-------------------

Further details and the actual implementation of both default and VCEP-specific predictors can be
accessed through our `GitHub <https://github.com/auto-acmg/src/vcep>`_ repository.

.. note::
   The implementation specifics for each VCEP can vary significantly based on the gene and
   associated conditions. It is recommended to consult the individual VCEP documentation and the
   source code for precise information on how each gene-specific predictor is implemented.
