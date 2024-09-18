.. _index:

===================
Welcome to AutoACMG
===================

Auto-acmg is a comprehensive software suite designed for the automatic prediction of ACMG (American
College of Medical Genetics and Genomics) criteria for the evaluation of sequence variants and copy
number variations (CNVs) in genetic data. Its primary function is to assist geneticists and
researchers by providing a systematic approach to classify genetic variants according to established
guidelines. More information to the guidelines can be found at the :ref:`acmg_seqvars_criteria` and
:ref:`acmg_seqvars_details` sections.

Key Features
------------

- **Variant Classification**: Auto-acmg classifies single nucleotide variants (SNVs) and CNVs based
    on ACMG guidelines, enhancing the accuracy and consistency of genetic interpretations.

- **Implementation of ACMG Criteria**: It supports a range of criteria including PVS1, PS1, PM1,
    PM2, PM4, PM5, PP2, PP3, BA1, BS1, BS2, BP1, BP3, BP4 and BP7, allowing users to
    assess the pathogenicity of variants in a comprehensive manner.

- **VCEP Specifications**: The software is designed to adhere to Variant Curation Expert Panel
    (VCEP) specifications, ensuring that each prediction is tailored to the specific genes and
    conditions studied by expert panels.

- **Extensible and Adaptable**: Auto-acmg can be customized and extended to accommodate new criteria
    and guidelines as genetic science evolves.


Getting Started
---------------

If you are new to auto-acmg, a good starting point is the :ref:`usage` section which provides
detailed instructions on how to use the software effectively.

If you are interested in contributing to auto-acmg, please refer to the :ref:`dev_quickstart`
section.

Also pay attention to the :ref:`auto_acmg_implementation` section to understand how auto-acmg
works internally.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   auto_acmg_implementation
   acmg_seqvars_details
   acmg_seqvars_criteria
   acmg_seqvar_implementation
   acmg_strucvar_implementation
   benchmarking

.. toctree::
    :maxdepth: 2
    :caption: Usage:

    usage

.. toctree::
   :maxdepth: 2
   :caption: Development:

   dev_quickstart

.. toctree::
   :maxdepth: 2
   :caption: Internal API Documentation:

   internal_api
