.. _usage:

=====
Usage
=====

The following section describe the usage of the AutoACMG package.

.. _installation:

Installation
------------
TODO

.. _configuration:

Configuration
-------------

To run the AutoACMG package, you'll need to create config. Example one can be found in
`.env.dev`. Change the settings to match your environment.

.. code-block:: bash

    cp .env.dev .env

.. _running:

Running
-------

To run the AutoACMG package, you'll need to run the following command:

.. code-block:: bash

    make run VAR=<variant_name> GR=<genome_release>

Where:
- ``<variant_name>`` is the name of the variant you want to analyze.
- ``<genome_release>`` is the genome release you want to use.

Alternatively you can use it inline:

.. code-block:: python

    from src.auto_acmg import AutoACMG

    autoacmg = AutoACMG(variant_name, genome_release)
    prediction = autoacmg.predict()

Note that the ``predict`` method will return a dictionary with the prediction.
This dictionary contains all ACMG criteria with their:
- ``name``: the name of the ACMG criteria.
- ``prediction``: the prediction of the ACMG criteria.
- ``summary``: the details of the prediction.
- ``description``: the description of the ACMG criteria.

The valid values for the ``prediction`` are:
- ``NotSet``: the ACMG criteria was not set (also due to errors).
- ``Positive``: the ACMG criteria is positive.
- ``Negative``: the ACMG criteria is negative.
- ``NotAutomated``: the ACMG criteria is not automated (e.g. due to the need of manual curation).
- ``NotApplicable``: the ACMG criteria is not applicable.
- ``Depricated``: the ACMG criteria is depricated (e.g. PP5 or BP6).

The evidence level and the prediction path for the PVS1 ACMG criteria is returned at the summary
field in regards to "PVS1 Paper". The values represent the following:

For evidence level:
- ``NotSet``: the evidence level was not set (also due to errors).
- ``PVS1``: the strongest evidence level.
- ``PVS1_Strong``: the strong evidence level.
- ``PVS1_Moderate``: the moderate evidence level.
- ``PVS1_Supporting``: the supporting evidence level.
- ``NotPVS1``: the evidence level is not PVS1.
- ``UnsupportedConsequence``: PVS1 is not applicable due to unsupported consequence.

For prediction path:

.. code-block:: python

    #: Mapping of PVS1 prediction path to description for sequence variant
    PVS1PredictionPathMapping: Dict[
        Union[PVS1PredictionSeqVarPath, PVS1PredictionStrucVarPath], str
    ] = {
        PVS1PredictionSeqVarPath.NotSet: "Not Set",
        PVS1PredictionSeqVarPath.PTEN: "Special guideline for PTEN -> Predicted to undergo NMD",
        PVS1PredictionSeqVarPath.NF1: (
            "Predicted to undergo NMD -> Exon is present in biologically-relevant transcript(s)"
        ),
        PVS1PredictionSeqVarPath.NF2: (
            "Predicted to undergo NMD -> Exon is absent from biologically-relevant transcript(s)"
        ),
        PVS1PredictionSeqVarPath.NF3: (
            "Not predicted to undergo NMD -> "
            "Truncated/altered region is critical to protein function"
        ),
        PVS1PredictionSeqVarPath.NF4: (
            "Not predicted to undergo NMD -> "
            "Role of region in protein function is unknown -> "
            "LoF variants in this exon are frequent in the general population and/or "
            "exon is absent from biologically-relevant transcript(s)"
        ),
        PVS1PredictionSeqVarPath.NF5: (
            "Not predicted to undergo NMD -> "
            "Role of region in protein function is unknown -> "
            "LoF variants in this exon are not frequent in the general population and "
            "exon is present in biologically-relevant transcript(s) -> "
            "Variant removes >10% of protein"
        ),
        PVS1PredictionSeqVarPath.NF6: (
            "Not predicted to undergo NMD -> "
            "Role of region in protein function is unknown -> "
            "LoF variants in this exon are not frequent in the general population and "
            "exon is present in biologically-relevant transcript(s) -> "
            "Variant removes <10% of protein"
        ),
        PVS1PredictionSeqVarPath.SS1: (
            "Exon skipping or use of a cryptic slice site disrupts reading frame and "
            "is predicted to undergo NMD -> "
            "Exon is present in biologically-relevant transcript(s)"
        ),
        PVS1PredictionSeqVarPath.SS2: (
            "Exon skipping or use of a cryptic slice site disrupts reading frame and "
            "is predicted to undergo NMD -> "
            "Exon is absent from biologically-relevant transcript(s)"
        ),
        PVS1PredictionSeqVarPath.SS3: (
            "Exon skipping or use of a cryptic slice site disrupts reading frame and "
            "is not predicted to undergo NMD -> "
            "Truncated/altered region is critical to protein function"
        ),
        PVS1PredictionSeqVarPath.SS4: (
            "Exon skipping or use of a cryptic slice site disrupts reading frame and "
            "is not predicted to undergo NMD -> "
            "Role of region in protein function is unknown -> "
            "LoF variants in this exon are frequent in the general population and/or "
            "exon is absent from biologically-relevant transcript(s)"
        ),
        PVS1PredictionSeqVarPath.SS5: (
            "Exon skipping or use of a cryptic slice site disrupts reading frame and "
            "is not predicted to undergo NMD -> "
            "Role of region in protein function is unknown -> "
            "LoF variants in this exon are not frequent in the general population and "
            "exon is present in biologically-relevant transcript(s) -> "
            "Variant removes >10% of protein"
        ),
        PVS1PredictionSeqVarPath.SS6: (
            "Exon skipping or use of a cryptic slice site disrupts reading frame and "
            "is not predicted to undergo NMD -> "
            "Role of region in protein function is unknown -> "
            "LoF variants in this exon are not frequent in the general population and "
            "exon is present in biologically-relevant transcript(s) -> "
            "Variant removes <10% of protein"
        ),
        PVS1PredictionSeqVarPath.SS7: (
            "Exon skipping or use of a cryptic slice site preserves reading frame -> "
            "Role of region in protein function is unknown -> "
            "LoF variants in this exon are frequent in the general population and/or "
            "exon is absent from biologically-relevant transcript(s)"
        ),
        PVS1PredictionSeqVarPath.SS8: (
            "Exon skipping or use of a cryptic slice site preserves reading frame -> "
            "Role of region in protein function is unknown -> "
            "LoF variants in this exon are not frequent in the general population and "
            "exon is present in biologically-relevant transcript(s) -> "
            "Variant removes >10% of protein"
        ),
        PVS1PredictionSeqVarPath.SS9: (
            "Exon skipping or use of a cryptic slice site preserves reading frame -> "
            "Role of region in protein function is unknown -> "
            "LoF variants in this exon are not frequent in the general population and "
            "exon is present in biologically-relevant transcript(s) -> "
            "Variant removes <10% of protein"
        ),
        PVS1PredictionSeqVarPath.SS10: (
            "Exon skipping or use of a cryptic slice site preserves reading frame -> "
            "Truncated/altered region is critical to protein function"
        ),
        PVS1PredictionSeqVarPath.IC1: (
            "No known alternative start codon in other transcripts -> "
            ">=1 pathogenic variant(s) upstream of closest potential in-frame start codon"
        ),
        PVS1PredictionSeqVarPath.IC2: (
            "No known alternative start codon in other transcripts -> "
            "No pathogenic variant(s) upstream of closest potential in-frame start codon"
        ),
        PVS1PredictionSeqVarPath.IC3: "Different functional transcript uses alternative start codon",
        PVS1PredictionStrucVarPath.NotSet: "Not Set",
        PVS1PredictionStrucVarPath.DEL1: "Full gene deletion",
        PVS1PredictionStrucVarPath.DEL2: (
            "Single to multi exon deletion disrupts reading frame and "
            "is predicted to undergo NMD -> "
            "Exon is present in biologically-relevant transcript(s)"
        ),
        PVS1PredictionStrucVarPath.DEL3: (
            "Single to multi exon deletion disrupts reading frame and "
            "is predicted to undergo NMD -> "
            "Exon is absent from biologically-relevant transcript(s)"
        ),
        PVS1PredictionStrucVarPath.DEL4: (
            "Single to multi exon deletion disrupts reading frame and "
            "is not predicted to undergo NMD -> "
            "Truncated/altered region is critical to protein function"
        ),
        PVS1PredictionStrucVarPath.DEL5_1: (
            "Single to multi exon deletion disrupts reading frame and "
            "is not predicted to undergo NMD -> "
            "Role of region in protein function is unknown -> "
            "LoF variants in this exon are frequent in the general population and/or "
            "exon is absent from biologically-relevant transcript(s)"
        ),
        PVS1PredictionStrucVarPath.DEL6_1: (
            "Single to multi exon deletion disrupts reading frame and "
            "is not predicted to undergo NMD -> "
            "Role of region in protein function is unknown -> "
            "LoF variants in this exon are not frequent in the general population and "
            "exon is present in biologically-relevant transcript(s) -> "
            "Variant removes >10% of protein"
        ),
        PVS1PredictionStrucVarPath.DEL7_1: (
            "Single to multi exon deletion disrupts reading frame and "
            "is not predicted to undergo NMD -> "
            "Role of region in protein function is unknown -> "
            "LoF variants in this exon are not frequent in the general population and "
            "exon is present in biologically-relevant transcript(s) -> "
            "Variant removes <10% of protein"
        ),
        PVS1PredictionStrucVarPath.DEL5_2: (
            "Single to multi exon deletion preserves reading frame -> "
            "Role of region in protein function is unknown -> "
            "LoF variants in this exon are frequent in the general population and/or "
            "exon is absent from biologically-relevant transcript(s)"
        ),
        PVS1PredictionStrucVarPath.DEL6_2: (
            "Single to multi exon deletion preserves reading frame -> "
            "Role of region in protein function is unknown -> "
            "LoF variants in this exon are not frequent in the general population and "
            "exon is present in biologically-relevant transcript(s) -> "
            "Variant removes >10% of protein"
        ),
        PVS1PredictionStrucVarPath.DEL7_2: (
            "Single to multi exon deletion preserves reading frame -> "
            "Role of region in protein function is unknown -> "
            "LoF variants in this exon are not frequent in the general population and "
            "exon is present in biologically-relevant transcript(s) -> "
            "Variant removes <10% of protein"
        ),
        PVS1PredictionStrucVarPath.DEL8: (
            "Single to multi exon deletion preserves reading frame -> "
            "Truncated/altered region is critical to protein function"
        ),
        PVS1PredictionStrucVarPath.DUP1: (
            "Proven in tandem -> " "Reading frame disrupted and NMD predicted to occur"
        ),
        PVS1PredictionStrucVarPath.DUP2_1: (
            "Proven in tandem -> " "No or unknown impact on reading frame and NMD"
        ),
        PVS1PredictionStrucVarPath.DUP2_2: (
            "Presumed in tandem -> " "No or unknown impact on reading frame and NMD"
        ),
        PVS1PredictionStrucVarPath.DUP3: (
            "Proven in tandem -> " "Reading frame presumed disrupted and NMD predicted to occur"
        ),
        PVS1PredictionStrucVarPath.DUP4: "Proven not in tandem",
    }



