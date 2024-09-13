.. _benchmarking:

============
Benchmarking
============

AutoACMG has undergone rigorous benchmarking against other established tools such as InterVar and
GeneBe to ensure accuracy and reliability. This section details the methodology and key findings
from our benchmarking exercises.

Methodology
-----------

We constructed a specialized dataset derived from the
`ClinGen Evidence Repository <https://clinicalgenome.org>`_. The dataset and the script used for the
benchmarking can be accessed from our GitHub repository:

- **Dataset**: `comparison_criteria_custom.csv <https://github.com/bihealth/auto-acmg/blob/main/src/bench/comparison_criteria_custom.csv>`_
- **Benchmarking Script**: Located in the `src/bench` directory of the repository.

The benchmarking process involved:

1. Comparing the predictions made by AutoACMG against those from InterVar and GeneBe using a custom
dataset.
2. Utilizing the APIs provided by InterVar and GeneBe to fetch predictions.
3. Computing statistical metrics such as kappa scores, F1 score, precision, and recall to evaluate
the performance.

How to Run the Benchmarking Script
----------------------------------

To run the benchmarking script, you'll need a development environment set up with the required
dependencies. Follow the guide in the `Development Environment Setup <dev_quickstart>`_ section to
set up the environment. Then run the following command:

.. code-block:: bash

    make bench

This command will execute the benchmarking script and generate a `stats.csv` file in the `src/bench`
directory. The file contains the statistical metrics computed during the benchmarking process. To
analyze the results, you can use the `results_analysis.ipynb` Jupyter notebook provided in the
`src/bench` directory. To run the notebook, just execute the following command:

.. code-block:: bash

    make jupyterlab

Results
-------

The results of the benchmarking are summarized in the `statistic_metrics.csv` file, which can be
viewed here:

- **Statistical Metrics**: `statistic_metrics.csv <https://github.com/bihealth/auto-acmg/blob/main/src/bench/statistic_metrics.csv>`_

Key Findings:

- **True Positives**: Number of instances where AutoACMG and the comparative tools agreed on a
    pathogenic variant.

- **False Positives**: Number of instances where AutoACMG predicted a variant as pathogenic, which
    was not confirmed by the other tools.

- **False Negatives**: Number of instances where AutoACMG did not predict a pathogenic variant, but
    the comparative tools did.

- **Kappa Score**: Measures the agreement between AutoACMG and the other tools beyond chance.

- **F1 Score, Precision, and Recall**: These metrics provide a more detailed insight into the
    accuracy of AutoACMG in identifying pathogenic variants compared to other tools.


Conclusion
----------

The benchmarking results indicate that AutoACMG performs comparably with, if not superior to, other
leading tools in the field. The precision and recall rates highlight AutoACMG's ability to reliably
identify pathogenic variants, making it a valuable tool for geneticists and researchers. Ongoing
improvements and updates will continue to refine its predictions and expand its utility in clinical
genomics.
