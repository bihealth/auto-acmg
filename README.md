# autopvs1
An automatic tool for the classification of PVS1 variants.

[![CI](https://github.com/bihealth/autopvs1/actions/workflows/main-ci.yml/badge.svg)](https://github.com/bihealth/autopvs1/actions/workflows/main-ci.yml)
[![Documentation Status](https://readthedocs.org/projects/autopvs1/badge/?version=latest)](https://autopvs1.readthedocs.io/en/latest/?badge=latest)

Inspirational paper: [Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6185798/)
Inspirational repo (AutoPVS1): https://github.com/JiguangPeng/autopvs1

Before the real documentation appears, please refer to the following information:

## Installation

For installaiton you need pyenv and pipeenv. If you don't have them, you have to install them first.

This project uses python 3.12, so add it to your pyenv:

```bash
pyenv install 3.12
```

Then, you can install the required packages with `pipenv` using the following command:

```bash
make deps
```

## Usage

Currently we provide CLI for the usage of the tool. You can use the following command run the prediction for a variant:

```bash
pipenv run python -m src.main <Variant-representation> --genome_release <Genome-release>
```

Alternatively, you can run the prediction for an example variant:

```bash
make example_run
```

For more information on the tool, you can use the following command:

```bash
pipenv run python -m src.main.py --help
```
