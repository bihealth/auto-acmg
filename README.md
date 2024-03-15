# autopvs1
An automatic tool for the classification of PVS1 variants.

[![CI](https://github.com/bihealth/autopvs1/actions/workflows/main-ci.yml/badge.svg)](https://github.com/bihealth/autopvs1/actions/workflows/main-ci.yml)

Inspirational repo (AutoPVS1): https://github.com/JiguangPeng/autopvs1

Before the real documentation appears, please refer to the following information:

## Installation

For installaiton you need pyenv and pipeenv. If you don't have them, please install them first.

This project uses python 3.12, so you need to install it first.

```bash
pyenv install 3.12
```

Then, you can install the required packages using the following command:

```bash
make deps
```

## Usage

Currently we provide CLI for the usage of the tool. You can use the following command to see the help message:

```bash
pipenv run python -m src.main 13-113803407-G-A --genome_release hg19
```

For more information, you can use the following command:

```bash
pipenv run python -m src.main.py --help
```
