.. _dev_quickstart:

==============
Dev Quickstart
==============

This document describes how to get started with developing and using the project.

--------
Overview
--------

You will need to have Python 3.12 and pipenv installed.
The next step is to checkout the repository and install the Python dependencies.
Then, you will be able to utilize the CLI and run the tests.
The following assumes a Debian/Ubuntu machine; your mileage may vary.

-------------
Prerequisites
-------------

You can use `pyenv <https://github.com/pyenv/pyenv>`__ for getting a specific python version.

.. code-block:: bash

    sudo apt-get update; sudo apt-get install make build-essential libssl-dev zlib1g-dev \
    libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm \
    libncursesw5-dev xz-utils tk-dev libxml2-dev libxmlsec1-dev libffi-dev liblzma-dev

    curl https://pyenv.run | bash

Append the following to your ```~/.bashrc``:

.. code-block:: bash

    export PATH="$HOME/.pyenv/bin:$PATH"
    eval "$(pyenv init --path)"
    eval "$(pyenv virtualenv-init -)"


... and ensure to execute/source this as well (``exec $SHELL``).

Now you can install a specific python version:

.. code-block:: bash

    pyenv install 3.12
    pyenv local 3.12

Install pipenv:

.. code-block:: bash

    pip install --user pipenv

-----------------
Clone Rsepository
-----------------

.. code-block:: bash

    git clone git@github.com:bihealth/auto-acmg.git

Then you need to pull large files from the LFS:

.. code-block:: bash

    git lfs pull

--------------------
Install Dependencies
--------------------

You can use the provided ``Makefile`` files to install the dependencies.

.. code-block:: bash

    make deps

------------------
Set up the SeqRepo
------------------

AutoACMG uses the SeqRepo to fetch the reference sequences.
To set up the SeqRepo, you'll need to have the RefSeq reference genomes of versions GRCh37 and
GRCh38. You can get them as follows:

1. Install dependencies for SeqRepo:
SeqRepo uses tabix and bgzip to index and compress the reference genomes. You can install them
using the following command for Linux:

.. code-block:: bash

    sudo apt-get install tabix bgzip

and for MacOS:

.. code-block:: bash

    brew install htslib

2. Initialize the SeqRepo:
First, you'll need to setup the SeqRepo instance. You can do this by running the following command:

.. code-block:: bash

    sudo chown $USER /usr/local/share/seqrepo
    sudo mkdir -p /usr/local/share/seqrepo

Then you can initialize the SeqRepo instance:

.. code-block:: bash

    pipenv run seqrepo init -i auto-acmg

.. note::

    The ``-i auto-acmg`` flag is used to set the SeqRepo instance name to ``auto-acmg``.
    If you want to use a different default seqrepo directory, you can set the ``SEQREPO_DIR``
    environment variable or provide the ``-r`` flag to the ``seqrepo`` command.


2. Download the reference genomes:

The following command downloads the reference sequences. Note, that it'll take some time.

**Important**: There might be a bug in the SeqRepo that causes the download to fail (refer to the
`github issue <https://github.com/biocommons/biocommons.seqrepo/issues/166>`__). If this happens,
modify the source code by running the following command:

.. code-block:: bash

    sed -i -e 's/if aliases_cur.fetchone() is not None/if next(aliases_cur, None) is not None/' \
    <your-path-to-lib>/biocommons/seqrepo/cli.py

And the actual download command:

.. code-block:: bash

    pipenv run seqrepo fetch-load -i auto-acmg -n RefSeq NC_000001.10 NC_000002.11 NC_000003.11 NC_000004.11 NC_000005.9 NC_000006.11 NC_000007.13 NC_000008.10 NC_000009.11 NC_000010.10 NC_000011.9 NC_000012.11 NC_000013.10 NC_000014.8 NC_000015.9 NC_000016.9 NC_000017.10 NC_000018.9 NC_000019.9 NC_000020.10 NC_000021.8 NC_000022.10 NC_000023.10 NC_000024.9 NC_012920.1 NC_000001.11 NC_000002.12 NC_000003.12 NC_000004.12 NC_000005.10 NC_000006.12 NC_000007.14 NC_000008.11 NC_000009.12 NC_000010.11 NC_000011.10 NC_000012.12 NC_000013.11 NC_000014.9 NC_000015.10 NC_000016.10 NC_000017.11 NC_000018.10 NC_000019.10 NC_000020.11 NC_000021.9 NC_000022.11 NC_000023.11 NC_000024.10 NC_012920.1

.. note::

    The above RefSeq identifiers are for all of the chromosomes and the mitochondrial genome. Note,
    that this is a large (and slow) download and will take some time.


----------------------
Set up the `.env` file
----------------------

You need to create a `.env` file in the root of the project. The default settings
can be found in the `.env.dev` file. Copy the contents with the following command:

.. code-block:: bash

    cp .env.dev .env

**Important**: You need to set the `SEQREPO_DIR` variable in the `.env` file to the path of the
SeqRepo instance you created in the previous step. It should look similar to this:

.. code-block:: bash

    SEQREPO_DIR=/usr/local/share/seqrepo/auto-acmg

---------------
Running the CLI
---------------

After you have set up the SeqRepo and the `.env` file, you can serve the API.
You can use the following commands:

.. code-block:: bash

    make serve

Then go to the `http://0.0.0.0:8080/api/v1/docs` to see the API documentation.
