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

    $ sudo apt-get update; sudo apt-get install make build-essential libssl-dev zlib1g-dev \
    libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm \
    libncursesw5-dev xz-utils tk-dev libxml2-dev libxmlsec1-dev libffi-dev liblzma-dev
    $ curl https://pyenv.run | bash

Append the following to your ```~/.bashrc``:

.. code-block:: bash

    $ pyenv
    export PATH="$HOME/.pyenv/bin:$PATH"
    eval "$(pyenv init --path)"
    eval "$(pyenv virtualenv-init -)"


... and ensure to execute/source this as well (``exec $SHELL``).

Now you can install a specific python version:

.. code-block:: bash

    $ pyenv install 3.12
    $ pyenv local 3.12

Install pipenv:

.. code-block:: bash

    $ pip install --user pipenv

-----------------
Clone Rsepository
-----------------

.. code-block:: bash

    $ git clone git@github.com:bihealth/auto-acmg.git
  
--------------------
Install Dependencies
--------------------

You can use the provided ``Makefile`` files to install the dependencies.

.. code-block:: bash

    $ make deps

----------------------
Set up the `.env` file
----------------------

You need to create a `.env` file in the root of the project. The default settings
can be found in the `.env.dev` file. Copy the contents with the following command:

.. code-block:: bash

    $ cp .env.dev .env

---------------
Running the CLI
---------------

You can run the CLI with the following command:

.. code-block:: bash

    $ make run VAR="NM_000152.4:c.1A>G" GR="GRCh37"

Also there's example for usage of CLI:

.. code-block:: bash

    $ make example_run