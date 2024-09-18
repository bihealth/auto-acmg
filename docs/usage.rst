.. _usage:

=====
Usage
=====

This section provides detailed instructions on how to set up and use the Docker container for the
AutoACMG system, configure environment variables, and access the API documentation.

Building the Docker Image
-------------------------

To build the Docker image, ensure Docker is installed on your system and navigate to the directory
containing the Dockerfile. Use the following command to build the image:

.. code-block:: bash

    docker build -t auto-acmg .

This command builds the Docker image with the tag ``auto-acmg``. You can replace ``auto-acmg`` with
any other tag suitable for your deployment or versioning schema.

Setting Up Volumes
------------------

The AutoACMG Docker container requires access to SeqRepo data directories for sequence information.
You must set up volumes that the Docker container can use to access this data without needing to
copy it into the container directly. Here's how you set up these volumes:

1. **SeqRepo Data Volume:** This volume stores the SeqRepo datasets.
2. **Custom Project Data Volume:** This volume stores the custom project seqrepo data.

Ensure these directories exist on your host and are populated with the necessary data:

.. code-block:: bash

    mkdir -p /usr/local/share/seqrepo
    chown -R root:root /usr/local/share/seqrepo

.. code-block:: bash

    pipenv run seqrepo init -i auto-acmg

.. code-block:: bash

    pipenv run seqrepo fetch-load -i auto-acmg -n RefSeq NC_000001.10 NC_000002.11 NC_000003.11 \
        NC_000004.11 NC_000005.9 NC_000006.11 NC_000007.13 NC_000008.10 NC_000009.11 NC_000010.10 \
        NC_000011.9 NC_000012.11 NC_000013.10 NC_000014.8 NC_000015.9 NC_000016.9 NC_000017.10 \
        NC_000018.9 NC_000019.9 NC_000020.10 NC_000021.8 NC_000022.10 NC_000023.10 NC_000024.9 \
        NC_012920.1 NC_000001.11 NC_000002.12 NC_000003.12 NC_000004.12 NC_000005.10 NC_000006.12 \
        NC_000007.14 NC_000008.11 NC_000009.12 NC_000010.11 NC_000011.10 NC_000012.12 NC_000013.11 \
        NC_000014.9 NC_000015.10 NC_000016.10 NC_000017.11 NC_000018.10 NC_000019.10 NC_000020.11 \
        NC_000021.9 NC_000022.11 NC_000023.11 NC_000024.10 NC_012920.1

.. note::

    The paths used in this example are for demonstration purposes! You should adjust the paths
    to match your actual directory structure and ensure that the directories have the correct
    permissions for Docker to access them. The ``/usr/local/share/seqrepo`` directory is the default
    location for Linux systems. The ``/home/auto-acmg/seqrepo/master`` directory is an example of
    a custom project data directory.


If the above doesn't work for you, you can try to download backups from the CUBI SharePoint. The
backups are located in the folder ``/Documents/Coding and Engineering/AutoACMG``. Then unarchive
them with the following command:

.. code-block:: bash

    tar -czvf seqrepo_local.tar.gz .dev/volumes/seqrepo/local --strip-components=1
    tar -czvf seqrepo_master.tar.gz .dev/volumes/seqrepo/master --strip-components=1

Finally, you should have the following directories structures in ``/usr/local/share/seqrepo`` (local):
and SEQREPO_DATA_DIR (master):


.. code-block:: bash

    seqrepo
    ├── master
    │   ├── aliases.sqlite3
    │   ├── sequences
    │          └── db.sqlite3
    │          └── 2024
    │                └── 1224
    │                └── ....
    │
    └── local
        ├── master
            ├── aliases.sqlite3
            ├── sequences
                └── db.sqlite3


Running the Docker Image
------------------------

Once the Docker image is built, you can run it using the following command to include the volumes:

.. code-block:: bash

    docker run -d -p 8080:8080 --env-file .env \
    -v /usr/local/share/seqrepo:/usr/local/share/seqrepo \
    -v /home/auto-acmg/seqrepo/master:/home/auto-acmg/seqrepo/master \
    auto-acmg

This command runs the container in detached mode (in the background), maps port 8080 of the
container to port 8080 on the host, making the application accessible via ``localhost:8080``.
The ``--env-file .env`` option tells Docker to use the environment variables defined in the
``.env`` file. Replace ``.env`` with the path to your actual environment file if different.
The ``-v`` flags map the local directories to their respective directories within the container,
ensuring that SeqRepo data is accessible.

.. note::

    You must configure the environment file before running the Docker container and ensure that
    the directories used for volumes are properly set up and have the correct permissions for
    Docker to access them. See the configuration details in the sections below.

Configuring the Environment File
--------------------------------

The application can be configured using environment variables. An example configuration file named
``.env.dev`` might look like this:

.. code-block:: none

    # Disable debug mode per default
    DEBUG=0
    # Disable cache to avoid memory issues
    USE_CACHE=0
    # Use the REEV API. Change it if you have other instance of REEV.
    API_REEV_URL=https://reev.cubi.bihealth.org/internal/proxy
    # Default path to seqrepo data for Docker. Change it to your local development path.
    # It can look like this: "/home/<username>/seqrepo/master"
    SEQREPO_DATA_DIR=/home/gromdimon/Custom/seqrepo/master

Adjust the values according to your environment. Here are brief descriptions of the variables. Note
that not all variables are required for the application to run. More info below.:

- ``DEBUG``: Enable or disable debug mode.
- ``USE_CACHE``: Enable or disable caching of API responses.
- ``CACHE_DIR``: Path to the cache directory.
- ``API_V1_STR``: Base path for API endpoints.
- ``API_REEV_URL``: URL of the REEV API.
- ``SEQREPO_DATA_DIR``: Path to the project-specific SeqRepo data directory.
- ``GENEBE_API_KEY``: API key for the GeneBE service. You'll need it for running the benchmarks.
- ``GENEBE_USERNAME``: Username for the GeneBE service. You'll need it for running the benchmarks.

You will most likely need to set the following variables:

- ``DEBUG``: Set to ``1`` to enable debug mode.
- ``USE_CACHE``: Set to ``1`` to enable caching. This is recommended only for development.
- ``SEQREPO_DATA_DIR``: Set to the path of the custom project SeqRepo data directory.

To pass this configuration to the Docker container, ensure the ``.env`` file is located where
you run the ``docker run`` command or specify the correct path to the file using the ``--env-file``
option.

.. note::

    Ensure that the environment variables are correctly set up and that the paths are valid and
    accessible by the Docker container.

Accessing the OpenAPI Documentation
------------------------------------

Once the application is running, you can access the OpenAPI documentation by navigating to:

.. code-block:: none

    http://localhost:8080/api/v1/docs

This URL provides interactive API documentation automatically generated from your OpenAPI specs. It
allows you to try out API calls directly from the browser.


API Endpoints
-------------

The API provides several endpoints for interacting with the AutoACMG system:

1. **Resolve Variant**
   Endpoint to resolve a variant based on its name and optionally specify the genome release.

   - **URL**: ``/api/v1/resolve``
   - **Method**: ``GET``
   - **Parameters**:
     - ``variant_name`` (required): The name or identifier of the variant.
     - ``genome_release`` (optional): The genome release version, defaults to ``GRCh38``.
   - **Success Response**: A JSON object containing resolved variant details.

   Example call:

   .. code-block:: none

       GET /api/v1/resolve?variant_name=chr1:228282272:G:A&genome_release=GRCh38

2. **Predict Sequence Variant**
   Endpoint to predict annotations for a sequence variant.

   - **URL**: ``/api/v1/predict/seqvar``
   - **Method**: ``GET``
   - **Parameters**:
     - ``variant_name`` (required): The name or identifier of the sequence variant.
   - **Success Response**: A JSON object containing prediction results.

   Example call:

   .. code-block:: none

       GET /api/v1/predict/seqvar?variant_name=chr1:228282272:G:A

3. **Predict Structural Variant**
   Endpoint to predict annotations for a structural variant.

   - **URL**: ``/api/v1/predict/strucvar``
   - **Method**: ``GET``
   - **Parameters**:
     - ``variant_name`` (required): The name or identifier of the structural variant.
     - ``duplication_tandem`` (optional): Specifies if the duplication is in tandem.
   - **Success Response**: A JSON object containing structural variant prediction results.

   Example call:

   .. code-block:: none

       GET /api/v1/predict/strucvar?variant_name=chr1:228282272:dup:Tandem


For more details on the API endpoints and their usage, refer to the OpenAPI documentation accessible
at the URL: ``http://localhost:8080/api/v1/docs``.
