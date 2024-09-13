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

Running the Docker Image
------------------------

Once the Docker image is built, you can run it using the following command:

.. code-block:: bash

    docker run -d -p 8000:8000 --env-file .env.dev auto-acmg

Here, ``-d`` runs the container in detached mode (in the background), ``-p 8000:8000`` maps port
8000 of the container to port 8000 on the host, making the application accessible via
``localhost:8000``. The ``--env-file .env.dev`` option tells Docker to use the environment variables d
efined in the ``.env.dev`` file. Replace ``.env.dev`` with the path to your actual environment
file if different.

.. note::

    You have to configure the environment file before running the Docker container. See the next
    section for details.

Configuring the Environment File
--------------------------------

The application can be configured using environment variables. An example configuration file named
``.env.dev`` might look like this:

.. code-block:: none

    API_V1_STR=/api/v1
    DATABASE_URL=postgres://user:password@localhost/dbname
    SECRET_KEY=your_secret_key

Adjust the values according to your environment. Here are brief descriptions of the variables:

- ``API_V1_STR``: Base path for API endpoints.
- ``USE_CACHE``: Enable or disable caching of API responses.
- ``API_REEV_URL``: URL of the REEV API.
- ``SEQREPO_DATA_DIR``: Path to the SeqRepo data directory.

To pass this configuration to the Docker container, ensure the ``.env.dev`` file is located where
you run the ``docker run`` command or specify the correct path to the file using the ``--env-file``
option.

.. note::

    The ``.env.dev`` file has example values for demonstration purposes. You will probably need to
    adjust SEQREPO_DATA_DIR to point to the correct location on your system.

Accessing the OpenAPI Documentation
------------------------------------

Once the application is running, you can access the OpenAPI documentation by navigating to:

.. code-block:: none

    http://localhost:8000/api/v1/docs

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
at the URL: ``http://localhost:8000/api/v1/docs``.
