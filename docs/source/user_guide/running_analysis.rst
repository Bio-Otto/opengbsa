Running Analysis
================

Running OpenGBSA is straightforward using the command-line interface (CLI).

Basic Usage
-----------

To run an analysis, simply pass your configuration file to the CLI:

.. code-block:: bash

    opengbsa run config.yaml

Or using python module syntax:

.. code-block:: bash

    python -m mmgbsa.cli config.yaml

Parallel Processing
-------------------

OpenGBSA automatically detects available CPU cores.

- **Energy Calculation**: Performed sequentially or in batches (GPU acceleration supported via OpenMM).
- **Decomposition**: The per-residue decomposition step is **highly parallelized**. It will use all available cores by default.

.. note::
   To limit the number of threads (e.g., on a shared cluster), set the ``OPENMM_CPU_THREADS`` environment variable before running.

   .. code-block:: bash
   
       export OPENMM_CPU_THREADS=4
       opengbsa run config.yaml

GPU Acceleration
----------------

If you have a CUDA-compatible GPU and OpenMM configured with CUDA support, OpenGBSA will automatically utilize it for the GBSA energy calculations, significantly speeding up the process.

Logging and Monitoring
----------------------

Real-time progress is displayed in the terminal. A detailed log file (``mmgbsa.log``) is also saved in the output directory, recording every step, parameter, and any warnings.
