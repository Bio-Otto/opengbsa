Troubleshooting
===============

This page lists common issues encountered with OpenGBSA and how to resolve them.

Installation & Environment Issues
---------------------------------

.. warning::
   **Problem:** Python 3.13 is installed, causing incompatibility with numpy or OpenMM
   
   **Error message:** ``No module named 'numpy.compat'`` or OpenMM/cudatoolkit errors
   
   **Cause:** Python 3.13 is not yet supported by numpy 1.x or OpenMM <8.3
   
   **Solution:**
     - Create a new conda environment with Python <3.13:
       .. code-block:: bash

          conda create -n opengbsa python=3.10
          conda activate opengbsa
          conda install -c conda-forge opengbsa

.. warning::
   **Problem:** numpy.compat error
   
   **Error message:** ``No module named 'numpy.compat'``
   
   **Cause:** numpy 2.x is installed, but some dependencies require numpy 1.x
   
   **Solution:**
     - Ensure numpy <2.0.0 is installed:
       .. code-block:: bash

          conda install numpy=1.26

.. warning::
   **Problem:** openmm/cudatoolkit version conflict
   
   **Error message:** openmm or cudatoolkit not found, or version mismatch
   
   **Cause:** openmm >=8.3 requires newer CUDA, which may not be available
   
   **Solution:**
     - Use openmm <8.3 (the conda package pins this automatically)

.. warning::
   **Problem:** CLI not found after install
   
   **Error message:** ``opengbsa: command not found``
   
   **Cause:** Environment not activated, or install failed
   
   **Solution:**
     - Activate the environment: ``conda activate opengbsa``
     - Reinstall if needed: ``conda install -c conda-forge opengbsa``

.. warning::
   **Problem:** YAML config errors
   
   **Error message:** YAML parsing error, missing keys, or invalid values
   
   **Cause:** Malformed YAML or missing required fields
   
   **Solution:**
     - Validate your YAML with an online linter
     - Use the example configs in the ``config/`` directory

.. warning::
   **Problem:** Missing input files
   
   **Error message:** FileNotFoundError or similar
   
   **Cause:** Input file path is incorrect or file is missing
   
   **Solution:**
     - Double-check all paths in your config file

General Tips
------------
.. tip::
   Always use a clean conda environment for OpenGBSA.

.. note::
   If you encounter persistent issues, try removing and recreating the environment.

.. important::
   Check the ``logs/`` directory for detailed error messages.

.. hint::
   For help, open an issue on GitHub: https://github.com/Bio-Otto/opengbsa/issues 