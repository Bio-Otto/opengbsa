Installation
============

Recommended: Conda Installation
------------------------------

OpenGBSA is best installed using conda to ensure all dependencies are handled correctly.

.. code-block:: bash

   conda create -n opengbsa python=3.10
   conda activate opengbsa
   conda install -c conda-forge opengbsa

.. tip::
   Always use a fresh conda environment for OpenGBSA to avoid dependency conflicts.

Troubleshooting & Compatibility
------------------------------

If you encounter issues, see the :doc:`troubleshooting` page for solutions to common problems (Python, NumPy, OpenMM, CUDA, etc.).

.. warning::
   Avoid pip installation unless you are an advanced user and understand the dependency requirements.

.. note::
   For CUDA support, ensure your system CUDA version matches the OpenMM requirements. 