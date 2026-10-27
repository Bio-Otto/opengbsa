Installation
============

OpenGBSA is distributed as a conda package for Linux (x86-64) and macOS (Apple
Silicon). The package brings OpenMM and every other dependency, so no
AmberTools installation is needed.

1. Requirements
---------------

* **Operating system:** Linux x86-64 or macOS on Apple Silicon. Windows is not
  supported natively.
* **Conda:** `Miniforge <https://github.com/conda-forge/miniforge>`_ is
  recommended. Miniconda and Mamba work as well.
* **GPU (optional):** an NVIDIA GPU with a recent driver for the CUDA platform,
  or any OpenCL device. Without a GPU, OpenGBSA runs on the CPU.

2. Install the package
----------------------

Create an environment with OpenGBSA from the ``bio-otto`` and ``conda-forge``
channels. The command below installs the version this documentation describes:

.. parsed-literal::

    conda create -n opengbsa -c bio-otto -c conda-forge python=3.10 opengbsa=\ |release|
    conda activate opengbsa

Leave out ``=``\ |release| to install the latest release.

3. Check the installation
-------------------------

Three commands check that OpenGBSA and OpenMM work:

.. code-block:: bash

    opengbsa --version
    opengbsa --check-deps
    python -m openmm.testInstallation

``opengbsa --version`` prints the version. ``opengbsa --check-deps`` imports
every dependency and ends with ``All dependencies are available!``.
``python -m openmm.testInstallation`` lists the OpenMM platforms that work on
this machine, for example ``Reference``, ``CPU``, ``CUDA`` and ``OpenCL``, and
checks that they compute the same forces.

4. GPU setup
------------

OpenGBSA uses the fastest platform that works: CUDA, then OpenCL, then CPU. A
requested platform that cannot be used is replaced by the next one, and the log
says so in a line that starts with ``PLATFORM FALLBACK``. Set
``platform_settings.strict_platform: true`` to make this an error instead.

**CUDA runtime and driver.** The CUDA platform needs a CUDA runtime that the
NVIDIA driver supports. Conda may install a newer runtime than the driver can
load. ``nvidia-smi`` shows the highest CUDA version the driver supports, top
right. Pin the runtime to that version, for example 12.9:

.. code-block:: bash

    nvidia-smi
    conda install -n opengbsa -c conda-forge "cuda-version=12.9"
    python -m openmm.testInstallation      # CUDA should now pass

**OpenCL** needs no CUDA runtime and also runs on non-NVIDIA GPUs.
**Several GPUs:** OpenGBSA uses the first visible device. Choose another with
``CUDA_VISIBLE_DEVICES``, for example ``CUDA_VISIBLE_DEVICES=1 opengbsa config.yaml``.

5. Optional components
----------------------

* **GROMACS .tpr input** uses `TprParser <https://pypi.org/project/TprParser/>`_
  (GPL-3.0). The conda package installs it; a source installation needs
  ``pip install TprParser``. A GROMACS ``.top`` works without it.
* **PyMOL** opens the hot-spot script that the per-residue decomposition
  writes: ``conda install -c conda-forge pymol-open-source``.
* **pytest** runs the unit tests (see :doc:`testing`):
  ``conda install -c conda-forge pytest``.

6. Install from source
----------------------

To run the tests or change the code, clone the repository and install it in
editable mode into an environment with the same dependencies as the package:

.. parsed-literal::

    git clone https://github.com/Bio-Otto/opengbsa.git
    cd opengbsa
    git checkout v\ |release|
    conda create -n opengbsa-dev -c conda-forge python=3.10 "openmm>=8.5.1,<8.6" \\
        mdtraj openff-toolkit openmmforcefields prolif rdkit mdanalysis parmed \\
        "numpy<2" scipy pandas matplotlib seaborn pyyaml plotly pathlib2 pooch \\
        "setuptools<81" pytest
    conda activate opengbsa-dev
    pip install --no-deps -e .
    opengbsa --version

7. Uninstall
------------

.. code-block:: bash

    conda env remove -n opengbsa
