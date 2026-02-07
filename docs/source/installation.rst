Installation
============

Quick Installation
------------------

OpenGBSA is distributed via Conda.

.. code-block:: bash

    # Create a new environment
    conda create -n mmgbsa python=3.10 -y
    conda activate mmgbsa

    # Install OpenGBSA and dependencies
    conda install -c conda-forge -c omnia mdtraj openmm openff-toolkit pymol-open-source
    pip install opengbsa

Requirements
------------

- Linux or macOS
- Python 3.9+
- **OpenMM** (with CUDA support recommended for performance)
- **MDTraj**
- **OpenForceField** (for ligand parameterization)
- **ParmEd**