# Installation Guide

Follow these steps to install OpenGBSA. The recommended method is using **Conda** (or Mamba).

---

## 1. Prerequisites

*   **Operating System:** Linux (x86-64) or macOS (Apple Silicon). Windows is not supported natively; WSL2 may work but is untested.
*   **Package Manager:** [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Mambaforge](https://github.com/conda-forge/miniforge).

---

## 2. Installing OpenGBSA

One command creates an environment with OpenGBSA and all its dependencies,
including OpenMM, from the `bio-otto` and `conda-forge` channels:

```bash
conda create -n opengbsa -c bio-otto -c conda-forge python=3.10 opengbsa
conda activate opengbsa
```

Packages are built and tested for Linux (x86-64) and macOS (Apple Silicon).
No AmberTools, sander or APBS installation is needed.

---

## 3. GPU (CUDA) setup

OpenMM's CUDA platform needs a CUDA runtime that your NVIDIA driver supports.
Conda may install a newer runtime than the driver can load; OpenMM then cannot
use the GPU and OpenGBSA falls back to the CPU with a `PLATFORM FALLBACK`
warning. Match the runtime to the driver: `nvidia-smi` shows the highest CUDA
version the driver supports (top right), for example 12.9:

```bash
conda install -c conda-forge "cuda-version=12.9"
```

To make a fallback an error instead of a warning, set
`platform_settings.strict_platform: true` in the configuration. The OpenCL
platform needs no CUDA runtime.

### Installing from source

For development, install the dependencies from `conda-forge` and the package in
editable mode:

```bash
conda create -n opengbsa-dev -c conda-forge python=3.10 "openmm>=8.5.1,<8.6" mdtraj \
    openff-toolkit openmmforcefields rdkit parmed prolif mdanalysis scipy pandas \
    matplotlib seaborn pyyaml plotly pytest
conda activate opengbsa-dev
pip install --no-deps -e .
```

---

## 4. Verifying Installation

To verify that the installation was successful, check the versions:

```bash
python -c "import openmm; print('OpenMM:', openmm.__version__)"
python -c "import mdtraj; print('MDTraj:', mdtraj.__version__)"
```

---

## 5. Running the Project

Perform a test run using the installed `opengbsa` command:

```bash
opengbsa --help
```

If you see the help message, the installation is complete.
