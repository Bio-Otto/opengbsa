# Installation Guide

Follow these steps to install OpenGBSA. The recommended method is using **Conda** (or Mamba).

---

## 1. Prerequisites

*   **Operating System:** Linux (Recommended: Ubuntu 20.04+) or macOS. WSL2 is recommended for Windows.
*   **Package Manager:** [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Mambaforge](https://github.com/conda-forge/miniforge).

---

## 2. Creating the Conda Environment

OpenGBSA relies on scientific libraries like `openmm`, `mdtraj`, and `openff-toolkit`. Create a clean environment to manage these dependencies:

```bash
# 1. Create a new environment (Python 3.10 is recommended)
conda create -n mmgbsa python=3.10

# 2. Activate the environment
conda activate mmgbsa
```

---

## 3. Installing Dependencies

Install the required scientific packages from the `conda-forge` channel:

```bash
# Core molecular dynamics and chemistry libraries
conda install -c conda-forge openmm mdtraj openff-toolkit rdkit parmed

# Helper libraries (Visualization and Data Processing)
pip install matplotlib seaborn pandas pyyaml Jinja2
```

> **Note:** If you plan to use GPU acceleration, ensure you install the OpenMM package version compatible with your system's CUDA version (Conda typically handles this automatically).

---

## 4. Verifying Installation

To verify that the installation was successful, check the versions:

```bash
python -c "import openmm; print('OpenMM:', openmm.__version__)"
python -c "import mdtraj; print('MDTraj:', mdtraj.__version__)"
```

---

## 5. Running the Project

Navigate to the directory where you downloaded the project files:

```bash
cd /path/to/mmgbsa_project/
```

Perform a test run:

```bash
python run_mmpbsa.py --help
```

If you see the help message, the installation is complete.
