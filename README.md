# OpenGBSA

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)
![OpenMM 8.5](https://img.shields.io/badge/OpenMM-8.5-green.svg)
[![Documentation](https://img.shields.io/badge/docs-readthedocs-blue.svg)](https://opengbsa.readthedocs.io/)

OpenGBSA computes MM/GBSA binding enthalpies from molecular dynamics
trajectories. All energies are computed with OpenMM, on CPUs or on GPUs
through CUDA or OpenCL. It reads Amber, GROMACS, CHARMM/NAMD and OpenMM input
as well as plain PDB files, and it needs no AmberTools installation.

**Documentation:** <https://opengbsa.readthedocs.io/>

## Features

* **Two interfaces:** a command-line interface (CLI) for templates, checks and
  test data, and a configuration-driven interface (CDI): one YAML file defines
  a whole calculation.
* **Input:** Amber `.prmtop`, GROMACS `.top`/`.tpr`, CHARMM/NAMD `.psf`,
  OpenMM `.xml`, or a PDB file of the complex.
* **Solvation:** the five GB models of Amber (HCT, OBC1, OBC2, GBn, GBn2) with
  salt screening, and LCPO or ACE for the non-polar term.
* **Amber-compatible mode:** `amber_compatible: true` reproduces Amber
  MMPBSA.py on the same frames.
* **Per-residue decomposition** with plots, an interactive HTML report and a
  PyMOL script that colours the hot spots.
* **Entropy estimates:** interaction entropy and C2 with reliability limits,
  quasi-harmonic and normal-mode analysis.

## Installation

```bash
conda create -n opengbsa -c bio-otto -c conda-forge python=3.10 opengbsa
conda activate opengbsa
opengbsa --check-deps
```

Packages exist for Linux x86-64 and macOS on Apple Silicon. For GPUs, match
the CUDA runtime to the driver (`conda install -c conda-forge "cuda-version=12.9"`
for a driver that supports CUDA 12.9); see the
[installation guide](https://opengbsa.readthedocs.io/en/latest/installation.html).

## Quickstart

```bash
git clone https://github.com/Bio-Otto/opengbsa.git
cd opengbsa
opengbsa examples/quickstart/quickstart.yaml
```

This evaluates 50 frames of a beta-lactamase-inhibitor complex that ships with
the repository. Expected: `Mean binding energy (Enthalpy): -24.99 ± 4.60
kcal/mol`. For your own system, start from a template:

```bash
opengbsa --create-config             # short template: mmgbsa_config.yaml
opengbsa --create-complete-config    # every option with its default
opengbsa mmgbsa_config.yaml
```

A minimal configuration:

```yaml
input_files:
  complex_pdb: complex.prmtop     # or .top, .tpr, .psf, .xml, or a .pdb (Coordinate Mode)
  trajectory: production.nc
analysis_settings:
  gb_model: OBC2
  max_frames: 100
  run_per_residue_decomposition: true
```

Every option is described in the
[configuration reference](https://opengbsa.readthedocs.io/en/latest/user_guide/configuration.html).

## Testing and validation

```bash
conda install -c conda-forge pytest
python -m pytest test/unit                             # 116-117 passed, 13-14 skipped without test data
opengbsa --fetch-test-data --only core-fixtures        # 10 GB of test trajectories
python test/configs/check_validation.py --list         # third-party validation cases
```

Six third-party validation cases in `test/configs/` cover protein-ligand,
protein-peptide, protein-protein, protein-RNA, protein-DNA and antibody-antigen
complexes and compare with recorded results and with Amber MMPBSA.py. The
[testing guide](https://opengbsa.readthedocs.io/en/latest/testing.html) walks
through every step with the expected output.

## Citing OpenGBSA

Please cite OpenGBSA with the metadata in [`CITATION.cff`](CITATION.cff) (also
available through GitHub's "Cite this repository" button). Every release is
archived on Zenodo.

## Licence

MIT. The optional GROMACS `.tpr` reader TprParser is GPL-3.0; see
[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).
