# OpenGBSA: Advanced MM/GBSA Analysis Tool

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)
![OpenMM](https://img.shields.io/badge/OpenMM-8.0+-green.svg)

OpenGBSA is a comprehensive and automated tool for **Binding Free Energy**, **Energy Decomposition**, and **Entropy** calculations for protein-ligand complexes derived from molecular dynamics simulations.

It supports **GROMACS** (`.tpr`, `.top`, `.gro`), **Amber** (`.prmtop`), and **raw PDB + trajectory** inputs.

---

## 🚀 Features

*   **Multiple Input Formats:** GROMACS `.tpr`/`.top`, Amber `.prmtop`, and raw PDB files
*   **Multiple GB Models:** OBC1, OBC2, HCT, GBn, GBn2
*   **Multiple Surface Area Models:**
    *   **ACE (Analytical Continuum Electrostatics):** Fast default model
    *   **LCPO (Linear Combinations of Pairwise Overlaps):** More accurate, physics-based
*   **Per-Residue Energy Decomposition:** Parallel CPU multi-process decomposition
*   **Entropy Calculations:** Interaction entropy, quasiharmonic analysis, and normal mode analysis
*   **Dimer Mode:** Protein-protein-ligand binding analysis with separate subunit definitions
*   **GPU Acceleration:** CUDA and OpenCL platform support for fast energy evaluation
*   **Interactive HTML Report:** Full decomposition heatmaps, 3D structure viewer, and energy plots
*   **Flexible YAML Configuration:** Complete control via a single config file
*   **Comprehensive Test Suite:** 48 automated test configurations covering all features

---

## 🛠️ Installation

The recommended installation method is using `conda` (or `mamba`):

```bash
# Create a new environment
conda create -n opengbsa python=3.10
conda activate opengbsa

# Install core dependencies
conda install -c conda-forge openmm mdtraj rdkit parmed
pip install openff-toolkit matplotlib seaborn pandas pyyaml Jinja2

# Install the package (requires OpenMM 8.5.0beta from openmm_rc)
conda install -c bio-otto -c conda-forge/label/openmm_rc -c conda-forge opengbsa=0.0.6
```

---

## 📖 Quick Start

### 1. Prepare Configuration File

Copy the master config as a starting point:

```bash
cp config_master.yaml my_analysis.yaml
# Edit my_analysis.yaml with your input file paths and settings
```

### 2. Run Analysis

```bash
# Using the CLI entry point
mmgbsa my_analysis.yaml

# Or reference the CLI module directly
python -m mmgbsa.cli my_analysis.yaml
```

### 3. Minimal Config Example

```yaml
input_files:
  complex_pdb: path/to/complex.pdb       # or .tpr, .prmtop
  ligand_mol: path/to/ligand.sdf
  trajectory: path/to/traj.xtc

analysis_settings:
  temperature: 300.0
  gb_model: OBC2
  salt_concentration: 0.15
  max_frames: 100
  run_per_residue_decomposition: true
  decomp_frames: 20
  parallel_processing: true

platform_settings:
  preferred_platform: CUDA          # or OpenCL, CPU
  decomposition_platform: CPU       # CPU multiprocessing for decomposition

output_settings:
  output_directory: results/my_analysis
```

---

## ⚠️ Notes for GROMACS Users

*   If your `.top` file includes other `.itp` files, **all referenced `.itp` files MUST be present in the same folder.**
*   Prefer using `.tpr` files as input — the runner will automatically strip solvent and select protein + ligand atoms.
*   The trajectory (`.xtc`) atom count must match the topology. If you get an "Atom count mismatch", ensure waters/ions are handled consistently.

---

## 📊 Output Files

In your output directory you will find:

| File | Description |
|---|---|
| `interactive_report.html` | Full report: energy plots, 3D viewer, decomposition heatmap |
| `fixed_enhanced_mmgbsa_results_obc2.csv` | Per-frame binding energies |
| `per_residue_detailed.csv` | Per-residue energy contributions |
| `binding_hot_spots.csv` | Top residue hot spots |
| `energy_analysis.png` | Rolling average and convergence plot |
| `per_residue_decomposition.png` | Residue contribution heatmap |
| `final_report.txt` | Plain-text summary of results |

---

## 🧪 Running Tests

```bash
# Run the full comprehensive test suite (48 configs)
python test/run_comprehensive_tests.py

# Run a single named test
python test/run_comprehensive_tests.py test_gb_obc2

# All configs are in:
# test/configs/comprehensive/
```

---

## 🔧 Common Errors

| Error Message | Cause | Solution |
|---|---|---|
| `Invalid file path` | Config points to a non-existent file | Check all paths in your YAML config |
| `Atom count mismatch` | Topology and trajectory atom counts differ | Use a `.tpr` file or a pre-stripped PDB |
| `Ligand residues not found` | Ligand residue name not recognized | Set `ligand_resname: LIG` (or the correct name) in `input_files` |
| `decomp_frames > max_frames` | Config validation error | Ensure `decomp_frames` ≤ `max_frames` |
| `Invalid parameter value: temperature` | Integer instead of float | Use `300.0` not `300` |

---

## 🧬 Surface Area Models

| Model | Speed | Accuracy | Use Case |
|-------|-------|----------|----------|
| ACE   | Fast  | Good     | Rapid screening, large datasets |
| LCPO  | Moderate | Better | Publication-quality, final analysis |

---

## 📞 Support

See `config_master.yaml` for a fully annotated config with all available options.
See `docs/` for detailed guides:
- [Configuration Guide](docs/CONFIGURATION.md)
- [GROMACS Guide](docs/GROMACS_GUIDE.md)
- [Installation Guide](docs/INSTALLATION.md)
- [Output Files](docs/OUTPUTS.md)
