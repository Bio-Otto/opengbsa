# OpenGBSA: Advanced MM/GBSA Analysis Tool

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)
![OpenMM](https://img.shields.io/badge/OpenMM-8.0+-green.svg)

OpenGBSA is a comprehensive and automated tool for **Binding Free Energy**, **Energy Decomposition**, and **Entropy** calculations for protein-ligand complexes derived from molecular dynamics simulations.

It is particularly optimized for **GROMACS** and **Amber** users.

---

## 🚀 Features

*   **Automatic Conversion:** Automatically converts GROMACS `.top` and `.xtc` files to Amber format and prepares them for MM/GBSA analysis.
*   **Multiple Surface Area Models:**
    *   **ACE (Analytical Continuum Electrostatics):** Fast, default model
    *   **LCPO (Linear Combinations of Pairwise Overlaps):** More accurate, physics-based surface area calculation
*   **Advanced Statistical Visualization:**
    *   **Rolling Average:** Shows whether the system has reached equilibrium.
    *   **Convergence Plot:** Analyzes the convergence status of binding energy.
    *   **Component Pie Chart:** Displays contribution ratios of VDW, Electrostatic, and Solvation energies.
    *   **Entropy Convergence:** Monitors the change of the entropy term over time.
*   **Detailed HTML Report:** Generates a professional report containing all charts and result tables.
*   **Per-Residue Decomposition:** Shows which amino acid contributes how much to binding (via Heatmap).
*   **Flexible Configuration:** The entire process is controlled by a single YAML file.

---

## 🛠️ Installation

The recommended installation method is using `conda` (or `mamba`):

```bash
# Create a new environment
conda create -n mmgbsa python=3.10
conda activate mmgbsa

# Install basic dependencies
conda install -c conda-forge openmm mdtraj openff-toolkit rdkit parmed
pip install matplotlib seaborn pandas pyyaml Jinja2
```

---

## 📖 Quick Start

### 1. Preparing Configuration File
To create a template containing all settings:

```bash
# Creates a sample config file
# (This file is also available as 'config_master.yaml' in the project root)
```

You can copy `config_master.yaml` to `my_config.yaml` and edit it.

### 2. Starting Analysis

```bash
# Run analysis
python run_mmpbsa.py my_config.yaml
```

---

## ⚠️ Critical Information for GROMACS Users

Things to consider when working with GROMACS files (`.top`, `.xtc`):

### 1. File Dependencies (.itp)
If your `topol.top` file references other files (e.g., `#include "ligand.itp"`), **all these files MUST be present in the analysis directory.**
If the program cannot find these references while reading the .top file, it gives a **"File not found"** error and analysis stops.

**Example Folder Structure:**
```text
my_project/
├── topol.top         <-- Main topology
├── ligand.itp        <-- Included in topology
├── protein.itp       <-- Included in topology
├── md_prod.xtc       <-- Trajectory
├── ligand.pdb        <-- Ligand structure (For parameterization)
└── config.yaml       <-- Configuration
```

### 2. Atom Count Mismatch
The number of atoms in the `trajectory` (xtc) and `topology` (top) file in your configuration file must match exactly.
*   If you get an "Atom count mismatch" error during analysis, ensure your `.xtc` file represents the same system as your `.top` file (Are waters removed? Are ions present?).

---

## 📊 Output Files

When analysis is complete, the following are created in the `results` folder:

1.  **Reports:**
    *   `advanced_analysis_report_LIG.html`: Interactive report containing all charts and summary.
    *   `fixed_enhanced_mmgbsa_results_obc2.csv`: Detailed energy values for each frame.

2.  **Plots (in `plots/` folder):**
    *   `rolling_average_LIG.png`: Energy stability chart.
    *   `convergence_plot_LIG.png`: Cumulative average chart.
    *   `energy_heatmap_LIG.png`: Amino acid contribution heatmap.
    *   `components_pie_LIG.png`: Energy components pie chart.

3.  **Visualization:**
    *   `view_binding.pml`: Ready-to-use session file for PyMOL (Shows important interactions in 3D).

---

## 🔧 Common Errors

| Error Message | Cause | Solution |
|-------------|-------|-------|
| `File not found in topology include` | .itp files in .top file are missing | Copy all .itp files to the working directory. |
| `Atom count mismatch` | Topology and Trajectory atom counts differ | Use an .xtc suitable for your .top file (e.g., water-free). |
| `Ligand Residue not found` | Ligand name is incorrect in config file | Check `ligand_resname` parameter in `config.yaml` (default: LIG). |

---

## 📞 Support

If you encounter issues, please review the detailed comments in the `config_master.yaml` file.
---

## 🧬 LCPO Surface Area Model

OpenGBSA now supports the **LCPO (Linear Combinations of Pairwise Overlaps)** method for more accurate surface area calculations.

### Why LCPO?

- **Physics-based:** More accurate than the analytical ACE approximation
- **AMBER-compatible:** Matches MMPBSA.py default behavior
- **Validated:** Tested against reference AMBER implementations

### Configuration Example

```yaml
analysis_settings:
  gb_model: OBC2
  sa_model: LCPO  # or ACE (default)
  salt_concentration: 0.150
```

### Performance Comparison

| Model | Computation Time | Accuracy | Use Case |
|-------|-----------------|----------|----------|
| ACE   | Fast            | Good     | Rapid screening, large datasets |
| LCPO  | Moderate        | Better   | Publication-quality, final analysis |

**Note:** LCPO requires OpenMM 8.1+ (built from source with LCPO support).

