# Configuration Guide

OpenGBSA uses a configuration file in **YAML** format to manage the analysis. This file consists of 4 main sections.

Main template file: [`config_master.yaml`](../config_master.yaml)

---

## 1. Input (Input Files)

This section defines the building blocks of the system.

| Parameter | Mandatory? | Description |
|-----------|------------|-------------|
| `topology` | **YES** | System topology (`.top`, `.prmtop`). |
| `trajectory` | **YES** | Simulation trajectory (`.xtc`, `.dcd`, `.nc`). |
| `ligand_resname` | **YES** | 3-letter code of the ligand (E.g., `LIG`). If not specified, the program may fail to find the ligand. |
| `ligand_pdb` | No | Isolated structure file for ligand parameterization. |

> **GROMACS Tip:** If providing a `.top` file as `topology`, ensure that any `.itp` files `#include`d in it are present in the same directory.

---

## 2. Params (General Parameters)

Determines the general behavior of the program and forcefields.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `output_directory` | `results_analysis` | Name of the folder where results will be saved. |
| `protein_forcefield` | `amber` | Forcefield for the protein (`amber`, `charmm`). |
| `ligand_forcefield` | `gaff` | Forcefield for the ligand (`gaff`, `openff`). `openff` is modern and recommended. |

---

## 3. Analysis Settings

Contains the scientific details of the calculation.

### Basic Settings
*   **`gb_model`**: Implicit solvent model.
    *   `OBC2`: (Recommended) The most balanced model.
    *   `GBn2`: Alternative model.
*   **`salt_concentration`**: Salt concentration (Molar). Typically 0.15 M.

### Performance and Frame Selection
*   **`max_frames`**: Total number of frames to analyze.
    *   50-100 frames: Suitable for preliminary analysis.
    *   500+ frames: Recommended for publication-quality statistics.
*   **`frame_start` / `frame_stride`**: Determines which frames to skip.

### Advanced Analyses
*   **`run_per_residue_decomposition` (`true`/`false`)**:
    *   If set to `true`, the energy contribution of each amino acid (Heatmap) is calculated.
    *   **Note:** Increases processing time.
*   **`entropy_method`**:
    *   `interaction`: Interaction Entropy method (Fast, recommended).
    *   `quasiharmonic`: Calculated via trajectory covariance (Slow).
    *   `none`: Entropy is not calculated.

---

## 4. Visualization

Controls output formats.

*   `generate_report`: Generates HTML report.
*   `generate_pymol`: Generates PyMOL (.pml) file.
*   `pymol.energy_threshold`: Energy threshold for bonds shown in PyMOL (e.g., -1.0 kcal/mol).

---

## Surface Area Model (`sa_model`)

OpenGBSA supports two surface area calculation methods:

| Model | Description | Performance | Accuracy | Use Case |
|-------|-------------|-------------|----------|----------|
| **ACE** | Analytical Continuum Electrostatics (default) | Fast | Good | Initial screening, large datasets |
| **LCPO** | Linear Combinations of Pairwise Overlaps | Moderate | Better | Publication-quality, final binding affinity calculations |

### Configuration Example

```yaml
analysis_settings:
  gb_model: OBC2
  sa_model: LCPO  # or ACE (default)
  salt_concentration: 0.150
```

### Technical Notes

- **LCPO** matches AMBER MMPBSA.py default behavior
- Requires OpenMM 8.1+ with LCPO support (built from source)
- Automatically assigns atom-specific parameters via topology analysis
- Uses 1.4 Å probe radius (standard water molecule size)

---

