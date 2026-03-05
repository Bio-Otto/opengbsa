# GROMACS User Guide

This guide details the process of performing MM/GBSA analysis with OpenGBSA for researchers using GROMACS (`.top`, `.gro`, `.xtc`) files.

---

## 📂 Required Files

Before starting the analysis, ensure the following files are in your project folder:

1.  **Topology File (`.top`):**
    *   Example: `topol.top`
    *   **IMPORTANT:** If your topology file contains lines like `#include "ligand.itp"` or `#include "posre.itp"`, these files **MUST** be in the same folder.
2.  **Trajectory File (`.xtc` or `.trr`):**
    *   Example: `md_complex_prod.xtc`
    *   **Tip:** This trajectory file must match the atoms in the `.top` file exactly (water counts, ion counts, etc.). Using a `noPBC` trajectory with water removed via `trjconv` is typically recommended.
3.  **Ligand File (PDB/MOL):**
    *   Example: `ligand.pdb`
    *   Required for ligand parameterization (GAFF/OpenFF) when converting your system to Amber force field.

---

## ⚙️ Configuration

Configure GROMACS settings in your `config.yaml` file as follows:

```yaml
input:
  topology: topol.top
  trajectory: md_noPBC_noWater.xtc
  ligand_pdb: ligand.pdb
  ligand_resname: LIG  # The residue name (RESNAME) of the ligand in your topology

params:
  output_directory: gromacs_analysis_results
  protein_forcefield: amber
  ligand_forcefield: openff # or gaff
```

---

## 🔄 Automatic Conversion Process

OpenGBSA processes GROMACS files in the following steps:

1.  `topol.top` and referenced `.itp` files are read using Parmed.
2.  The GROMACS system is converted to Amber (`.prmtop`) format.
3.  Non-periodic (improper) torsions are fixed to be Amber-compatible during conversion (`fix_torsions`).
4.  MM/GBSA analysis (GB Model: OBC2) is run on the generated Amber topology.

---

## 🐛 Common Issues

### 1. "Included file not found" Error
**Cause:** Your `.top` file calls an `.itp` file, but that file is not in the folder.
**Solution:** Copy the missing `.itp` file to the analysis folder.

### 2. "Atom count mismatch" / "Coordinate mismatch"
**Cause:** The number of atoms in the `.top` file does not match the `.xtc` file.
**Solution:** This usually happens if your `.top` file has water (SOL) but you removed water from the trajectory.
*   Either use an `.xtc` that contains water (but this can be huge).
*   Or create a copy of your `.top` file, remove the SOL (water) lines, and use that copy.

### 3. Ligand Not Recognized
**Cause:** The `ligand_resname` parameter is entered incorrectly.
**Solution:** Open your `.pdb` or `.gro` file, check the 3-letter name of the ligand (e.g., UNL, LIG, DRG), and correct the config file.
