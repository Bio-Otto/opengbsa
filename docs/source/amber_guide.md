# AMBER User Guide

This guide details the process of performing MM/GBSA analysis with OpenGBSA for researchers using **AMBER** (`.prmtop`, `.inpcrd`, `.nc`) files.

---

## 📂 Required Files

For AMBER systems, OpenGBSA works natively without conversion. Ensure you have:

1.  **Topology File (`.prmtop`):**
    *   Contains molecular parameters and topology.
    *   Example: `complex.prmtop`
2.  **Trajectory File (`.nc`, `.mdcrd`, `.dcd`):**
    *   Contains the simulation coordinates.
    *   Example: `production.nc`
3.  **Ligand Residue Name:**
    *   Know the 3-letter code of your ligand (e.g., `LIG`).

---

## ⚙️ Configuration

Configure the `input` section of your `config.yaml` to point to your AMBER files:

```yaml
input:
  topology: complex.prmtop
  trajectory: production.nc
  ligand_resname: LIG

analysis_settings:
  gb_model: OBC2 # Recommended for AMBER
  sa_model: LCPO # Recommended (matches MMPBSA.py)
  salt_concentration: 0.150
  protein_forcefield: amber
  ligand_forcefield: openff # Not used for native AMBER but required by schema
```

---

## 🚀 Running Analysis

```bash
mmgbsa config.yaml
```

---

## 💡 Tips for AMBER Users

### 1. Stripping Water/Ions
While OpenGBSA can handle solvated systems by masking, it is **highly recommended** to use a "clean" (dry) trajectory for faster processing.
- You can strip water/ions using `cpptraj`:
  ```bash
  cpptraj -p solvated.prmtop -y solvated.nc -x complex.nc -s complex.prmtop << EOF
  strip :WAT,Na+,Cl-
  autoimage
  EOF
  ```
- Then use `complex.prmtop` and `complex.nc` in your config.

### 2. Consistency
Ensure your `.prmtop` matches the atoms in your trajectory. If you use a solvated topology with a stripped trajectory, the analysis will fail with an atom count mismatch.

### 3. Periodic Boundary Conditions (PBC)
OpenGBSA handles basic PBC imaging if the trajectory information is correct, but pre-imaging with `cpptraj` (autoimage) is safer.
