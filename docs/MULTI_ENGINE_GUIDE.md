# Multi-Engine Support Guide

OpenGBSA is designed to be engine-agnostic by leveraging OpenMM's flexibility. However, each simulation engine (GROMACS, AMBER, NAMD, etc.) uses different file formats. This guide explains how to use OpenGBSA with outputs from various MD packages.

---

## 1. AMBER (Native Support)

AMBER is the native format for OpenGBSA. No conversion is required.

**Required Files:**
*   Topology: `.prmtop` (or `.parm7`)
*   Trajectory: `.nc` (NetCDF) or `.mdcrd` (ASCII, not recommended due to size)

**Configuration:**
```yaml
input:
  topology: complex.prmtop
  trajectory: prod.nc
  ligand_resname: LIG
```

---

## 2. GROMACS (Automatic Conversion)

OpenGBSA includes a built-in preprocessor for GROMACS files involving an automatic conversion to Amber format.

**Required Files:**
*   Topology: `.top` (and all included `.itp` files)
*   Trajectory: `.xtc` or `.trr`

**Configuration:**
```yaml
input:
  topology: topol.top
  trajectory: md_prod.xtc
  ligand_resname: LIG
```

> **Note:** See [GROMACS User Guide](GROMACS_GUIDE.md) for details on handling `.itp` dependencies.

---

## 3. OPENMM (Direct Support)

If you are already using OpenMM, you can provide the serialized `System` XML.

**Required Files:**
*   Topology: `system.xml` (Serialized OpenMM System) PLUS `system.pdb` (for atom names/residues).
    *   *Note:* The PDB file must have the same basename as the XML (e.g., `complex.xml` and `complex.pdb`) and be in the same folder.
*   Trajectory: `.dcd` or `.xtc`

**Configuration:**
```yaml
input:
  topology: complex.xml
  trajectory: trajectory.dcd
  ligand_resname: LIG
```

---

## 4. NAMD / CHARMM (Via Conversion)

OpenGBSA creates OpenMM Systems for analysis. While OpenMM supports CHARMM `.psf` files, accurate loading often requires multiple parameter files (`.prm`, `.str`) which can be complex to manage automatically.

**Recommended Strategy:**
We recommend converting your NAMD system to Amber format using **ParmEd** before analysis. This ensures all parameters are self-contained in a single file.

**Conversion Script (Python):**
```python
import parmed as pmd

# Load CHARMM system
# You must provide the PSF and all parameter files used in simulation
params = pmd.charmm.CharmmParameterSet('par_all36_prot.prm', 'toppar_water_ions.str', 'ligand.prm')
structure = pmd.load_file('step3_pbcsetup.psf', xyz='step3_pbcsetup.coor')
structure.load_parameters(params)

# Save as Amber
structure.save('complex.prmtop')
structure.save('complex.inpcrd')
```

**Configuration (after conversion):**
```yaml
input:
  topology: complex.prmtop
  trajectory: step5_production.dcd  # OpenMM can read DCD directly
  ligand_resname: LIG
```

---

## 5. DESMOND (Via Conversion)

Desmond outputs (`.cms`) are not natively supported by OpenMM. You must convert the system.

**Strategy:**
Use VMD or Schrödinger's tools to convert the structure to PDB, or use `intermol` / `parmed` to convert to GROMACS or AMBER format.

**Example (Files to generate):**
*   Convert `.cms` to `complex.pdb` (with correct connectivity).
*   Use standard OpenGBSA pipeline with `ligand_mol` and `ligand_pdb` to re-parameterize the system using OpenFF/Amber force fields.

```yaml
input:
  topology: complex.pdb
  trajectory: trajectory.dcd
  ligand_resname: LIG
params:
  ligand_forcefield: openff  # Re-assign parameters
```
