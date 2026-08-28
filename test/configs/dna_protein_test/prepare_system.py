#!/usr/bin/env python3
"""
Builds the OpenGBSA-ready complex PDB + trajectory for the IRF11-MDA5
promoter DNA complex (Zenodo 14377950) from the raw multi-model PDB
trajectory fetched by fetch_data.sh, and (optionally) an independent Amber
reference topology for cross-checking OpenGBSA's own result.

Unlike the GROMACS-sourced PPI datasets (1GCQ/2OOB), this dataset ships no
separate coordinate/trajectory pair -- the single downloaded file IS the
full 101-frame trajectory (a multi-model PDB). It also needed no bond-
inference workaround (this file has no CONECT records and neither mdtraj
nor OpenMM produced spurious bonds loading it, unlike the GROMACS .gro
conversions). The one fix still needed: like every Amber-family force
field, amber14-all.xml has no generic "HIS" residue template -- the raw
file's 4 histidines are relabeled HID/HIE/HIP based on which of HD1/HE2
they actually contain (all 4 turn out to be HIE here).

See ../../../../OXA-MD/Publishing/dna_protein_validation/README.md for the
full validation write-up this script's output was used to produce.

Usage:
    python3 prepare_system.py
    python3 prepare_system.py --build-amber-reference
"""
import argparse
import sys
from pathlib import Path

import mdtraj as md

HERE = Path(__file__).resolve().parent
RAW_PDB = HERE / "raw" / "IRF11-MDA5promoter_complex-traj.pdb"


def fix_his_protonation(pdb_path):
    lines = open(pdb_path).readlines()
    his_atoms = {}
    for line in lines:
        if line.startswith("ATOM") and line[17:20].strip() in ("HIS", "HID", "HIE", "HIP"):
            resnum = line[22:26].strip()
            atom_name = line[12:16].strip()
            his_atoms.setdefault(resnum, set()).add(atom_name)
    resnum_to_form = {}
    for resnum, atoms in his_atoms.items():
        has_hd1, has_he2 = "HD1" in atoms, "HE2" in atoms
        resnum_to_form[resnum] = "HIP" if (has_hd1 and has_he2) else "HID" if has_hd1 else "HIE"
    out = []
    for line in lines:
        if line.startswith("ATOM") and line[17:20].strip() in ("HIS", "HID", "HIE", "HIP"):
            resnum = line[22:26].strip()
            line = line[:17] + resnum_to_form[resnum] + line[20:]
        out.append(line)
    open(pdb_path, "w").writelines(out)
    return resnum_to_form


def prepare(out_dir):
    if not RAW_PDB.exists():
        sys.exit(f"Missing raw input {RAW_PDB} -- run ./fetch_data.sh first.")

    out_dir.mkdir(parents=True, exist_ok=True)
    ref_out = out_dir / "complex_ref.pdb"
    dcd_out = out_dir / "complex_traj.dcd"

    t = md.load(str(RAW_PDB))
    print(f"Loaded {t.n_frames} frames, {t.n_atoms} atoms, "
          f"{t.topology.n_chains} chains from raw trajectory.")

    t[0].save_pdb(str(ref_out))
    his_forms = fix_his_protonation(ref_out)
    print(f"HIS residues relabeled: {his_forms or '(none present)'}")

    t.save_dcd(str(dcd_out))
    print(f"Wrote {ref_out.name} and {dcd_out.name} ({t.n_frames} frames).")
    return ref_out, dcd_out


def build_amber_reference_openmm(ref_pdb, out_dir):
    """OpenMM ForceField.createSystem + ParmEd openmm.load_topology,
    exported directly to .prmtop -- same route used for 2OOB, chosen here
    since no atom-reorder was needed (see README)."""
    import parmed as pmd
    from openmm import app, unit

    topology = app.PDBFile(str(ref_pdb)).topology
    positions = app.PDBFile(str(ref_pdb)).positions
    ff = app.ForceField("amber14-all.xml")
    system = ff.createSystem(topology, nonbondedMethod=app.NoCutoff,
                              constraints=None, rigidWater=False)
    struct = pmd.openmm.load_topology(topology, system, xyz=positions.value_in_unit(unit.angstrom))
    pmd.tools.changeRadii(struct, "mbondi2").execute()
    struct.save(str(out_dir / "complex.prmtop"), overwrite=True)
    struct.save(str(out_dir / "complex.inpcrd"), overwrite=True)
    print(f"Amber reference (OpenMM+ParmEd): {out_dir / 'complex.prmtop'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-amber-reference", action="store_true",
                     help="Also build an independent Amber reference topology "
                          "(OpenMM+ParmEd) for cross-checking.")
    args = ap.parse_args()

    out_dir = HERE / "prepared"
    ref_pdb, _ = prepare(out_dir)

    if args.build_amber_reference:
        ref_dir = out_dir / "amber_reference"
        ref_dir.mkdir(parents=True, exist_ok=True)
        build_amber_reference_openmm(ref_pdb, ref_dir)

    print("\nDone. Config: irf11_dna_config.yaml (already points at prepared/).")


if __name__ == "__main__":
    main()
