#!/usr/bin/env python3
"""
Builds the OpenGBSA-ready complex PDB + trajectory for a Zenodo-6638504
protein-protein complex (1GCQ or 2OOB) from the raw GROMACS .gro/.xtc files
fetched by fetch_data.sh, and (optionally) an independent Amber reference
topology for cross-checking OpenGBSA's own result.

Both systems ship only coordinates (.gro) and a trajectory (.xtc) -- no
topology (.top/.itp/.tpr) -- so every downstream tool (OpenGBSA's own
OpenMM-based Coordinate Mode, and any independent Amber reference) must
parameterize the structure from scratch via a residue-template force field
(here, Amber ff14SB / amber14-all.xml). Two complications apply uniformly to
both consumers and are handled once, here, rather than separately per tool:

1. Naively converting .gro -> .pdb with mdtraj (or loading directly with
   OpenMM's app.PDBFile) lets bond inference run on interatomic distance,
   which can spuriously bond atoms that are merely close together in space
   (confirmed for 1GCQ: GLN11-O to the next residue's amide H). Writing the
   PDB with NO CONECT/bond records and letting the residue-template matcher
   infer standard amino-acid connectivity avoids this entirely.
2. GROMACS' generic "HIS" residue name doesn't encode protonation state.
   Amber's ff14SB has no generic HIS template -- it needs HID/HIE/HIP
   depending on which of HD1/HE2 are actually present in the structure.

See ../../../../OXA-MD/Publishing/ppi_validation/README.md (this repo's
sibling analysis directory) for the full validation write-up this script's
output was used to produce.

Usage:
    python3 prepare_system.py --system 1gcq
    python3 prepare_system.py --system 2oob
    python3 prepare_system.py --system 1gcq --build-amber-reference
"""
import argparse
import subprocess
import sys
from pathlib import Path

import mdtraj as md
import numpy as np

HERE = Path(__file__).resolve().parent

SYSTEMS = {
    "1gcq": {
        "raw_dir": "raw/1GCQ/complex/without_water",
        "n_frames_total": 7501,
    },
    "2oob": {
        "raw_dir": "raw/2OOB/complex/without_water",
        "n_frames_total": 7501,
    },
}

N_LAST_FRAMES = 300


def build_nobonds_pdb(gro_path, out_pdb):
    """Writes a PDB with no CONECT/bond records, so downstream tools infer
    bonds from residue templates instead of interatomic distance.

    The source .gro has no chain markers of its own (mdtraj loads it as a
    single chain) -- the two protein chains are distinguished only by
    residue numbering resetting back to 1 partway through (confirmed: 1GCQ's
    129 residues are numbered 1-59, then 1-70). That reset is used here to
    split into two chains, since OpenGBSA's `chainid N` selection strings
    (used by this system's config) need real chain boundaries to select on.
    """
    t = md.load(str(gro_path))
    top = md.Topology()
    chain = top.add_chain()
    prev_resSeq = None
    for r in t.topology.residues:
        if prev_resSeq is not None and r.resSeq <= prev_resSeq:
            chain = top.add_chain()
        new_res = top.add_residue(r.name, chain, resSeq=r.resSeq)
        for a in r.atoms:
            top.add_atom(a.name, a.element, new_res)
        prev_resSeq = r.resSeq
    t2 = md.Trajectory(t.xyz[:1], top)
    t2.save_pdb(str(out_pdb))
    return t.n_atoms, top.n_chains


def fix_his_protonation(pdb_path):
    """Relabels generic HIS residues to HID/HIE/HIP based on which of
    HD1/HE2 are actually present, since Amber force fields have no plain
    'HIS' template."""
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


def slice_last_frames(gro_path, xtc_path, pdb_path, out_dcd, n_last):
    """Writes a .dcd trajectory containing only the last `n_last` frames,
    with coordinates matching the no-bond PDB's topology/atom order."""
    t = md.load(str(xtc_path), top=str(gro_path))
    sliced = t[-n_last:]
    ref_top = md.load(str(pdb_path)).topology
    sliced = md.Trajectory(sliced.xyz, ref_top, unitcell_lengths=sliced.unitcell_lengths,
                            unitcell_angles=sliced.unitcell_angles)
    sliced.save_dcd(str(out_dcd))
    return sliced.n_frames


def prepare_common(system_name, out_dir):
    """Shared first stage for both systems: raw .gro/.xtc -> no-bond PDB
    (with HIS relabeled) + a sliced last-300-frame .dcd. This is the pair of
    files OpenGBSA's own config (<system>_config.yaml) points to."""
    cfg = SYSTEMS[system_name]
    raw_dir = HERE / cfg["raw_dir"]
    gro = raw_dir / "start.gro"
    xtc = raw_dir / "md.xtc"
    if not gro.exists() or not xtc.exists():
        sys.exit(f"Missing raw input {gro} / {xtc} -- run ./fetch_data.sh first.")

    out_dir.mkdir(parents=True, exist_ok=True)
    pdb_out = out_dir / "complex_chains.pdb"
    dcd_out = out_dir / "complex_chains_last300.dcd"

    n_atoms, n_chains = build_nobonds_pdb(gro, pdb_out)
    his_forms = fix_his_protonation(pdb_out)
    print(f"[{system_name}] {n_atoms} atoms in {n_chains} chains; "
          f"HIS residues relabeled: {his_forms or '(none present)'}")

    n_sliced = slice_last_frames(gro, xtc, pdb_out, dcd_out, N_LAST_FRAMES)
    print(f"[{system_name}] wrote {pdb_out.name} and {dcd_out.name} "
          f"({n_sliced} frames, last {N_LAST_FRAMES})")
    return pdb_out, dcd_out


# --- Independent Amber reference build (optional, --build-amber-reference) ---
# Two different topology-construction routes were used across the two
# systems during the original validation, kept here as distinct code paths
# rather than unified, since the README documents *why* each was chosen:
# tleap needed explicit HIS-tautomer/N-terminal-atom fixes and reorders
# atoms to its own canonical order; OpenMM+ParmEd needs no reorder step
# since it preserves the input PDB's atom order exactly.

def fix_nterm_h(pdb_path):
    lines = open(pdb_path).readlines()
    out, fixed = [], False
    for line in lines:
        if line.startswith("ATOM") and not fixed:
            atom_name = line[12:16].strip()
            resnum = line[22:26].strip()
            if atom_name == "H" and resnum == "1":
                line = line[:12] + " H1 " + line[16:]
                fixed = True
        out.append(line)
    open(pdb_path, "w").writelines(out)
    return fixed


def run_tleap(script_path):
    result = subprocess.run(["tleap", "-f", script_path.name], capture_output=True, text=True,
                             cwd=str(script_path.parent))
    if "Errors = 0" not in result.stdout:
        print(result.stdout[-3000:])
        raise RuntimeError(f"tleap failed for {script_path}")
    return result.stdout


def build_amber_reference_tleap(pdb_path, out_dir):
    """1GCQ's route: split into per-chain PDBs, fix each chain's N-terminal
    H1, run tleap, combine. tleap reorders atoms to its own canonical
    per-residue order, so a name-based reindex map is also written out for
    use when reordering the trajectory to match."""
    # mdtraj infers bonds by residue-template/geometry on *load* even from a
    # PDB with no CONECT records, and a plain atom_slice()+save_pdb() writes
    # those inferred bonds back out as real CONECT records -- which can
    # include spurious ones (seen here: an intra-residue CA-HB3 "bond") that
    # make tleap abort. Rebuilding a fresh, bond-free Topology per chain
    # (same technique as build_nobonds_pdb) avoids this.
    t = md.load(str(pdb_path))
    chain_pdbs = []
    for i, chain in enumerate(t.topology.chains):
        sub_top = md.Topology()
        new_chain = sub_top.add_chain()
        atom_idx = []
        for r in chain.residues:
            new_res = sub_top.add_residue(r.name, new_chain, resSeq=r.resSeq)
            for a in r.atoms:
                sub_top.add_atom(a.name, a.element, new_res)
                atom_idx.append(a.index)
        sub = md.Trajectory(t.xyz[:, atom_idx], sub_top)
        chain_pdb = out_dir / f"chain_{i}.pdb"
        sub.save_pdb(str(chain_pdb))
        # mdtraj normalizes HID/HIE/HIP back to the generic 'HIS' residue
        # name on load (it treats them as synonyms), so that tautomer info
        # doesn't survive the load-and-resave above -- re-derive it here
        # from each atom's actual HD1/HE2 content, same as build_nobonds_pdb.
        fix_his_protonation(chain_pdb)
        fix_nterm_h(chain_pdb)
        chain_pdbs.append(chain_pdb)

    tleap_script = out_dir / "tleap_complex.in"
    tleap_script.write_text(
        "source leaprc.protein.ff14SB\n"
        f"molA = loadpdb {chain_pdbs[0].name}\n"
        f"molB = loadpdb {chain_pdbs[1].name}\n"
        "combined = combine {molA molB}\n"
        "saveamberparm combined complex.prmtop complex.inpcrd\n"
        "quit\n"
    )
    run_tleap(tleap_script)
    print(f"Amber reference (tleap): {out_dir / 'complex.prmtop'}")


def build_amber_reference_openmm(pdb_path, out_dir, his68_hie_to_hid=False):
    """2OOB's route: OpenMM ForceField.createSystem + ParmEd
    openmm.load_topology, exported directly to .prmtop -- preserves the
    input PDB's atom order (no reindexing needed, unlike tleap)."""
    import parmed as pmd
    from openmm import app, unit

    work_pdb = pdb_path
    if his68_hie_to_hid:
        # Cross-simulation tautomer mismatch fix (2OOB-specific): if this
        # complex's His68 was built as HIE elsewhere (e.g. a standalone
        # receptor/ligand reference topology) but this complex uses HID,
        # convert HE2 -> HD1 at an estimated position along the CG-ND1
        # bond direction so both topologies agree (see README for why this
        # matters: MMPBSA.py rejects inconsistent tautomers across its
        # complex/receptor/ligand topologies).
        work_pdb = out_dir / "complex_his_fixed.pdb"
        _convert_hie_to_hid(pdb_path, work_pdb, resnum="68")

    topology = app.PDBFile(str(work_pdb)).topology
    positions = app.PDBFile(str(work_pdb)).positions
    ff = app.ForceField("amber14-all.xml")
    system = ff.createSystem(topology, nonbondedMethod=app.NoCutoff,
                              constraints=None, rigidWater=False)
    struct = pmd.openmm.load_topology(topology, system, xyz=positions.value_in_unit(unit.angstrom))
    pmd.tools.changeRadii(struct, "mbondi2").execute()
    struct.save(str(out_dir / "complex.prmtop"), overwrite=True)
    struct.save(str(out_dir / "complex.inpcrd"), overwrite=True)
    print(f"Amber reference (OpenMM+ParmEd): {out_dir / 'complex.prmtop'}")


def _convert_hie_to_hid(pdb_path, out_path, resnum):
    lines = open(pdb_path).readlines()
    he2 = None
    cg = nd1 = None
    for line in lines:
        if line.startswith("ATOM") and line[22:26].strip() == resnum:
            name = line[12:16].strip()
            coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if name == "HE2":
                he2 = coord
            elif name == "CG":
                cg = coord
            elif name == "ND1":
                nd1 = coord
    if he2 is None:
        raise RuntimeError(f"Residue {resnum} has no HE2 -- not in HIE form, nothing to convert.")
    if cg is None or nd1 is None:
        raise RuntimeError(f"Residue {resnum} missing CG/ND1 -- cannot estimate HD1 position.")
    direction = (cg - nd1)
    direction /= np.linalg.norm(direction)
    hd1_pos = nd1 + direction * 1.0  # ~1.0 A N-H bond length, geometrically approximate

    out = []
    for line in lines:
        if line.startswith("ATOM") and line[22:26].strip() == resnum:
            name = line[12:16].strip()
            resname = line[17:20].strip()
            if name == "HE2":
                line = (line[:12] + " HD1" + line[16:17] + "HID" + line[20:30]
                        + f"{hd1_pos[0]:8.3f}{hd1_pos[1]:8.3f}{hd1_pos[2]:8.3f}" + line[54:])
            elif resname == "HIE":
                line = line[:17] + "HID" + line[20:]
        out.append(line)
    open(out_path, "w").writelines(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--system", choices=sorted(SYSTEMS), required=True)
    ap.add_argument("--build-amber-reference", action="store_true",
                     help="Also build an independent Amber reference topology "
                          "(tleap for 1gcq, OpenMM+ParmEd for 2oob) for cross-checking.")
    args = ap.parse_args()

    out_dir = HERE / "prepared" / args.system
    pdb_out, dcd_out = prepare_common(args.system, out_dir)

    if args.build_amber_reference:
        ref_dir = out_dir / "amber_reference"
        ref_dir.mkdir(parents=True, exist_ok=True)
        if args.system == "1gcq":
            build_amber_reference_tleap(pdb_out, ref_dir)
        else:
            build_amber_reference_openmm(pdb_out, ref_dir)

    print(f"\nDone. Config for this system: {args.system}_config.yaml "
          f"(already points at prepared/{args.system}/).")


if __name__ == "__main__":
    main()
