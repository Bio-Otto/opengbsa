#!/usr/bin/env python
"""Build Amber topologies for the VRC01 system so MMPBSA.py can be run on it.

MMPBSA.py needs four topologies (complex, receptor, ligand, plus the solvated
complex it strips from) and they must describe exactly the same atoms, in the
same order, as the trajectory OpenGBSA was given -- otherwise the two tools are
not comparing the same system.

Approach: load the CHARMM .psf with ParmEd, attach the CHARMM36m parameters,
strip water and ions, then write Amber-format topologies via ParmEd's
`save(..., format='amber')`, which performs the CHARMM->Amber conversion
(the `chamber` operation). CMAP terms survive the conversion; Amber's sander
supports them.

Run from this directory:  python prepare_amber.py
"""
import sys
import pathlib
import parmed as pmd

HERE = pathlib.Path(__file__).parent
PARAMS = [HERE / "params" / "par_all36m_prot.prm",
          HERE / "params" / "toppar_water_ions_prot.str"]

# segid -> role.  Verified against the deposit: contiguous atom blocks.
ANTIGEN_SEG = {"G120"}
ANTIBODY_SEG = {"HC", "LC"}
SOLVENT_SEG = {"WAT", "ION1", "ION2"}


def build(system: str, outdir: pathlib.Path) -> None:
    psf = HERE / "raw" / f"{system}.psf"
    pdb = HERE / "raw" / f"{system}.pdb"
    if not psf.exists():
        sys.exit(f"missing {psf}")

    print(f"\n=== {system} ===")
    st = pmd.load_file(str(psf))
    st.load_parameters(pmd.charmm.CharmmParameterSet(*[str(p) for p in PARAMS]))
    # coordinates must come from the PDB; a PSF carries none
    st.coordinates = pmd.load_file(str(pdb)).coordinates
    print(f"  loaded {len(st.atoms)} atoms, {len(st.residues)} residues")

    # index the three groups by segid before any stripping reorders things
    def mask(segs):
        return [i for i, a in enumerate(st.atoms) if a.residue.segid in segs]

    ag, ab, solv = mask(ANTIGEN_SEG), mask(ANTIBODY_SEG), mask(SOLVENT_SEG)
    print(f"  antigen {len(ag)}  antibody {len(ab)}  solvent {len(solv)}")

    outdir.mkdir(parents=True, exist_ok=True)

    # Dry complex only.  We deliberately do NOT write a solvated topology for
    # MMPBSA.py's -sp argument: cpptraj's strip mask is Amber-centric
    # (:WAT,Cl*,CIO,...,NA,CL) and does not match this system's CHARMM residue
    # names (TIP3, SOD, CLA), so all 10,224 solvent atoms survive the strip and
    # sander dies with "NATOM mismatch in coord and topology files". Instead the
    # trajectory is written solvent-free up front (see make_dry_traj.py) so it
    # matches the dry topology exactly.
    dry = st[[i for i in range(len(st.atoms)) if i not in set(solv)]]

    # receptor = antigen, ligand = antibody Fv, indexed within the dry complex
    n_ag = len(ag)
    rec = dry[list(range(n_ag))]
    lig = dry[list(range(n_ag, n_ag + len(ab)))]

    # A CHARMM .psf carries no GB radii, and ParmEd's CHARMM->Amber conversion
    # does not invent them: every atom comes out with solvent_radius = 0 and
    # screen = 0, which makes sander's igb=5 return EGB = NaN (silently -- the
    # run completes and only DELTA TOTAL shows nan). Assign the radius set the
    # GB model requires, exactly as OpenGBSA does internally. OBC2 (igb=5)
    # needs mbondi2.
    # changeRadii needs an AmberParm, not a CharmmPsfFile, so write first and
    # reload the Amber topology before assigning radii.
    for name, struct in (("complex", dry), ("receptor", rec), ("ligand", lig)):
        path = outdir / f"{system}_{name}.prmtop"
        struct.save(str(path), format="amber", overwrite=True)
        parm = pmd.load_file(str(path))
        pmd.tools.changeRadii(parm, "mbondi2").execute()
        n_zero = sum(1 for a in parm.atoms if a.solvent_radius == 0)
        assert n_zero == 0, f"{name}: {n_zero} atoms still have zero GB radius"
        parm.save(str(path), overwrite=True)
        print(f"  {name}: {len(parm.atoms)} atoms, radii=mbondi2")

    # sanity: the split must be exact, not approximate
    assert len(rec.atoms) + len(lig.atoms) == len(dry.atoms), "split lost atoms"
    print("  OK")


if __name__ == "__main__":
    systems = sys.argv[1:] or ["wt"]
    for s in systems:
        build(s, HERE / "amber" / s)
