#!/usr/bin/env python
"""Generate per-system OpenGBSA configs and solvent-free Amber trajectories.

One config per mutant, derived from vrc01_wt_config.yaml. The residue-index
selections are re-derived per system rather than copied, because an alanine
mutation removes side-chain atoms and the atom counts therefore differ between
systems -- assuming wt's indices would silently mis-split the mutants.
"""
import pathlib
import sys

import mdtraj as md
import parmed as pmd

HERE = pathlib.Path(__file__).parent
SYSTEMS = ["wt", "g54a", "d99a", "n100aa", "w100ba"]
# experimental dG, Conti 2022 Antibodies 11(3):51, Table 1, RSC3/VRC01 column
EXP_DG = {"wt": -10.51, "g54a": -11.78, "d99a": -10.45,
          "n100aa": -9.40, "w100ba": -6.09}

ANTIGEN, ANTIBODY = {"G120"}, {"HC", "LC"}
SOLVENT = {"WAT", "ION1", "ION2"}
FRAME_START, FRAME_END = 4301, 4601


def ranges(system):
    """Residue-index ranges and atom counts, measured from this system's files."""
    st = pmd.load_file(str(HERE / "raw" / f"{system}.psf"))
    top = md.load(str(HERE / "raw" / f"{system}.pdb")).topology

    seg_res = {}
    for i, a in enumerate(st.atoms):
        sg = a.residue.segid
        r = top.atom(i).residue.index
        d = seg_res.setdefault(sg, [r, r])
        d[0], d[1] = min(d[0], r), max(d[1], r)

    ag = seg_res["G120"]
    ab = [seg_res["HC"][0], seg_res["LC"][1]]
    keep = [i for i, a in enumerate(st.atoms) if a.residue.segid not in SOLVENT]

    n_ag = sum(1 for a in st.atoms if a.residue.segid in ANTIGEN)
    n_ab = sum(1 for a in st.atoms if a.residue.segid in ANTIBODY)

    # the selection must reproduce the segid-derived counts exactly
    sel_ag = f"resid {ag[0]} to {ag[1]}"
    sel_ab = f"resid {ab[0]} to {ab[1]}"
    got_ag, got_ab = len(top.select(sel_ag)), len(top.select(sel_ab))
    assert got_ag == n_ag, f"{system}: antigen {got_ag} != {n_ag}"
    assert got_ab == n_ab, f"{system}: antibody {got_ab} != {n_ab}"
    return sel_ag, sel_ab, n_ag, n_ab, keep


def main():
    wt_cfg = (HERE / "vrc01_wt_config.yaml").read_text()
    for s in (sys.argv[1:] or SYSTEMS):
        sel_ag, sel_ab, n_ag, n_ab, keep = ranges(s)
        print(f"{s:8} antigen {sel_ag:18} {n_ag:5} atoms | "
              f"antibody {sel_ab:20} {n_ab:5} atoms | dry {len(keep)}")

        if s != "wt":
            cfg = (wt_cfg
                   .replace("raw/wt.psf", f"raw/{s}.psf")
                   .replace("raw/wt.pdb", f"raw/{s}.pdb")
                   .replace("raw/wt.dcd", f"raw/{s}.dcd")
                   .replace("results/wt", f"results/{s}")
                   .replace("receptor_selection: 'resid 0 to 350'",
                            f"receptor_selection: '{sel_ag}'")
                   .replace("ligand_selection: 'resid 351 to 570'",
                            f"ligand_selection: '{sel_ab}'")
                   .replace("dG = -10.51 kcal/mol for wild type",
                            f"dG = {EXP_DG[s]} kcal/mol for {s}"))
            (HERE / f"vrc01_{s}_config.yaml").write_text(cfg)

        # solvent-free NetCDF for the Amber arm
        out = HERE / "amber_run" / s
        out.mkdir(parents=True, exist_ok=True)
        nc = out / f"{s}_dry_300.nc"
        if not nc.exists():
            t = md.load_dcd(str(HERE / "raw" / f"{s}.dcd"),
                            top=str(HERE / "raw" / f"{s}.pdb"))
            t[FRAME_START:FRAME_END].atom_slice(keep).save_netcdf(str(nc))
            print(f"         wrote {nc.name} ({nc.stat().st_size/1e6:.0f} MB)")


if __name__ == "__main__":
    main()
