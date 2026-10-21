#!/usr/bin/env python3
"""Regression tests for the per-residue decomposition paths.

1. The parallel CPU path (decomposition_platform CPU with more than one worker)
   must report the per-residue non-polar (surface-area) term, as the serial path
   does. Up to 0.0.17 it computed van der Waals, electrostatic and polar terms
   only, so per_residue_detailed.csv had no `nonpolar_solvation` column.
   Fixture: test/data/6t1h_1151_comp/ (Amber topology, protein-ligand), 3 frames.

2. A CHARMM protein-protein system with SOD/CLA ions and a selection-defined
   ligand must decompose on both paths with the receptor GB self-energy
   correction. Up to 0.0.17 the decomposition stripped only Amber/GROMACS ion
   names, so the trajectory kept 29 ions the system did not have: the serial
   path failed on every frame and the parallel path fell back to the distance
   heuristic for the polar term. It also took the receptor as "every residue
   not named like the guessed ligand" (a CHARMM HSD), so the self-energy
   correction was skipped. Fixture: test/data/charmm_namd_test/ (CTLA4 and a
   peptide, NAMD .psf), 2 frames.

The trajectories are not tracked in git, so the tests are skipped where they
are absent.
"""
import os

import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

T1H = os.path.join(DATA, "6t1h_1151_comp")
T1H_FILES = (os.path.join(T1H, "complex.prmtop"), os.path.join(T1H, "md_complex_prod.xtc"),
             os.path.join(T1H, "replica_1", "gro", "complex_solv_ions.gro"))

CT = os.path.join(DATA, "charmm_namd_test")
CT_FILES = (os.path.join(CT, "ctla4_P16_wat.psf"), os.path.join(CT, "min_ctla4_P16.0.dcd"),
            os.path.join(CT, "ctla4_P16_wat.pdb"))
CT_PARAMS = [os.path.join(CT, "par_all36m_prot.prm"), os.path.join(CT, "toppar_water_ions_prot.str")]

TERMS = ["vdw", "electrostatic", "solvation", "nonpolar_solvation", "total"]


def _present(files):
    return all(map(os.path.exists, files))


def _run(out, files, platform, n_jobs, frames, analysis=None, forcefield=None):
    from mmgbsa.runner import MMGBSARunner

    complex_top, traj, solvated = files
    config = {
        "input_files": {"complex_pdb": complex_top, "trajectory": traj, "solvated_topology": solvated},
        "output_settings": {},
        "analysis_settings": {
            "temperature": 300.0, "gb_model": "OBC2", "salt_concentration": 0.15,
            "max_frames": frames, "decomp_frames": frames, "use_cache": False,
            "run_entropy_analysis": False, "run_per_residue_decomposition": True, **(analysis or {}),
        },
        "forcefield_settings": forcefield or {"protein_forcefield": "amber"},
        "platform_settings": {"preferred_platform": "CPU", "decomposition_platform": platform},
        "advanced_settings": {"n_jobs": n_jobs},
    }
    MMGBSARunner(config, output_dir=str(out)).run_analysis()
    files = list(out.rglob("per_residue_detailed.csv"))
    assert files, "no per_residue_detailed.csv written"
    return pd.read_csv(files[0]).set_index("residue_id").sort_index()


@pytest.mark.skipif(not _present(T1H_FILES), reason="6T1H trajectory fixture not present")
def test_parallel_and_serial_paths_report_the_same_nonpolar_term(tmp_path):
    parallel = _run(tmp_path / "parallel", T1H_FILES, "CPU", 2, 3)
    serial = _run(tmp_path / "serial", T1H_FILES, "Reference", 1, 3)
    for df in (parallel, serial):
        assert "nonpolar_solvation" in df.columns
    assert (serial.nonpolar_solvation.abs() > 1e-6).any()   # residues at the binding site
    common = parallel.index.intersection(serial.index)
    assert len(common) == len(serial)
    pd.testing.assert_series_equal(parallel.loc[common, "nonpolar_solvation"],
                                   serial.loc[common, "nonpolar_solvation"], atol=1e-6, check_names=False)
    # The serial path also writes the non-polar term frame by frame, and its frame mean
    # is the per-residue value.
    frames = pd.read_csv(next((tmp_path / "serial").rglob("frame_by_frame_decomposition.csv")))
    for rid, value in serial.nonpolar_solvation.items():
        name, num, chain = rid.split("_")
        assert frames[f"{name}{num}_{chain}_nonpolar_solvation"].mean() == pytest.approx(value, abs=1e-6)


@pytest.mark.skipif(not _present(CT_FILES), reason="CHARMM/NAMD test fixture not present")
def test_charmm_ppi_decomposition_on_both_paths(tmp_path, capfd):
    kw = dict(analysis={"binding_mode": "ppi", "receptor_selection": "chainid 1", "ligand_selection": "chainid 0"},
              forcefield={"protein_forcefield": "charmm", "charmm_params": CT_PARAMS,
                          "charmm_coordinates": CT_FILES[2]})
    serial = _run(tmp_path / "serial", CT_FILES, "Reference", 1, 2, **kw)
    parallel = _run(tmp_path / "parallel", CT_FILES, "CPU", 2, 2, **kw)
    log = capfd.readouterr()
    log = log.out + log.err
    assert "keeping 2030" in log                                    # 29 SOD/CLA ions stripped
    assert log.count("Protein-only GB context prepared (1758 atoms)") == 2
    assert "Rigorous GB decomposition failed" not in log
    assert len(serial) == 118                                       # CTLA4 residues
    pd.testing.assert_frame_equal(serial[TERMS], parallel.loc[serial.index, TERMS], atol=1e-6)
