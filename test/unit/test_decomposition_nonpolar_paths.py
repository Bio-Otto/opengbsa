#!/usr/bin/env python3
"""Regression test: the parallel CPU decomposition path must report the per-residue
non-polar (surface-area) term, as the serial path does.

Up to 0.0.17 the parallel path (decomposition_platform CPU with more than one worker)
computed only van der Waals, electrostatic and polar terms, so its
per_residue_detailed.csv had no `nonpolar_solvation` column and per-residue solvation
was polar-only. Both paths now give the same non-polar term.

Fixture: test/data/6t1h_1151_comp/ (Amber topology, protein-ligand), 3 frames; the
trajectory is not tracked in git, so the test is skipped where it is absent.
"""
import os

import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(HERE, "..", "data", "6t1h_1151_comp")
PRMTOP = os.path.join(FIXTURE_DIR, "complex.prmtop")
XTC = os.path.join(FIXTURE_DIR, "md_complex_prod.xtc")
GRO = os.path.join(FIXTURE_DIR, "replica_1", "gro", "complex_solv_ions.gro")

pytestmark = pytest.mark.skipif(
    not all(map(os.path.exists, (PRMTOP, XTC, GRO))), reason="6T1H trajectory fixture not present"
)


def _run(tmp_path, platform, n_jobs):
    from mmgbsa.runner import MMGBSARunner

    out = tmp_path / f"{platform}_{n_jobs}"
    config = {
        "input_files": {"complex_pdb": PRMTOP, "trajectory": XTC, "solvated_topology": GRO},
        "output_settings": {},
        "analysis_settings": {
            "temperature": 300.0, "gb_model": "OBC2", "salt_concentration": 0.15,
            "max_frames": 3, "decomp_frames": 3, "use_cache": False,
            "run_entropy_analysis": False, "run_per_residue_decomposition": True,
        },
        "forcefield_settings": {"protein_forcefield": "amber"},
        "platform_settings": {"preferred_platform": "CPU", "decomposition_platform": platform},
        "advanced_settings": {"n_jobs": n_jobs},
    }
    MMGBSARunner(config, output_dir=str(out)).run_analysis()
    files = list(out.rglob("per_residue_detailed.csv"))
    assert files, "no per_residue_detailed.csv written"
    return pd.read_csv(files[0]).set_index("residue_id").sort_index()


def test_parallel_and_serial_paths_report_the_same_nonpolar_term(tmp_path):
    parallel = _run(tmp_path, "CPU", 2)        # parallel worker path
    serial = _run(tmp_path, "Reference", 1)    # serial path
    for df in (parallel, serial):
        assert "nonpolar_solvation" in df.columns
    assert (serial.nonpolar_solvation.abs() > 1e-6).any()   # residues at the binding site
    common = parallel.index.intersection(serial.index)
    assert len(common) == len(serial)
    pd.testing.assert_series_equal(parallel.loc[common, "nonpolar_solvation"],
                                   serial.loc[common, "nonpolar_solvation"], atol=1e-6, check_names=False)
