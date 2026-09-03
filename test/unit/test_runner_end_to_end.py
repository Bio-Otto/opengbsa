#!/usr/bin/env python3
"""
End-to-end regression test for mmgbsa.runner.MMGBSARunner -- the class the
CLI (`opengbsa config.yaml`) actually instantiates, one layer above
GBSACalculator itself. test_charmm_loading.py's
test_charmm_run_end_to_end_produces_physical_energy already exercises the
same fixture directly through GBSACalculator.run(); this test exercises
the same fixture through the runner/config layer (binding_mode resolution,
YAML-shaped config dict, run_comprehensive dispatch) to catch regressions
in that layer specifically -- e.g. the runner.py:23-24 duplicate import,
config-key plumbing (forcefield_settings.charmm_params/charmm_coordinates
-> GBSACalculator kwargs), or binding_mode-to-selection resolution.

Fixture: test/data/charmm_namd_test/, see test_charmm_loading.py's
docstring for provenance.
"""
import math
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(HERE, "..", "data", "charmm_namd_test")
PSF_PATH = os.path.join(FIXTURE_DIR, "ctla4_P16_wat.psf")
PDB_PATH = os.path.join(FIXTURE_DIR, "ctla4_P16_wat.pdb")
PRM_PATH = os.path.join(FIXTURE_DIR, "par_all36m_prot.prm")
STR_PATH = os.path.join(FIXTURE_DIR, "toppar_water_ions_prot.str")
DCD_PATH = os.path.join(FIXTURE_DIR, "min_ctla4_P16.0.dcd")

pytestmark = pytest.mark.skipif(
    not os.path.exists(PSF_PATH), reason="CHARMM/NAMD test fixture not present"
)


def test_runner_end_to_end_via_config_dict(tmp_path):
    from mmgbsa.runner import MMGBSARunner

    config = {
        "input_files": {
            "complex_pdb": PSF_PATH,
            "trajectory": DCD_PATH,
            "solvated_topology": PDB_PATH,
        },
        "output_settings": {},
        "analysis_settings": {
            "binding_mode": "ppi",
            "receptor_selection": "chainid 1",
            "ligand_selection": "chainid 0",
            "temperature": 300.0,
            "gb_model": "OBC2",
            "salt_concentration": 0.15,
            "max_frames": 4,
            "use_cache": False,
            "run_entropy_analysis": False,
            "run_per_residue_decomposition": False,
        },
        "forcefield_settings": {
            "protein_forcefield": "charmm",
            "charmm_params": [PRM_PATH, STR_PATH],
            "charmm_coordinates": PDB_PATH,
        },
        "advanced_settings": {},
    }

    runner = MMGBSARunner(config, output_dir=str(tmp_path))
    results = runner.run_analysis()

    assert results is not None
    assert "mmgbsa" in results
    mean_binding = results["mmgbsa"]["mean_binding_energy"]
    assert math.isfinite(mean_binding)
    # Same generous bound as test_charmm_loading.py's direct-GBSACalculator
    # equivalent -- this isn't re-verifying the physics (that test already
    # does), just that the runner/config layer wires everything through
    # to the same non-blown-up result.
    assert abs(mean_binding) < 5000
