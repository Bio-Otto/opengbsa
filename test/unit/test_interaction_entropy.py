#!/usr/bin/env python3
"""
Interaction-energy entropy estimates (mmgbsa.interaction_entropy).

Up to 0.0.18 the interaction entropy was computed from the fluctuation of the
binding enthalpy, whose GB and SA terms damp the fluctuation. Duan et al.
(J Am Chem Soc 2016) and gmx_MMPBSA use the gas-phase interaction energy
dE_vdW + dE_elec. These tests pin the formulas to gmx_MMPBSA's own
implementation and check that the pipeline feeds them that energy.
"""
import math
import os

import numpy as np
import pandas as pd
import pytest

from mmgbsa.interaction_entropy import (
    C2_SIGMA_LIMIT, GAS_CONSTANT, IE_SIGMA_LIMIT, c2_entropy, entropy_estimates,
    entropy_result, interaction_energies, interaction_entropy,
)

# Deterministic test series; the reference values come from gmx_MMPBSA's
# InteractionEntropyCalc (ie_value) and C2EntropyCalc (c2data) at 300 K on
# the same array. They differ from ours in the 8th digit only through the
# gas constant (gmx_MMPBSA: 0.00198720425864083 kcal/(mol K)).
SERIES = -50.0 + 2.5 * np.sin(np.arange(200) * 0.37) + 1.5 * np.cos(np.arange(200) * 1.91)
GMX_MMPBSA_IE = 2.2378411997895853
GMX_MMPBSA_C2 = 3.57258426657035
GMX_MMPBSA_SIGMA = 2.0638974783959343


def test_matches_gmx_mmpbsa():
    result = entropy_estimates(SERIES, 300.0)
    assert result['interaction_entropy'] == pytest.approx(GMX_MMPBSA_IE, rel=1e-6)
    assert result['c2_entropy'] == pytest.approx(GMX_MMPBSA_C2, rel=1e-6)
    assert result['sigma_interaction_energy'] == pytest.approx(GMX_MMPBSA_SIGMA, rel=1e-12)


def test_constant_energy_has_no_entropy():
    E = np.full(50, -42.0)
    assert interaction_entropy(E, 300.0) == pytest.approx(0.0, abs=1e-12)
    assert c2_entropy(E, 300.0) == 0.0


def test_gaussian_ie_approaches_c2():
    # For Gaussian energies the cumulant series ends at second order, so
    # IE converges to C2 = sigma^2 / (2 kT)
    E = np.random.default_rng(7).normal(-30.0, 1.0, 400_000)
    assert interaction_entropy(E, 300.0) == pytest.approx(1.0 / (2 * GAS_CONSTANT * 300.0), rel=0.02)


def test_reliability_limits_are_ekberg_ryde():
    assert IE_SIGMA_LIMIT == pytest.approx(15 / 4.184)
    assert C2_SIGMA_LIMIT == pytest.approx(25 / 4.184)
    small = entropy_estimates(SERIES, 300.0)
    assert small['interaction_entropy_reliable'] and small['c2_entropy_reliable']
    large = entropy_estimates(SERIES * 3.0, 300.0)  # sigma about 6.2 kcal/mol
    assert not large['interaction_entropy_reliable'] and not large['c2_entropy_reliable']


def test_interaction_energy_excludes_solvation():
    frames = pd.DataFrame({
        'delta_vdw': [-40.0, -42.0, -38.0, -41.0],
        'delta_elec': [-10.0, -12.0, -9.0, -11.0],
        'delta_gb': [20.0, 31.0, 12.0, 25.0],
        'binding_energy': [-30.0, -23.0, -35.0, -27.0],
    })
    np.testing.assert_allclose(interaction_energies(frames), [-50.0, -54.0, -47.0, -52.0])
    result = entropy_result(interaction_energies(frames), 'c2', 300.0)
    assert result['label'] == 'C2 Entropy'
    assert result['value'] == pytest.approx(np.var([-50.0, -54.0, -47.0, -52.0]) / (2 * GAS_CONSTANT * 300.0))
    with pytest.raises(ValueError, match='delta_elec'):
        interaction_energies(frames.drop(columns='delta_elec'))


def test_unknown_method_and_too_few_frames():
    with pytest.raises(ValueError, match='Unknown'):
        entropy_result(SERIES, 'quasiharmonic', 300.0)
    with pytest.raises(ValueError, match='at least 2'):
        interaction_entropy([1.0], 300.0)


HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(HERE, "..", "data", "charmm_namd_test")
PSF_PATH = os.path.join(FIXTURE_DIR, "ctla4_P16_wat.psf")


@pytest.mark.skipif(not os.path.exists(PSF_PATH), reason="CHARMM/NAMD test fixture not present")
@pytest.mark.parametrize("entropy_method", ["c2", None])
def test_pipeline_uses_interaction_energy(tmp_path, entropy_method):
    """Both entropy paths -- GBSACalculator.save_results (entropy_method) and
    the runner's entropy step (run_entropy_analysis, interaction entropy by
    default) -- report the estimate of dE_vdW + dE_elec from the written
    frame table, and the MM/GBSA results that feed the HTML report carry it."""
    from mmgbsa.runner import MMGBSARunner

    config = {
        "input_files": {
            "complex_pdb": PSF_PATH,
            "trajectory": os.path.join(FIXTURE_DIR, "min_ctla4_P16.0.dcd"),
            "solvated_topology": os.path.join(FIXTURE_DIR, "ctla4_P16_wat.pdb"),
        },
        "output_settings": {},
        "analysis_settings": {
            "binding_mode": "ppi",
            "receptor_selection": "chainid 1",
            "ligand_selection": "chainid 0",
            "temperature": 310.0,
            "gb_model": "OBC2",
            "salt_concentration": 0.15,
            "max_frames": 4,
            "use_cache": False,
            "run_entropy_analysis": True,
            "run_per_residue_decomposition": False,
        },
        "forcefield_settings": {
            "protein_forcefield": "charmm",
            "charmm_params": [os.path.join(FIXTURE_DIR, "par_all36m_prot.prm"),
                              os.path.join(FIXTURE_DIR, "toppar_water_ions_prot.str")],
            "charmm_coordinates": os.path.join(FIXTURE_DIR, "ctla4_P16_wat.pdb"),
        },
        "advanced_settings": {},
    }
    if entropy_method:
        config["analysis_settings"]["entropy_method"] = entropy_method
    results = MMGBSARunner(config, output_dir=str(tmp_path)).run_analysis()

    frames = pd.read_csv(results["mmgbsa"]["output_file"])
    E = frames["delta_vdw"] + frames["delta_elec"]
    method = entropy_method or "interaction"
    expected = (c2_entropy if method == "c2" else interaction_entropy)(E, 310.0)
    assert math.isfinite(expected) and expected > 0

    core = results["mmgbsa"]["entropy"]
    assert core["method"] == method and core["value"] == pytest.approx(expected, rel=1e-9)
    assert results["mmgbsa"]["delta_g"] == pytest.approx(
        results["mmgbsa"]["mean_binding_energy"] + expected, rel=1e-9)
    assert results["entropy"]["entropy_penalty"] == pytest.approx(expected, rel=1e-9)
    assert results["entropy"]["sigma_interaction_energy"] == pytest.approx(float(np.std(E)), rel=1e-9)
