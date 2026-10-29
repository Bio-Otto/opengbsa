#!/usr/bin/env python3
"""
Entropy results and the systems they are computed from.

- entropy_method quasiharmonic is stored in results['entropy'] (and so in
  results_summary.yaml and final_report.txt) like the other methods.
- QHA uses real atomic masses in the atom order of the sliced trajectory
  (Coordinate Mode had hydrogen-repartitioned masses; the complex slice is
  receptor-then-ligand, which is not the System's order when the ligand is
  not the last block of a Native-Mode topology).
- Coordinate-Mode systems have no bond constraints and real masses, so
  normal-mode entropy works there (it stopped with 'needs a System without
  bond constraints').
"""
import math
from pathlib import Path

import mdtraj as md
import numpy as np
import openmm
import pytest
from openmm import app, unit

from mmgbsa.quasi_harmonic import QuasiHarmonicAnalysis

ROOT = Path(__file__).resolve().parents[2]
PRMTOP = ROOT / "test" / "data" / "6t1h_1151_comp" / "complex.prmtop"
XTC = ROOT / "examples" / "quickstart" / "trajectory.xtc"


def _small_trajectory(n_frames=40):
    top = md.Topology()
    res = top.add_residue("MOL", top.add_chain())
    for element in ("C", "H", "H", "O", "N"):
        top.add_atom(element, md.element.get_by_symbol(element), res)
    rng = np.random.default_rng(1)
    base = rng.normal(0, 0.15, (5, 3))
    xyz = base[None] + rng.normal(0, 0.01, (n_frames, 5, 3))
    return md.Trajectory(xyz.astype(np.float32), top)


def _system(masses):
    system = openmm.System()
    for m in masses:
        system.addParticle(m)
    return system


def test_qha_uses_the_masses_given():
    traj = _small_trajectory()
    real = [a.element.mass for a in traj.topology.atoms]
    repartitioned = [m + 3.0 if a.element.symbol == "H" else m for m, a in zip(real, traj.topology.atoms)]
    from_system = QuasiHarmonicAnalysis(traj[:], _system(real), 300.0, verbose=False).calculate_entropy()
    from_masses = QuasiHarmonicAnalysis(traj[:], _system(repartitioned), 300.0, verbose=False,
                                        masses=real).calculate_entropy()
    hmr = QuasiHarmonicAnalysis(traj[:], _system(repartitioned), 300.0, verbose=False).calculate_entropy()
    assert from_masses == pytest.approx(from_system, rel=1e-9)
    assert abs(hmr - from_system) > 1e-3
    with pytest.raises(ValueError):
        QuasiHarmonicAnalysis(traj[:], None, 300.0, verbose=False, masses=real[:3]).calculate_entropy()


@pytest.mark.skipif(not (PRMTOP.exists() and XTC.exists()), reason="6T1H quickstart fixture not present")
def test_quasiharmonic_entropy_is_stored(tmp_path):
    import yaml
    from mmgbsa.runner import MMGBSARunner
    config = {"input_files": {"complex_pdb": str(PRMTOP), "trajectory": str(XTC)},
              "analysis_settings": {"max_frames": 4, "run_entropy_analysis": True,
                                    "entropy_method": "quasiharmonic"},
              "platform_settings": {"preferred_platform": "CPU"}}
    out = tmp_path / "run"
    results = MMGBSARunner(config, output_dir=str(out)).run_analysis()
    entropy = results["entropy"]
    assert entropy["status"] == "completed" and entropy["method"] == "Quasi-harmonic entropy"
    assert entropy["entropy_penalty"] is None              # ligand only: no binding entropy
    assert math.isfinite(entropy["TS_ligand"]) and entropy["n_frames"] == 4
    summary = yaml.safe_load((out / "results_summary.yaml").read_text())
    assert summary["entropy"]["TS_ligand"] == pytest.approx(entropy["TS_ligand"])
    assert "Ligand TS" in (out / "final_report.txt").read_text()

    # the same (re-imaged: the raw ligand is split by the box in frame 3) ligand
    # frames with real masses, computed directly
    from mmgbsa.trajectory import reimage_molecules
    traj = reimage_molecules(md.load(str(XTC), top=str(PRMTOP))[:4])
    ligand = traj.atom_slice(traj.topology.select("resname LIG"))
    masses = [a.element.mass for a in ligand.topology.atoms]
    direct = QuasiHarmonicAnalysis(ligand, None, 300.0, verbose=False, masses=masses).calculate_entropy()
    assert entropy["TS_ligand"] == pytest.approx(direct, abs=1e-3)


def _peptide_dimer(tmp_path):
    """Two ACE-ALA-NH2 peptides 6 A apart, chains A and B, with a 3-frame trajectory."""
    import openff.toolkit
    src = Path(openff.toolkit.__file__).parent / "data" / "proteins" / "ace-ala-nh2.pdb"
    if not src.exists():
        pytest.skip("OpenFF test peptide not available")
    pdb = app.PDBFile(str(src))
    top, pos = app.Topology(), []
    for shift, cid in ((0.0, "A"), (0.6, "B")):
        chain = top.addChain(cid)
        atoms = {}
        for res in pdb.topology.residues():
            r = top.addResidue(res.name, chain)
            for a in res.atoms():
                atoms[a] = top.addAtom(a.name, a.element, r)
        for a, b in pdb.topology.bonds():
            top.addBond(atoms[a], atoms[b])
        pos += [openmm.Vec3(p[0] + shift, p[1], p[2]) for p in pdb.positions.value_in_unit(unit.nanometer)]
    path = tmp_path / "dimer.pdb"
    with open(path, "w") as f:
        app.PDBFile.writeFile(top, pos * unit.nanometer, f, keepIds=True)
    t = md.load(str(path))
    rng = np.random.default_rng(0)
    xyz = np.repeat(t.xyz, 3, axis=0) + rng.normal(0, 0.002, (3,) + t.xyz.shape[1:])
    md.Trajectory(xyz, t.topology, unitcell_lengths=np.full((3, 3), 5.0),
                  unitcell_angles=np.full((3, 3), 90.0)).save_dcd(str(tmp_path / "dimer.dcd"))
    return path, tmp_path / "dimer.dcd"


def test_coordinate_mode_normal_modes(tmp_path):
    from mmgbsa.runner import MMGBSARunner
    pdb, dcd = _peptide_dimer(tmp_path)
    config = {"input_files": {"complex_pdb": str(pdb), "trajectory": str(dcd)},
              "analysis_settings": {"binding_mode": "ppi", "max_frames": 2, "run_entropy_analysis": True,
                                    "entropy_method": "nmode", "nmode_frames": 1},
              "platform_settings": {"preferred_platform": "CPU"}}
    runner = MMGBSARunner(config, output_dir=str(tmp_path / "run"))
    results = runner.run_analysis()
    for name in ("complex_system", "protein_system", "ligand_system"):
        system = runner.calculator.systems[name]
        assert system.getNumConstraints() == 0, name
        masses = [system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(system.getNumParticles())]
        assert min(masses) == pytest.approx(1.008, abs=1e-3), name       # hydrogen, not 4 amu
    entropy = results["entropy"]
    assert entropy["method"] == "Normal-mode entropy" and entropy["n_converged"] == 1
    assert math.isfinite(entropy["entropy_penalty"])
