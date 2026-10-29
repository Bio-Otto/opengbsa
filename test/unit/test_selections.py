#!/usr/bin/env python3
"""
Receptor and ligand selections.

- ligand_selection alone: the receptor is every other atom (until 0.0.22 the
  receptor positions included the ligand and every frame failed).
- receptor_selection alone: the ligand is every other atom.
- Both: they must not overlap; atoms in neither are left out of the
  calculation, in Native Mode by reducing the topology and the trajectory.
- binding_mode ppi needs exactly two chains; a .parm7 receptor/ligand topology
  is used and any non-Amber file is an error instead of being ignored.
"""
import math
import os
from pathlib import Path

import mdtraj as md
import numpy as np
import pytest

from mmgbsa.mmgbsa_core import GBSACalculator, index_selection

ROOT = Path(__file__).resolve().parents[2]
PRMTOP = ROOT / "test" / "data" / "6t1h_1151_comp" / "complex.prmtop"
XTC = ROOT / "examples" / "quickstart" / "trajectory.xtc"
needs_6t1h = pytest.mark.skipif(not (PRMTOP.exists() and XTC.exists()),
                                reason="6T1H quickstart fixture not present")


def _topology(chains):
    """An MDTraj topology with `chains` = [(n_residues, resname), ...], 2 atoms per residue."""
    top = md.Topology()
    for n_res, resname in chains:
        chain = top.add_chain()
        for _ in range(n_res):
            res = top.add_residue(resname, chain)
            top.add_atom("C", md.element.carbon, res)
            top.add_atom("H", md.element.hydrogen, res)
    return top


def test_index_selection_is_exact():
    top = _topology([(10, "ALA")])
    for idx in ([0], [3, 4, 5, 9], [0, 2, 4, 6], list(range(20))):
        sel = index_selection(idx)
        assert list(top.select(sel)) == idx, sel
    assert index_selection([5, 1, 2, 3]) == "index 1 to 3 or index 5"


def test_overlapping_selections_are_an_error():
    calc = GBSACalculator()
    top = _topology([(5, "ALA"), (3, "GLY")])
    with pytest.raises(ValueError, match="overlap in 6 atoms"):
        calc.split_by_selections(top, "chainid 1", "all", "the test complex")


def test_empty_selection_is_an_error():
    calc = GBSACalculator()
    top = _topology([(5, "ALA"), (3, "GLY")])
    with pytest.raises(ValueError, match="ligand_selection 'chainid 4' matches no atom"):
        calc.split_by_selections(top, "chainid 4", "chainid 0", "the test complex")
    with pytest.raises(ValueError, match="lost its chain records"):
        calc.split_by_selections(top, "chainid 1", "chainid 7", "x.pdb (converted from x.tpr)")


def test_atoms_in_neither_selection_are_reported_not_water():
    calc = GBSACalculator()
    top = _topology([(5, "ALA"), (3, "GLY"), (2, "SER"), (4, "HOH")])
    lig, rec, left_out = calc.split_by_selections(top, "chainid 1", "chainid 0", "the test complex")
    assert lig == list(range(10, 16)) and rec == list(range(10))
    assert left_out == list(range(16, 20))          # the SER chain; water is not counted
    assert any("left out" in a and "chain 2" in a for a in calc.physics_assumptions)


@pytest.mark.parametrize("chains,match", [([(5, "ALA"), (3, "HOH")], "single chain"),
                                          ([(5, "ALA"), (3, "GLY"), (2, "SER")], "3 chains")])
def test_ppi_needs_exactly_two_chains(chains, match):
    with pytest.raises(ValueError, match=match):
        GBSACalculator().check_two_chains(_topology(chains), "the test complex")
    GBSACalculator().check_two_chains(_topology([(5, "ALA"), (3, "GLY"), (4, "HOH")]), "ok")


@pytest.mark.parametrize("name,ok", [("rec.prmtop", True), ("rec.parm7", True), ("REC.PARM7", True),
                                     ("rec.new", True), ("rec.top", False), ("rec.pdb", False)])
def test_explicit_receptor_topology_extensions(name, ok):
    if ok:
        assert GBSACalculator.explicit_amber_topology(name, "receptor_topology") == name
    else:
        with pytest.raises(ValueError, match="receptor_topology"):
            GBSACalculator.explicit_amber_topology(name, "receptor_topology")
    assert GBSACalculator.explicit_amber_topology(None, "receptor_topology") is None


def test_runner_maps_a_lone_receptor_selection_to_its_complement():
    from mmgbsa.runner import MMGBSARunner
    runner = MMGBSARunner({"input_files": {"complex_pdb": "c.prmtop", "trajectory": "t.xtc"},
                           "analysis_settings": {"receptor_selection": "chainid 0"}})
    assert runner._resolve_binding_mode_selections() == (None, "not (chainid 0)", False)
    runner.config["analysis_settings"].update(binding_mode="ppi", receptor_selection=None)
    assert runner._resolve_binding_mode_selections() == ("chainid 0", "chainid 1", True)
    runner.config["analysis_settings"].update(receptor_selection="chainid 1", ligand_selection="chainid 0")
    assert runner._resolve_binding_mode_selections() == ("chainid 1", "chainid 0", False)


def _run(tmp_path, name, complex_pdb=PRMTOP, trajectory=XTC, **analysis):
    from mmgbsa.runner import MMGBSARunner
    config = {"input_files": {"complex_pdb": str(complex_pdb), "trajectory": str(trajectory)},
              "analysis_settings": dict(max_frames=1, **analysis),
              "platform_settings": {"preferred_platform": "CPU"}}
    results = MMGBSARunner(config, output_dir=str(tmp_path / name)).run_analysis()
    assert results, name
    return results["mmgbsa"]["binding_energies"]


@needs_6t1h
def test_ligand_selection_alone_gives_the_standard_result(tmp_path):
    standard = _run(tmp_path, "standard")
    by_selection = _run(tmp_path, "ligand_only", ligand_selection="resname LIG")
    assert np.allclose(by_selection, standard, atol=1e-3)


@needs_6t1h
def test_ppi_on_a_single_chain_topology_stops(tmp_path):
    from mmgbsa.runner import MMGBSARunner
    config = {"input_files": {"complex_pdb": str(PRMTOP), "trajectory": str(XTC)},
              "analysis_settings": {"binding_mode": "ppi", "max_frames": 1}}
    with pytest.raises(ValueError, match="single chain"):
        MMGBSARunner(config, output_dir=str(tmp_path / "ppi")).run_analysis()


@needs_6t1h
def test_atoms_in_neither_selection_match_a_reduced_topology(tmp_path):
    """Receptor = residues 0-119, ligand = LIG: residues 120-240 are left out. The
    result must equal a run on a topology and trajectory reduced by hand."""
    import parmed as pmd
    reduced_run = _run(tmp_path, "selections", ligand_selection="resname LIG",
                       receptor_selection="resid 0 to 119", sa_model="ACE")

    structure = pmd.load_file(str(PRMTOP))
    keep = [a.idx for a in structure.atoms if a.residue.idx < 120 or a.residue.name == "LIG"]
    structure.strip("!(:1-120,LIG)")
    assert len(structure.atoms) == len(keep)
    reduced_top = tmp_path / "reduced.prmtop"
    structure.save(str(reduced_top))
    traj = md.load(str(XTC), top=str(PRMTOP))[:1].atom_slice(keep)
    reduced_dcd = tmp_path / "reduced.dcd"
    traj.save_dcd(str(reduced_dcd))
    reference = _run(tmp_path, "reference", complex_pdb=reduced_top, trajectory=reduced_dcd,
                     sa_model="ACE")
    assert math.isfinite(reduced_run[0])
    assert np.allclose(reduced_run, reference, atol=1e-3), (reduced_run, reference)
    assert (tmp_path / "selections" / "distilled.xtc").exists()
