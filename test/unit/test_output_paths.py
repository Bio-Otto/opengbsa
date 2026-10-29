#!/usr/bin/env python3
"""
A run writes only into its own folder: the .trj -> .mdcrd alias and the
sanitized copy of a GROMACS topology go there, and input files are never
changed (before, the alias was created and the .top rewritten next to the
inputs).
"""
from pathlib import Path

from mmgbsa.conversion import GromacsPreprocessor
from mmgbsa.trajectory import TrajectoryProcessor


def test_trj_alias_is_created_in_the_run_folder(tmp_path):
    inputs, run = tmp_path / "inputs", tmp_path / "run"
    inputs.mkdir()
    trj = inputs / "dry_md.trj"
    trj.write_text("Cpptraj Generated trajectory\n")
    alias = Path(TrajectoryProcessor.resolve_mdtraj_loadable_path(str(trj), run))
    assert alias == run / "dry_md.mdcrd"
    assert alias.is_symlink() and alias.resolve() == trj.resolve()
    assert sorted(p.name for p in inputs.iterdir()) == ["dry_md.trj"]

    # without a run folder the alias goes to a temporary folder, not next to the input
    other = Path(TrajectoryProcessor.resolve_mdtraj_loadable_path(str(trj)))
    assert other.parent != inputs and other.suffix == ".mdcrd"
    assert sorted(p.name for p in inputs.iterdir()) == ["dry_md.trj"]
    # any other trajectory is returned unchanged
    assert TrajectoryProcessor.resolve_mdtraj_loadable_path("a.xtc", run) == "a.xtc"


def test_gromacs_topology_is_sanitized_into_the_run_folder(tmp_path):
    inputs, run = tmp_path / "inputs", tmp_path / "run"
    inputs.mkdir()
    run.mkdir()
    (inputs / "ligand.itp").write_text("; ligand\n")
    top = inputs / "topol.top"
    original = '  #include "ligand.itp"\n#include "amber99sb.ff/forcefield.itp"\n[ system ]\nx\n'
    top.write_text(original)

    readable = GromacsPreprocessor.sanitize_topology(str(top), run)
    assert Path(readable).parent == run
    text = Path(readable).read_text().splitlines()
    assert text[0] == f'#include "{(inputs / "ligand.itp").resolve()}"'
    assert text[1] == '#include "amber99sb.ff/forcefield.itp"'      # not next to the input: unchanged
    assert top.read_text() == original                               # the input is not touched
    assert sorted(p.name for p in inputs.iterdir()) == ["ligand.itp", "topol.top"]

    clean = inputs / "clean.top"
    clean.write_text('#include "ligand.itp"\n')
    assert GromacsPreprocessor.sanitize_topology(str(clean), run) == clean
