#!/usr/bin/env python3
"""
PyMOL hot-spot script of the per-residue decomposition.

Up to 0.0.19 a PyMOL script was written only by the in-calculator "fast"
decomposition, and with an empty chain identifier its selections were invalid
("chain  and ..."). Every per-residue decomposition now writes hotspots.pdb and
view_binding.pml: the receptor atoms carry their residue's contribution as
B-factor, receptor and ligand are told apart by segment identifier, and the
residues keep the numbers and chains of the decomposition tables.
"""
import numpy as np
from openmm import app

from mmgbsa.pymol_viz import write_hotspot_pymol


def _complex():
    topology = app.Topology()
    receptor = topology.addChain()
    for name in ("PHE", "ALA"):
        residue = topology.addResidue(name, receptor)
        for atom in ("N", "CA", "C"):
            topology.addAtom(atom, app.element.get_by_symbol(atom[0]), residue)
    ligand = topology.addResidue("LIG", topology.addChain())
    topology.addAtom("C1", app.element.carbon, ligand)
    return topology


def test_bfactors_segments_and_numbers(tmp_path):
    topology = _complex()
    positions = np.arange(21, dtype=float).reshape(7, 3) * 0.1
    residue_atoms = {"PHE_111_0": [0, 1, 2], "ALA_126_0": [3, 4, 5]}
    values = {"PHE_111_0": -2.72, "ALA_126_0": 0.35}

    script = write_hotspot_pymol(topology, positions, values, residue_atoms, [6], tmp_path,
                                 energy_threshold=-0.5)

    atoms = [l for l in (tmp_path / "hotspots.pdb").read_text().splitlines()
             if l.startswith(("ATOM", "HETATM"))]
    assert len(atoms) == 7
    assert [float(l[60:66]) for l in atoms] == [-2.72] * 3 + [0.35] * 3 + [0.0]
    assert [l[72:76].strip() for l in atoms] == ["REC"] * 6 + ["LIG"]
    assert {(l[17:20], int(l[22:26]), l[21]) for l in atoms[:6]} == {("PHE", 111, "A"), ("ALA", 126, "A")}

    text = script.read_text()
    assert "load hotspots.pdb, complex" in text
    assert "segi REC and b < -0.5" in text
    assert "chain  " not in text
