#!/usr/bin/env python3
"""
Re-imaging in a triclinic cell (GROMACS rhombic dodecahedron).

A real 6T1H frame is placed in a rhombic dodecahedron, the ligand and parts
of the protein are moved by lattice vectors (as a wrapped trajectory has
them), and re-imaging must restore the receptor-ligand interaction energy
and the complex energy. Imaging with the box lengths only (what earlier
versions did) does not. Rectangular boxes keep the original arithmetic.
"""
import copy
from pathlib import Path

import mdtraj as md
import numpy as np
import pytest

from mmgbsa.trajectory import reimage_molecules, _lattice_shift

ROOT = Path(__file__).resolve().parents[2]
PRMTOP = ROOT / "test" / "data" / "6t1h_1151_comp" / "complex.prmtop"
XTC = ROOT / "examples" / "quickstart" / "trajectory.xtc"
pytestmark = pytest.mark.skipif(not (PRMTOP.exists() and XTC.exists()),
                                reason="6T1H quickstart fixture not present")


def _dodecahedron(d):
    """Box vectors (rows, nm) of a rhombic dodecahedron with xy-square of side d."""
    return np.array([[d, 0, 0], [0, d, 0], [d / 2, d / 2, d * np.sqrt(2) / 2]])


def _energies(structure, xyz_nm):
    """(complex energy, receptor-ligand interaction energy) in kcal/mol, no cutoff."""
    import openmm
    from openmm import app, unit
    out = []
    for mask in (None, ":LIG", "!:LIG"):
        s = copy.deepcopy(structure)
        keep = np.arange(len(structure.atoms))
        if mask is not None:
            sel = [a.idx for a in s.atoms if (a.residue.name == "LIG") == (mask == "!:LIG")]
            s.strip(mask)
            keep = np.asarray(sel)
        system = s.createSystem(nonbondedMethod=app.NoCutoff)
        context = openmm.Context(system, openmm.VerletIntegrator(0.001),
                                 openmm.Platform.getPlatformByName("Reference"))
        context.setPositions(xyz_nm[keep].astype(np.float64) * unit.nanometer)
        out.append(context.getState(getEnergy=True).getPotentialEnergy()
                   .value_in_unit(unit.kilocalorie_per_mole))
    return out[0], out[0] - out[1] - out[2]


@pytest.fixture(scope="module")
def setup():
    import parmed as pmd
    structure = pmd.load_file(str(PRMTOP))
    frame = reimage_molecules(md.load_frame(str(XTC), 0, top=str(PRMTOP)))   # whole, cubic box
    return structure, frame


def test_lattice_shift_is_the_minimum_image():
    B = _dodecahedron(8.0)
    B_inv = np.linalg.inv(B)
    rng = np.random.default_rng(3)
    for _ in range(200):
        d = rng.uniform(-12, 12, 3)
        remainder = d - _lattice_shift(d, B, B_inv)[0]
        # no lattice vector brings it closer
        n_shifted = [np.linalg.norm(remainder - (np.array(n) - 1) @ B) for n in np.ndindex(3, 3, 3)]
        assert np.linalg.norm(remainder) <= min(n_shifted) + 1e-9


def test_triclinic_frame_is_restored(setup):
    structure, frame = setup
    xyz0 = frame.xyz[0].astype(np.float64)
    e_complex, e_interaction = _energies(structure, xyz0)

    B = _dodecahedron(9.0)
    lig = frame.topology.select("resname LIG")
    block = np.arange(500, 1500)                 # part of the protein chain
    other = np.arange(2600, 2900)
    moved = xyz0.copy()
    moved[lig] += B[0] - B[2]
    moved[block] += B[1] + B[2]
    moved[other] -= B[0]
    scrambled = md.Trajectory(moved[None].astype(np.float32), frame.topology)
    scrambled.unitcell_vectors = B[None].astype(np.float32)
    assert not np.allclose(scrambled.unitcell_angles, 90.0)

    # the scrambled frame really is broken
    e_broken, _ = _energies(structure, scrambled.xyz[0])
    assert abs(e_broken - e_complex) > 1e3

    imaged = reimage_molecules(scrambled)
    e2_complex, e2_interaction = _energies(structure, imaged.xyz[0])
    assert e2_complex == pytest.approx(e_complex, abs=0.05)
    assert e2_interaction == pytest.approx(e_interaction, abs=0.05)
    assert e_interaction < -20.0                    # a bound ligand, not a separated one

    # imaging with the box lengths only (rectangular, as before) does not restore it
    as_rectangle = scrambled[:]
    as_rectangle.unitcell_lengths = scrambled.unitcell_lengths.copy()
    as_rectangle.unitcell_angles = np.full((1, 3), 90.0, dtype=np.float32)
    old = reimage_molecules(as_rectangle)
    e3_complex, e3_interaction = _energies(structure, old.xyz[0])
    assert abs(e3_complex - e_complex) > 10.0 or abs(e3_interaction - e_interaction) > 1.0


def test_rectangular_box_is_unchanged(setup):
    """The rectangular path keeps the original arithmetic: frames of the quickstart
    trajectory image to the same coordinates as with the rule of 0.0.22."""
    _, frame = setup
    traj = md.load(str(XTC), top=str(PRMTOP))[:5]
    new = reimage_molecules(traj)
    old = _reimage_rectangular_reference(traj)
    assert np.array_equal(new.xyz, old.xyz)


def _reimage_rectangular_reference(traj):
    """The 0.0.22 algorithm (box lengths only), for comparison."""
    from collections import deque
    top = traj.topology
    molecules = list(top.find_molecules())
    box = traj.unitcell_lengths
    result = traj[:]
    xyz = result.xyz
    bonds = np.asarray([[b0.index, b1.index] for b0, b1 in top.bonds], dtype=np.int64)
    mol_indices = [np.asarray([a.index for a in m], dtype=np.int64) for m in molecules]
    largest = max(len(m) for m in molecules)
    anchors = [m for m in molecules if len(m) >= max(1, largest // 2)]
    anchor_indices = sorted({a.index for m in anchors for a in m})
    anchor_mols = sorted((np.asarray([a.index for a in m], dtype=np.int64) for m in anchors),
                         key=len, reverse=True)
    anchor_set = set(anchor_indices)
    adj = [[] for _ in range(xyz.shape[1])]
    for b0, b1 in bonds:
        adj[b0].append(b1)
        adj[b1].append(b0)
    roots = [int(mi[0]) for mi in mol_indices]
    for f in range(xyz.shape[0]):
        L = box[f]
        for root, mi in zip(roots, mol_indices):
            visited, queue = {root}, deque([root])
            while queue:
                i = queue.popleft()
                xi = xyz[f, i]
                for j in adj[i]:
                    if j not in visited:
                        visited.add(j)
                        delta = xyz[f, j] - xi
                        xyz[f, j] -= np.round(delta / L) * L
                        queue.append(j)
        placed = [anchor_mols[0]]
        for mi in anchor_mols[1:]:
            ref = xyz[f, np.concatenate(placed)].mean(axis=0)
            com = xyz[f, mi].mean(axis=0)
            xyz[f, mi] -= np.round((com - ref) / L) * L
            placed.append(mi)
        anchor_com = xyz[f, anchor_indices].mean(axis=0)
        for mi in mol_indices:
            if int(mi[0]) in anchor_set:
                continue
            com = xyz[f, mi].mean(axis=0)
            xyz[f, mi] -= np.round((com - anchor_com) / L) * L
    com = xyz[:, anchor_indices, :].mean(axis=1)
    xyz += (box / 2.0 - com)[:, np.newaxis, :]
    return result
