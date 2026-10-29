"""Regression test: CHARMM NBFIX pair corrections reach the energy.

ParmEd 4.3.1's createSystem builds a Lennard-Jones lookup table for structures
with NBFIX and merges an atom type into an earlier type with the same
rmin/epsilon whenever the earlier one has no NBFIX, even if the later one has.
CHARMM36m's arginine NC2 then falls into backbone NH1's class and carboxylate
OC into carbonyl O's, and the NC2-OC NBFIX vanishes. Up to OpenGBSA 0.0.13 every
CHARMM .psf and GROMACS CHARMM36 run was affected (CTLA4 -0.86, VRC01 -2.7,
pembrolizumab -1.8 kcal/mol in van der Waals against Amber). The structure below
reproduces the merge: NH1 (no NBFIX) comes first and NC2 has the same LJ
parameters.
"""
import numpy as np
import openmm
import parmed as pmd
import pytest
from openmm import app, unit

from mmgbsa.mmgbsa_core import StructureManager

NBFIX_RMIN, NBFIX_EPS = 3.637, 0.154919          # CHARMM36m NC2-OC
X = {"N": 0.0, "NH1": 10.0, "OD1": 13.2}          # x positions, A


def _structure():
    nh1 = pmd.AtomType("NH1", None, 14.007, atomic_number=7)
    nc2 = pmd.AtomType("NC2", None, 14.007, atomic_number=7)
    oc = pmd.AtomType("OC", None, 15.999, atomic_number=8)
    nh1.set_lj_params(0.20, 1.85)
    nc2.set_lj_params(0.20, 1.85)                 # same LJ as NH1: ParmEd merges them
    oc.set_lj_params(0.12, 1.70)
    nc2.add_nbfix("OC", NBFIX_RMIN, NBFIX_EPS)
    oc.add_nbfix("NC2", NBFIX_RMIN, NBFIX_EPS)
    s = pmd.Structure()
    for i, (name, t) in enumerate((("N", nh1), ("NH1", nc2), ("OD1", oc))):
        a = pmd.Atom(name=name, type=t.name, charge=0.0, mass=t.mass, atomic_number=t.atomic_number)
        a.atom_type = t
        s.add_atom(a, "RES", i + 1)
    s.coordinates = np.array([[X[a.name], 0.0, 0.0] for a in s.atoms])
    return s


def _lj(eps, rmin, r):
    return eps * ((rmin / r) ** 12 - 2 * (rmin / r) ** 6)


def _vdw_energy(system, coords):
    for f in system.getForces():
        f.setForceGroup(1 if isinstance(f, openmm.CustomNonbondedForce) else 0)
    ctx = openmm.Context(system, openmm.VerletIntegrator(1e-3), openmm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(coords * 0.1)
    return ctx.getState(getEnergy=True, groups=1 << 1).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)


def _lb(e1, r1, e2, r2, r):
    return _lj(np.sqrt(e1 * e2), r1 + r2, r)


EXPECTED = (_lb(0.20, 1.85, 0.20, 1.85, X["NH1"] - X["N"])              # N - NC2
            + _lb(0.20, 1.85, 0.12, 1.70, X["OD1"] - X["N"])            # N - OC
            + _lj(NBFIX_EPS, NBFIX_RMIN, X["OD1"] - X["NH1"]))          # NC2 - OC, NBFIX
WITHOUT_NBFIX = EXPECTED - _lj(NBFIX_EPS, NBFIX_RMIN, 3.2) + _lb(0.20, 1.85, 0.12, 1.70, 3.2)


def test_nbfix_is_applied():
    s = _structure()
    assert s.has_NBFIX()
    system = StructureManager.create_openmm_system(s)
    got = _vdw_energy(system, s.coordinates)
    assert got == pytest.approx(EXPECTED, abs=1e-6)
    assert abs(got - WITHOUT_NBFIX) > 1e-3     # the test can tell the two apart


def test_amber_prmtop_keeps_its_own_table():
    """An Amber prmtop stores the full pair table (NBFIX included); it is used as is."""
    s = _structure()
    amber = pmd.amber.AmberParm.from_structure(s)
    assert amber.has_NBFIX()
    system = StructureManager.create_openmm_system(amber)
    assert _vdw_energy(system, s.coordinates) == pytest.approx(EXPECTED, abs=1e-6)


def test_parmed_alone_drops_it():
    """Documents the upstream behaviour the fix works around (ParmEd 4.3.1)."""
    s = _structure()
    got = _vdw_energy(s.createSystem(nonbondedMethod=app.NoCutoff), s.coordinates)
    fixed = _vdw_energy(StructureManager.create_openmm_system(_structure()),
                        s.coordinates)
    if abs(got - fixed) < 1e-9:
        pytest.skip("this ParmEd version keeps NBFIX on its own")
    assert got != pytest.approx(fixed)
