"""Regression test: one evaluation per force group gives the same per-frame
quantities as the original per-frame calls.

`_frame_energy_terms` (mmgbsa_core) replaced 32 getState() calls per frame, in
which the total energy and the expensive GB / non-polar / nonbonded groups were
evaluated two or three times, with one evaluation per force group whose sum is
the total. This test builds a small system with forces in the code's group
layout (0 electrostatic, 1 GB, 3 van der Waals, 4 non-polar, 10 bonds,
11 angles) and checks every derived quantity against direct getState() calls on
the Reference platform (double precision), where the two must agree exactly up
to summation rounding.
"""
import numpy as np
import openmm
import pytest
from openmm import unit

from mmgbsa.mmgbsa_core import _force_groups, _frame_energy_terms, _group_energies

KCAL = unit.kilocalories_per_mole


def _system(n, seed):
    rng = np.random.default_rng(seed)
    s = openmm.System()
    for _ in range(n):
        s.addParticle(12.0)
    q = rng.uniform(-0.6, 0.6, n)
    nb = openmm.NonbondedForce()                       # electrostatics only
    nb.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
    for qi in q:
        nb.addParticle(qi, 0.3, 0.0)
    vdw = openmm.CustomNonbondedForce("4*epsilon*((sigma/r)^12-(sigma/r)^6); "
                                      "sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")
    vdw.addPerParticleParameter("sigma")
    vdw.addPerParticleParameter("epsilon")
    for _ in range(n):
        vdw.addParticle([0.3, 0.4])
    gb = openmm.GBSAOBCForce()
    gb.setSurfaceAreaEnergy(0.0)
    for qi in q:
        gb.addParticle(qi, 0.15, 0.8)
    sa = openmm.GBSAOBCForce()                         # ACE-style non-polar: zero charges
    for _ in range(n):
        sa.addParticle(0.0, 0.15, 0.8)
    bonds = openmm.HarmonicBondForce()
    angles = openmm.HarmonicAngleForce()
    for i in range(n - 1):
        bonds.addBond(i, i + 1, 0.15, 2.0e5)
    for i in range(n - 2):
        angles.addAngle(i, i + 1, i + 2, 1.9, 300.0)
    for f, g in ((nb, 0), (gb, 1), (vdw, 3), (sa, 4), (bonds, 10), (angles, 11)):
        f.setForceGroup(g)
        s.addForce(f)
    return s


def _context(system):
    return openmm.Context(system, openmm.VerletIntegrator(0.001),
                          openmm.Platform.getPlatformByName("Reference"))


def _e(ctx, groups=None):
    st = ctx.getState(getEnergy=True) if groups is None else ctx.getState(getEnergy=True, groups=groups)
    return st.getPotentialEnergy().value_in_unit(KCAL)


@pytest.fixture(scope="module")
def frame():
    rng = np.random.default_rng(7)
    n_rec, n_lig = 9, 5
    pos = np.cumsum(rng.normal(0, 0.12, (n_rec + n_lig, 3)), axis=0)
    ctxs = []
    for system, p in ((_system(n_rec + n_lig, 1), pos), (_system(n_rec, 1), pos[:n_rec]),
                      (_system(n_lig, 2), pos[n_rec:])):
        c = _context(system)
        c.setPositions(p * unit.nanometer)
        ctxs.append(c)
    return ctxs


def test_group_energies_sum_to_total(frame):
    for ctx in frame:
        groups = _force_groups(ctx.getSystem())
        assert sum(_group_energies(ctx, groups).values()) == pytest.approx(_e(ctx), abs=1e-8)


def test_frame_terms_match_direct_calls(frame):
    c, p, l = frame
    groups = tuple(_force_groups(x.getSystem()) for x in frame)
    t = _frame_energy_terms(c, p, l, groups)
    internal = lambda ctx: _e(ctx, 1 << 10) + _e(ctx, 1 << 11)
    expected_binding = ((_e(c) - internal(c)) - (_e(p) - internal(p)) - (_e(l) - internal(l)))
    assert t["binding_e"] == pytest.approx(expected_binding, abs=1e-8)
    assert t["complex_e"] == pytest.approx(_e(c), abs=1e-8)
    assert t["c_elec"] == pytest.approx(_e(c, 1 << 0), abs=1e-10)
    assert t["c_vdw"] == pytest.approx(_e(c, 1 << 3), abs=1e-10)
    assert t["c_gb"] == pytest.approx(_e(c, {1, 2}), abs=1e-10)
    assert t["c_sa"] == pytest.approx(_e(c, 1 << 4), abs=1e-10)
    assert t["p_sa"] == pytest.approx(_e(p, (1 << 4) | (1 << 16)), abs=1e-10)
    assert t["l_gb"] == pytest.approx(_e(l, {1, 2}), abs=1e-10)
    assert t["e_bond"] == pytest.approx(_e(c, 1 << 10), abs=1e-10)
