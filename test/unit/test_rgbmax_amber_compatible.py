"""Born-radius cutoff (rgbmax) and the MMPBSA.py-compatible preset.

Amber (sander, NAB) sums the descreening integral of an atom only over atoms
within rgbmax: an atom j beyond rgbmax + s_j is skipped, and one straddling the
cutoff is integrated up to rgbmax with sander's closed form

    I_ij = (1/8 d) [1 + 2d/(d - s) + (d^2 - 4 rgbmax d - s^2)/rgbmax^2 + 2 ln((d - s)/rgbmax)].

MMPBSA.py fixes rgbmax = 25 A. The reference below implements OBC2 with that
rule directly (numpy, no OpenMM); OpenGBSA's GB force with rgbmax must match it,
and without rgbmax must match the uncut integral.
"""
import numpy as np
import openmm
import pytest
from openmm import unit

from mmgbsa.mmgbsa_core import GBSACalculator, GBSAForceManager

ONE_4PI_EPS0 = 138.935485          # kJ nm / (mol e^2), as in OpenMM's GB expressions
OFFSET = 0.009                     # nm
RGBMAX = 0.6                       # nm (6 A), small so that the test system straddles it
# charge (e), radius (nm), screen; positions (nm) chosen so pairs fall inside,
# across and beyond 0.6 nm
ATOMS = [(-0.40, 0.155, 0.79, (0.00, 0.00, 0.00)),
         (+0.30, 0.120, 0.85, (0.11, 0.00, 0.00)),
         (+0.55, 0.170, 0.72, (0.52, 0.10, 0.00)),
         (-0.65, 0.150, 0.85, (0.66, 0.00, 0.05)),
         (+0.20, 0.180, 0.96, (1.10, 0.20, 0.00)),
         (+0.00, 0.120, 0.85, (1.62, 0.00, 0.10))]


def _hct(d, sj, ri, rgbmax=None):
    """Descreening of atom i (offset radius ri) by atom j (scaled radius sj) at distance d."""
    if rgbmax is not None and d > rgbmax + sj:
        return 0.0
    if rgbmax is not None and d > rgbmax - sj:             # sander's straddling formula
        return 0.125 / d * (1 + 2 * d / (d - sj) + (d * d - 4 * rgbmax * d - sj * sj) / rgbmax ** 2
                            + 2 * np.log((d - sj) / rgbmax))
    if d + sj <= ri:
        return 0.0
    L, U = max(ri, abs(d - sj)), d + sj
    return 0.5 * (1 / L - 1 / U + 0.25 * (d - sj * sj / d) * (1 / U ** 2 - 1 / L ** 2) + 0.5 * np.log(L / U) / d)


def _reference_energy(rgbmax, eps_out=78.5):
    q = np.array([a[0] for a in ATOMS]); x = np.array([a[3] for a in ATOMS])
    orad = np.array([a[1] for a in ATOMS]) - OFFSET
    sr = np.array([a[2] for a in ATOMS]) * orad
    n = len(ATOMS)
    born = np.empty(n)
    for i in range(n):
        I = sum(_hct(np.linalg.norm(x[i] - x[j]), sr[j], orad[i], rgbmax) for j in range(n) if j != i)
        psi = I * orad[i]
        born[i] = 1 / (1 / orad[i] - np.tanh(psi - 0.8 * psi ** 2 + 4.85 * psi ** 3) / (orad[i] + OFFSET))
    pre = -ONE_4PI_EPS0 * (1 - 1 / eps_out)
    e = sum(0.5 * pre * q[i] ** 2 / born[i] for i in range(n))
    for i in range(n):
        for j in range(i + 1, n):
            r2 = np.sum((x[i] - x[j]) ** 2)
            f = np.sqrt(r2 + born[i] * born[j] * np.exp(-r2 / (4 * born[i] * born[j])))
            e += pre * q[i] * q[j] / f
    return e / 4.184


def _openmm_energy(rgbmax):
    from openmm.app.internal import customgbforces
    force = customgbforces.GBSAOBC2Force(solventDielectric=78.5, soluteDielectric=1.0, SA=None, cutoff=None, kappa=0.0)
    for qq, rad, screen, _ in ATOMS:
        force.addParticle([qq, rad, screen])
    force.finalize()
    manager = GBSAForceManager(gb_model="OBC2")
    manager.rgbmax = None if rgbmax is None else rgbmax * 10.0          # A
    manager._apply_rgbmax(force)
    system = openmm.System()
    for _ in ATOMS:
        system.addParticle(1.0)
    system.addForce(force)
    ctx = openmm.Context(system, openmm.VerletIntegrator(1e-3), openmm.Platform.getPlatformByName("Reference"))
    ctx.setPositions([a[3] for a in ATOMS])
    return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)


def test_no_cutoff_matches_the_full_integral():
    assert _openmm_energy(None) == pytest.approx(_reference_energy(None), abs=1e-6)


def test_rgbmax_matches_ambers_rule():
    got, ref, full = _openmm_energy(RGBMAX), _reference_energy(RGBMAX), _reference_energy(None)
    assert got == pytest.approx(ref, abs=1e-6)
    assert abs(ref - full) > 1e-3            # the system does straddle the cutoff


def test_preset_settings():
    calc = GBSACalculator(amber_compatible=True)
    assert (calc.sa_model, calc.rgbmax, calc.gbsa_manager.rgbmax) == ("LCPO", 25.0, 25.0)
    assert calc.gbsa_manager.dielectric_from_preset
    default = GBSACalculator()
    assert (default.sa_model, default.rgbmax, default.solvent_dielectric) == ("LCPO", None, 78.5)
    explicit = GBSACalculator(amber_compatible=True, solvent_dielectric=70.0, rgbmax=20.0)
    assert (explicit.solvent_dielectric, explicit.rgbmax) == (70.0, 20.0)
    assert not explicit.gbsa_manager.dielectric_from_preset


class _Tagged:
    def __init__(self, charmm):
        self.opengbsa_charmm_like = charmm


@pytest.mark.parametrize("charmm, eps", [(True, 80.0), (False, 78.3)])
def test_preset_dielectric_follows_mmpbsa_path(charmm, eps):
    manager = GBSAForceManager()
    manager.dielectric_from_preset = True
    manager.resolve_mmpbsa_dielectric(_Tagged(charmm))
    assert manager.solvent_dielectric == eps
    manager.resolve_mmpbsa_dielectric(_Tagged(not charmm))    # decided once, from the complex
    assert manager.solvent_dielectric == eps
