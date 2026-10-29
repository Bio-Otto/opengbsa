"""The configured dielectrics and salt screening reach the GB force of every model.

Without salt the GB energy is proportional to (1/eps_in - 1/eps_out), because the
Born radii do not depend on the dielectrics. A solute dielectric of 2 must
therefore scale the GB energy by (1/2 - 1/78.5) / (1 - 1/78.5) for every model, and,
as sander's intdiel does, halve the Coulomb energy (ParmEd used to build HCT, GBn
and GBn2 with its default solute dielectric of 1, and no model scaled Coulomb). The
screening parameter must follow the configured temperature (Amber's constant with
amber_compatible), and the force must carry it.
"""
import copy
import math
from pathlib import Path

import openmm
import parmed
import pytest
from openmm import unit
from parmed.tools import changeRadii

from mmgbsa.mmgbsa_core import GBSAForceManager, StructureManager

DATA = Path(__file__).resolve().parents[1] / "data" / "nmode_baicalein"
MODELS = ["HCT", "OBC1", "OBC2", "GBn", "GBn2"]
RADII = {"HCT": "mbondi", "OBC1": "mbondi2", "OBC2": "mbondi2", "GBn": "bondi", "GBn2": "mbondi3"}


@pytest.fixture(scope="module")
def ligand():
    return parmed.load_file(str(DATA / "ligand.prmtop"), xyz=str(DATA / "ligand_frame29701.pdb"))


def energies(struct, gb_model, solute_dielectric=1.0, salt=0.0, temperature=300.0, amber_salt=False):
    """(Coulomb, GB) energies in kcal/mol of the refined system."""
    struct = copy.copy(struct)
    changeRadii(struct, RADII[gb_model]).execute()
    manager = GBSAForceManager(gb_model=gb_model, salt_concentration=salt, solute_dielectric=solute_dielectric,
                               solvent_dielectric=78.5, temperature=temperature)
    manager.amber_salt = amber_salt
    manager.gb_radii_set = None
    system = StructureManager.create_openmm_system(
        struct, implicitSolvent=manager.current_app_model, implicitSolventKappa=manager.debye_kappa(),
        solventDielectric=78.5)
    system = manager.refine_gbsa_forces(system, struct)
    for f in system.getForces():         # ParmEd's groups are reassigned in the run as well
        f.setForceGroup(0 if isinstance(f, openmm.NonbondedForce) else 2 if isinstance(f, openmm.CustomGBForce) else 5)
    ctx = openmm.Context(system, openmm.VerletIntegrator(1e-3), openmm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(struct.positions)
    e = lambda g: ctx.getState(getEnergy=True, groups={g}).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    return e(0), e(2)


def gb_energy(*args, **kwargs):
    return energies(*args, **kwargs)[1]


def test_debye_kappa():
    m = GBSAForceManager(salt_concentration=0.15, solvent_dielectric=78.5, temperature=300.0)
    assert m.debye_kappa() == pytest.approx(7.3 * 50.33355 * math.sqrt(0.15 / (78.5 * 300.0)))
    m.temperature = 310.0 * unit.kelvin
    assert m.debye_kappa() == pytest.approx(7.3 * 50.33355 * math.sqrt(0.15 / (78.5 * 310.0)))
    m.amber_salt = True                       # sander: no temperature or dielectric in kappa
    assert m.debye_kappa() == pytest.approx(7.3 * math.sqrt(0.10806 * 0.15))
    m.solvent_dielectric = 80.0
    assert m.debye_kappa() == pytest.approx(7.3 * math.sqrt(0.10806 * 0.15))
    assert GBSAForceManager(salt_concentration=0.0).debye_kappa() == 0.0


@pytest.mark.parametrize("gb_model", MODELS)
def test_solute_dielectric_reaches_every_model(ligand, gb_model):
    coul1, gb1 = energies(ligand, gb_model, solute_dielectric=1.0)
    coul2, gb2 = energies(ligand, gb_model, solute_dielectric=2.0)
    assert gb2 / gb1 == pytest.approx((0.5 - 1 / 78.5) / (1 - 1 / 78.5), rel=1e-9)
    assert coul2 / coul1 == pytest.approx(0.5, rel=1e-9)


@pytest.mark.parametrize("gb_model", MODELS)
def test_salt_screening_follows_temperature(ligand, gb_model):
    e300 = gb_energy(ligand, gb_model, salt=0.15, temperature=300.0)
    e350 = gb_energy(ligand, gb_model, salt=0.15, temperature=350.0)
    e0 = gb_energy(ligand, gb_model, salt=0.0)
    assert e300 != pytest.approx(e0, abs=1e-6)
    assert e350 != pytest.approx(e300, abs=1e-6)
    # amber_compatible: the screening of sander, independent of the temperature
    a300 = gb_energy(ligand, gb_model, salt=0.15, temperature=300.0, amber_salt=True)
    a350 = gb_energy(ligand, gb_model, salt=0.15, temperature=350.0, amber_salt=True)
    assert a350 == pytest.approx(a300, abs=1e-9)
