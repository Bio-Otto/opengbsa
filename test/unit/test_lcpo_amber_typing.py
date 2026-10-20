"""LCPO parameter assignment as Amber MMPBSA.py does it (lcpo_typing='amber').

Amber's gbsa=1 code chooses LCPO parameters by element and then by atom-type
name, and gives any other element sp2-carbon parameters; OpenGBSA's default
assigns them by element and bonding. The reference value below is sander's
ESURF (igb=5, gbsa=1, surften=0.0072) for the GAFF2 ligand of the 6T1H_1151
fixture at the coordinates in complex.inpcrd; its chlorine (type 'cl') is the
atom on which the two rules differ.
"""
import pathlib

import openmm
import parmed
import pytest
from openmm import unit

from mmgbsa.lcpo_amber import amber_lcpo_key, type_names_from_structure
from mmgbsa.mmgbsa_core import GBSACalculator, GBSAForceManager

DATA = pathlib.Path(__file__).resolve().parent.parent / "data" / "6t1h_1151_comp"
SANDER_ESURF = 5.3744          # kcal/mol, sander single point on this ligand


@pytest.mark.parametrize("element, type_name, n_heavy, n_total, key", [
    (17, "cl", 1, 1, "C_sp2_2"),        # GAFF chlorine: "carbon SA parms for unknown atom"
    (8, "OC", 1, 1, "O_sp3_1"),         # CHARMM carboxylate O: not 'O' or 'O2'
    (8, "O2", 1, 1, "O_carboxylate"),
    (8, "o", 1, 1, "O_sp2_1"),          # names are compared uppercased
    (7, "NH3", 1, 4, "N_sp2_1"),        # CHARMM lysine N: only 'N3*' is sp3
    (7, "N3", 1, 4, "N_sp3_1"),
    (16, "S", 1, 2, "S_2"),             # CHARMM cysteine S: only 'SH' gets S_1
    (16, "SH", 1, 2, "S_1"),
    (6, "CT", 2, 4, "C_sp3_2"),
    (6, "C", 3, 3, "C_sp2_3"),
    (1, "HC", 1, 1, "none"),
    (30, "Zn", 0, 0, "none"),
    (12, "MG", 0, 0, "Mg"),
    (9, "F", 1, 1, "F"),
])
def test_amber_rule(element, type_name, n_heavy, n_total, key):
    assert amber_lcpo_key(element, type_name, n_heavy, n_total) == key


def _ligand():
    lig = parmed.load_file(str(DATA / "ligand.prmtop"))
    cx = parmed.load_file(str(DATA / "complex.prmtop"), str(DATA / "complex.inpcrd"))
    lig.coordinates = cx[f":{len(cx.residues)}"].coordinates
    return lig


def _esurf(lig, typing):
    manager = GBSAForceManager(sa_model="LCPO")
    manager.lcpo_typing = typing
    system = openmm.System()
    for _ in lig.atoms:
        system.addParticle(1.0)
    system.addForce(manager._setup_lcpo_force(system, lig.topology, type_names_from_structure(lig)))
    ctx = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(lig.positions)
    return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)


def test_amber_typing_reproduces_sander():
    lig = _ligand()
    assert any(a.type == "cl" for a in lig.atoms)
    assert _esurf(lig, "amber") == pytest.approx(SANDER_ESURF, abs=1e-3)
    # the element-based rule gives chlorine its own parameters, so it differs
    assert abs(_esurf(lig, "element") - SANDER_ESURF) > 0.1


def test_preset_sets_amber_typing():
    assert GBSACalculator().lcpo_typing == "element"
    calc = GBSACalculator(amber_compatible=True)
    assert calc.lcpo_typing == "amber" and calc.gbsa_manager.lcpo_typing == "amber"
    assert GBSACalculator(amber_compatible=True, lcpo_typing="element").lcpo_typing == "element"
