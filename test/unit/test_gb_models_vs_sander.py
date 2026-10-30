"""The GB force of every model against Amber's sander, on a DNA-peptide complex.

test/data/gb_decomp_dna: a six-base-pair DNA strand (DG5 ... DC3) and a peptide
built with tleap (OL15, ff14SB, mbondi3 radii), 289 atoms. The reference values
are the EGB of sander (Amber 24, imin=5, cut=999, rgbmax=999, extdiel=78.5,
gbsa=0) for igb 1, 2, 5, 7 and 8, without salt and with saltcon=0.15. OpenGBSA
builds the GB force from the ParmEd structure, as in Native Mode; with salt it
uses amber_compatible's screening parameter, which is sander's. The GB energies
of the two programs differ by about 0.1-0.2 kcal/mol out of 1,100 here; giving
the terminal residues DG5 and DC3 protein parameters in GBn2 moved it by 131.
"""
from pathlib import Path

import openmm
import parmed
import pytest
from openmm import unit

from mmgbsa.mmgbsa_core import GBSAForceManager, StructureManager

DATA = Path(__file__).resolve().parents[1] / "data" / "gb_decomp_dna"
SANDER = {  # model: (EGB saltcon=0, EGB saltcon=0.15), kcal/mol
    "HCT": (-1069.0019, -1072.4911),
    "OBC1": (-1087.8774, -1091.3646),
    "OBC2": (-1046.6438, -1050.1225),
    "GBn": (-1092.3387, -1095.8347),
    "GBn2": (-1110.9368, -1114.4277),
}


@pytest.fixture(scope="module")
def complex_structure():
    return parmed.load_file(str(DATA / "complex.prmtop"), xyz=str(DATA / "complex.inpcrd"))


def gb_energy(struct, gb_model, salt):
    manager = GBSAForceManager(gb_model=gb_model, salt_concentration=salt, solvent_dielectric=78.5)
    manager.amber_salt = True
    manager.gb_radii_set = None            # keep the topology's mbondi3 radii, as sander
    system = StructureManager.create_openmm_system(
        struct, implicitSolvent=manager.current_app_model, implicitSolventKappa=manager.debye_kappa(),
        solventDielectric=78.5)
    system = manager.refine_gbsa_forces(system, struct)
    for f in system.getForces():
        f.setForceGroup(2 if isinstance(f, openmm.CustomGBForce) else 5)
    ctx = openmm.Context(system, openmm.VerletIntegrator(1e-3), openmm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(struct.positions)
    return ctx.getState(getEnergy=True, groups={2}).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)


@pytest.mark.parametrize("gb_model", list(SANDER))
@pytest.mark.parametrize("salt", [0.0, 0.15])
def test_gb_energy_matches_sander(complex_structure, gb_model, salt):
    reference = SANDER[gb_model][0 if salt == 0.0 else 1]
    assert gb_energy(complex_structure, gb_model, salt) == pytest.approx(reference, abs=0.3)


@pytest.mark.parametrize("gb_model", list(SANDER))
def test_salt_effect_matches_sander(complex_structure, gb_model):
    """The screening of amber_compatible reproduces sander's salt effect closely."""
    sander = SANDER[gb_model][1] - SANDER[gb_model][0]
    ours = gb_energy(complex_structure, gb_model, 0.15) - gb_energy(complex_structure, gb_model, 0.0)
    assert ours == pytest.approx(sander, abs=0.01)
