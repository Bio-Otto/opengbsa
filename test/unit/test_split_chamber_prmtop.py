"""Receptor and ligand split from a CHARMM (chamber) complex prmtop.

Without receptor_topology/ligand_topology, OpenGBSA builds receptor and ligand by
stripping the complex topology (StructureManager.split_components). Up to 0.0.17
it stripped a generic ParmEd Structure copy, which loses a ChamberParm's CHARMM
1-4 Lennard-Jones parameters and its pair-specific LJ table: on the 6T1H_1151
fixture (CHARMM36 protein, CGenFF-style ligand) the van der Waals binding term
came out +217 kcal/mol too high. The split parts must give the energies of the
fixture's own receptor and ligand prmtops.
"""
import pathlib

import openmm
import parmed
import pytest
from openmm import app, unit

from mmgbsa.mmgbsa_core import StructureManager

DATA = pathlib.Path(__file__).resolve().parent.parent / "data" / "6t1h_1151_comp"
N_RECEPTOR = 3840


def _nonbonded(structure, xyz):
    system = structure.createSystem(nonbondedMethod=app.NoCutoff)
    names = {"NonbondedForce", "CustomNonbondedForce"}
    for i, force in enumerate(system.getForces()):
        force.setForceGroup(1 if type(force).__name__ in names else 0)
    context = openmm.Context(system, openmm.VerletIntegrator(1e-3), openmm.Platform.getPlatformByName("Reference"))
    context.setPositions(xyz * unit.angstrom)
    return context.getState(getEnergy=True, groups={1}).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)


@pytest.mark.parametrize("by_indices", [False, True])
def test_chamber_split_matches_the_component_prmtops(by_indices):
    complex_ = parmed.load_file(str(DATA / "complex.prmtop"), str(DATA / "complex.inpcrd"))
    assert isinstance(complex_, parmed.amber.ChamberParm)
    xyz = complex_.coordinates
    kw = {"ligand_indices": range(N_RECEPTOR, len(complex_.atoms))} if by_indices else {"ligand_resname": "LIG"}
    receptor, ligand = StructureManager.split_components(complex_, **kw)
    ref_receptor = parmed.load_file(str(DATA / "receptor.prmtop"))
    ref_ligand = parmed.load_file(str(DATA / "ligand.prmtop"))
    assert len(receptor.atoms) == N_RECEPTOR and len(ligand.atoms) == len(ref_ligand.atoms)
    assert _nonbonded(receptor, xyz[:N_RECEPTOR]) == pytest.approx(_nonbonded(ref_receptor, xyz[:N_RECEPTOR]), abs=1e-3)
    assert _nonbonded(ligand, xyz[N_RECEPTOR:]) == pytest.approx(_nonbonded(ref_ligand, xyz[N_RECEPTOR:]), abs=1e-3)
