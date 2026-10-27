#!/usr/bin/env python3
"""
OpenMM System input (.xml with a .pdb of the same name).

Before 0.0.20 a serialised System given as `complex_pdb` failed in input
validation. It is now converted with ParmEd to an Amber topology and runs in
Native Mode. The System carries no GB radii; the radii set of the GB model is
assigned, as for GROMACS input. A force ParmEd cannot express in an Amber
topology (CHARMM's pair-specific Lennard-Jones terms as a CustomNonbondedForce)
is refused instead of being dropped.
"""
import os

import numpy as np
import openmm
import parmed as pmd
import pytest
from openmm import app, unit

from mmgbsa.conversion import OpenMMSystemPreprocessor

HERE = os.path.dirname(os.path.abspath(__file__))
RNA = os.path.join(HERE, "..", "configs", "rna_protein_test", "raw")
PRMTOP = os.path.join(RNA, "3dd2.prmtop")
DCD = os.path.join(RNA, "3dd2_suMD1.dcd")


def _energy(prmtop_path, positions):
    prmtop = app.AmberPrmtopFile(prmtop_path)
    system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, implicitSolvent=app.OBC2)
    context = openmm.Context(system, openmm.VerletIntegrator(0.001),
                             openmm.Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)


@pytest.mark.skipif(not (os.path.exists(PRMTOP) and os.path.exists(DCD)),
                    reason="3dd2 data not present (test/configs/rna_protein_test/fetch_data.sh)")
def test_converted_system_has_the_energy_of_the_original(tmp_path):
    import mdtraj as md
    frame = md.load_frame(DCD, 0, top=PRMTOP)
    positions = frame.xyz[0] * unit.nanometer
    prmtop = app.AmberPrmtopFile(PRMTOP)
    system = prmtop.createSystem(nonbondedMethod=app.NoCutoff)
    xml = tmp_path / "complex.xml"
    xml.write_text(openmm.XmlSerializer.serialize(system))
    with open(tmp_path / "complex.pdb", "w") as fh:
        app.PDBFile.writeFile(prmtop.topology, positions, fh)

    converted = OpenMMSystemPreprocessor.convert_to_amber(
        xml, OpenMMSystemPreprocessor.find_pdb(xml), tmp_path / "out", gb_model="OBC2")

    # the reference: the original topology with the same radii set
    from parmed.tools import changeRadii
    original = pmd.load_file(PRMTOP)
    changeRadii(original, "mbondi2").execute()
    reference = str(tmp_path / "reference.prmtop")
    original.save(reference)

    conv = pmd.load_file(converted)
    np.testing.assert_allclose([a.charge for a in conv.atoms], [a.charge for a in original.atoms], atol=1e-6)
    np.testing.assert_allclose([a.solvent_radius for a in conv.atoms],
                               [a.solvent_radius for a in original.atoms], atol=1e-9)
    # full potential energy with OBC2, so bonded, nonbonded and GB terms all enter
    assert _energy(converted, positions) == pytest.approx(_energy(reference, positions), abs=1e-3)


def _two_atom_input(tmp_path, force):
    system = openmm.System()
    for _ in range(2):
        system.addParticle(12.0)
    system.addForce(force)
    topology = app.Topology()
    residue = topology.addResidue("LIG", topology.addChain())
    for name in ("C1", "C2"):
        topology.addAtom(name, app.element.carbon, residue)
    xml = tmp_path / "pair.xml"
    xml.write_text(openmm.XmlSerializer.serialize(system))
    with open(tmp_path / "pair.pdb", "w") as fh:
        app.PDBFile.writeFile(topology, [openmm.Vec3(0, 0, 0), openmm.Vec3(0.4, 0, 0)] * unit.nanometer, fh)
    return xml


def test_custom_nonbonded_force_is_refused(tmp_path):
    force = openmm.CustomNonbondedForce("4*eps*((0.3/r)^12-(0.3/r)^6); eps=0.2")
    force.addParticle([])
    force.addParticle([])
    xml = _two_atom_input(tmp_path, force)
    with pytest.raises(ValueError, match="CustomNonbondedForce"):
        OpenMMSystemPreprocessor.convert_to_amber(xml, OpenMMSystemPreprocessor.find_pdb(xml), tmp_path)


def test_missing_pdb_is_reported(tmp_path):
    xml = tmp_path / "lonely.xml"
    xml.write_text(openmm.XmlSerializer.serialize(openmm.System()))
    with pytest.raises(FileNotFoundError, match="lonely.pdb"):
        OpenMMSystemPreprocessor.find_pdb(xml)
