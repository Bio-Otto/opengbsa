"""Regression test: Coordinate Mode GB radii follow the selected radii set.

Coordinate Mode (a PDB or GROMACS coordinates without a topology file) builds
the system through an OpenMM ForceField, so the topology reaching the GB force
is an openmm.app.Topology whose atoms carry no radii. Up to 0.0.12 every such
atom got one radius per element (every hydrogen 1.2 A) whatever gb_model was
selected, while mbondi2 gives N-bonded hydrogens 1.3 A and mbondi changes more.
On 1GCQ this moved the binding enthalpy by up to 1.5 kcal/mol depending on the
model, and on the IRF11 protein-DNA complex by 22 kcal/mol, silently. The atoms
now get the gb_model's set (or gb_radii_set), assigned by ParmEd's changeRadii
as tleap does; with that, OpenGBSA reproduces Amber on the same frames.
"""
import pytest

pytest.importorskip("openff.toolkit")

import parmed  # noqa: E402
from openff.toolkit.utils import get_data_file_path  # noqa: E402
from openmm import app  # noqa: E402
from parmed.tools import changeRadii  # noqa: E402

from mmgbsa.mmgbsa_core import GBSAForceManager, StructureManager  # noqa: E402

OFFSET_NM = {"HCT": 0.009, "OBC1": 0.009, "OBC2": 0.009, "GBn": 0.009, "GBn2": 0.0195141}


@pytest.fixture(scope="module")
def peptide():
    pdb = app.PDBFile(get_data_file_path("proteins/ace-ala-nh2.pdb"))
    system = app.ForceField("amber14-all.xml").createSystem(
        pdb.topology, nonbondedMethod=app.NoCutoff, constraints=None)
    return pdb.topology, system


def _built_radii(topology, system, gb_model, gb_radii_set="auto"):
    manager = GBSAForceManager(gb_model=gb_model, salt_concentration=0.15)
    manager.gb_radii_set = gb_radii_set
    force = manager._create_fallback_obc_force(system, topology)
    return [round((force.getParticleParameters(i)[1] + OFFSET_NM[gb_model]) * 10, 4)
            for i in range(force.getNumParticles())]


def _parmed_radii(topology, radii_set):
    struct = parmed.openmm.load_topology(topology)
    changeRadii(struct, radii_set).execute()
    return [round(a.solvent_radius, 4) for a in struct.atoms]


@pytest.mark.parametrize("gb_model", ["HCT", "OBC1", "OBC2", "GBn", "GBn2"])
def test_openmm_topology_gets_model_radii_set(peptide, gb_model):
    topology, system = peptide
    expected = _parmed_radii(topology, StructureManager._GB_MODEL_RADII_SET[gb_model])
    assert _built_radii(topology, system, gb_model) == expected


def test_n_bonded_hydrogens_are_not_element_defaults(peptide):
    topology, system = peptide
    radii = _built_radii(topology, system, "OBC2")
    hn = [radii[a.index] for b in topology.bonds() for a, other in ((b[0], b[1]), (b[1], b[0]))
          if a.element.symbol == "H" and other.element.symbol == "N"]
    assert hn and set(hn) == {1.3}          # mbondi2; the old element fallback gave 1.2


def test_explicit_radii_set_overrides(peptide):
    topology, system = peptide
    assert _built_radii(topology, system, "OBC2", gb_radii_set="mbondi") == _parmed_radii(topology, "mbondi")


def test_structure_with_radii_is_left_alone(peptide):
    topology, _ = peptide
    struct = parmed.openmm.load_topology(topology)
    changeRadii(struct, "mbondi").execute()
    manager = GBSAForceManager(gb_model="OBC2")
    assert manager._assigned_radii(struct, list(struct.atoms)) == (None, None)
