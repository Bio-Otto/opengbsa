"""Regression test: an Amber .prmtop keeps its own GB radii by default.

0.0.8-0.0.10 assigned the radii set matching gb_model to every structure,
including Amber topologies built with another set. Amber MMPBSA.py always uses
a topology's own radii, so results for, e.g., an mbondi topology run with OBC2
(whose usual set is mbondi2) no longer matched Amber: 3CLpro and 3dd2 moved by
up to 5 kcal/mol. With gb_radii_set='auto' (default) the topology's radii are
kept; an explicit gb_radii_set still overrides; non-Amber inputs still get the
gb_model set; GBn/GBn2 fall back to their own set when the topology's radii are
outside OpenMM's 1-2 A neck range.

mbondi and mbondi2 differ on hydrogens bonded to carbon (1.3 vs 1.2 A), which is
what the tests look at.
"""
import pytest

pytest.importorskip("openff.toolkit")

import parmed  # noqa: E402
from openff.toolkit.utils import get_data_file_path  # noqa: E402
from openmm import app  # noqa: E402
from parmed.tools import changeRadii  # noqa: E402

from mmgbsa.mmgbsa_core import StructureManager  # noqa: E402


@pytest.fixture(scope="module")
def mbondi_peptide(tmp_path_factory):
    d = tmp_path_factory.mktemp("prmtop")
    pdb = app.PDBFile(get_data_file_path("proteins/ace-ala-nh2.pdb"))
    system = app.ForceField("amber14-all.xml").createSystem(
        pdb.topology, nonbondedMethod=app.NoCutoff, constraints=None)
    amber = parmed.amber.AmberParm.from_structure(
        parmed.openmm.load_topology(pdb.topology, system, xyz=pdb.positions))
    changeRadii(amber, "mbondi").execute()
    amber.write_parm(str(d / "pep.prmtop"))
    with open(d / "pep.pdb", "w") as fh:
        app.PDBFile.writeFile(pdb.topology, pdb.positions, fh)
    return str(d / "pep.pdb"), str(d / "pep.prmtop")


def _hc_radii(struct):
    return {round(a.solvent_radius, 2) for a in struct.atoms
            if a.element == 1 and any(b.element == 6 for b in a.bond_partners)}


def test_auto_keeps_topology_radii(mbondi_peptide):
    st = StructureManager.load_complex(*mbondi_peptide, gb_model="OBC2")
    assert _hc_radii(st) == {1.3}          # mbondi kept, not mbondi2's 1.2
    assert st.opengbsa_radii_set is None   # receptor/ligand follow: keep


def test_explicit_set_overrides(mbondi_peptide):
    st = StructureManager.load_complex(*mbondi_peptide, gb_model="OBC2", gb_radii_set="mbondi2")
    assert _hc_radii(st) == {1.2}
    assert st.opengbsa_radii_set == "mbondi2"


def test_non_amber_input_gets_model_set(mbondi_peptide):
    st = parmed.load_file(mbondi_peptide[1])
    assert StructureManager.gb_radii_decision(st, "GBn2", amber_native=False) == "mbondi3"


def test_gbn_falls_back_outside_neck_range(mbondi_peptide):
    st = parmed.load_file(mbondi_peptide[1])
    st.atoms[0].solvent_radius = 0.8      # e.g. an mbondi hydroxyl H
    assert StructureManager.gb_radii_decision(st, "GBn", amber_native=True) == "bondi"
    assert StructureManager.gb_radii_decision(st, "OBC2", amber_native=True) is None


def test_converters_use_the_radii_sets_of_the_calculator():
    """GROMACS .top and OpenMM .xml input get the radii set of each GB model that every
    other input gets (GBn: bondi; the converters used mbondi before 0.0.23)."""
    from mmgbsa.conversion import GromacsPreprocessor
    from mmgbsa.mmgbsa_core import StructureManager
    assert GromacsPreprocessor._GB_MODEL_RADII_SET == StructureManager._GB_MODEL_RADII_SET
    assert StructureManager._GB_MODEL_RADII_SET["GBn"] == "bondi"
