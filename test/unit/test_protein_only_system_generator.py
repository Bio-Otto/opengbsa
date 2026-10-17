"""Regression test: the protein-only SystemGenerator skips the small-molecule
force field without changing the System it builds.

Protein-only (PPI) systems used `SystemGenerator(forcefields=...)`, which loads
its default OpenFF small-molecule force field on every call even though no small
molecule is ever added. Passing `small_molecule_forcefield=None` avoids that, but
openmmforcefields' `__init__` then calls `add_molecules(None)`, which raises.
`_ProteinOnlySystemGenerator` skips that call. These tests check that it
constructs with no template generator and that, for a capped peptide, it builds
a System byte-identical (serialized XML) to the default SystemGenerator's.
"""
import pytest

pytest.importorskip("openmmforcefields")
pytest.importorskip("openff.toolkit")

from openff.toolkit.utils import get_data_file_path  # noqa: E402
from openmm import XmlSerializer, app  # noqa: E402

from mmgbsa.mmgbsa_core import SystemGenerator, _ProteinOnlySystemGenerator  # noqa: E402

KWARGS = dict(
    forcefields=["amber/protein.ff14SB.xml"],
    forcefield_kwargs={"constraints": None, "rigidWater": False},
    nonperiodic_forcefield_kwargs={"nonbondedMethod": app.NoCutoff},
)


def test_constructs_without_template_generator():
    generator = _ProteinOnlySystemGenerator(small_molecule_forcefield=None, **KWARGS)
    assert generator.template_generator is None


@pytest.mark.parametrize("pdb", ["proteins/ace-ala-nh2.pdb", "proteins/MainChain_TRP.pdb"])
def test_same_system_as_default_generator(pdb):
    topology = app.PDBFile(get_data_file_path(pdb)).topology
    try:
        default_generator = SystemGenerator(**KWARGS)
    except ImportError as e:
        # The reference loads OpenFF's small-molecule stack, which imports
        # pkg_resources (absent with recent setuptools, e.g. in CI). The
        # protein-only generator never imports it, which is part of the point.
        pytest.skip(f"default SystemGenerator cannot load its OpenFF force field here: {e}")
    default = default_generator.create_system(topology)
    protein_only = _ProteinOnlySystemGenerator(small_molecule_forcefield=None, **KWARGS).create_system(topology)
    assert XmlSerializer.serialize(protein_only) == XmlSerializer.serialize(default)
