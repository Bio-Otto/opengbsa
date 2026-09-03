#!/usr/bin/env python3
"""
Tests for TopologyLoader._load_amber and TopologyLoader._load_gromacs --
complements test_charmm_loading.py (which covers _load_charmm) so the
three real, currently-used topology-loading branches all have coverage.
_load_openmm_xml and _load_generic are not covered here: the former needs
a hand-serialized System XML fixture, the latter is exercised indirectly
by the ppi_1gcq_test/dna_protein_test Coordinate-Mode validation cases
under test/configs/.

Fixtures: test/data/6t1h_prmtop_test/ (Amber .prmtop-native) and
test/data/7khz_tpr_test/ (GROMACS .tpr-native), both fetched via
`opengbsa --fetch-test-data` (see mmgbsa/test_data_manifest.py).
"""
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))

AMBER_DIR = os.path.join(HERE, "..", "data", "6t1h_prmtop_test")
AMBER_COMPLEX_PRMTOP = os.path.join(AMBER_DIR, "complex.prmtop")

GROMACS_DIR = os.path.join(HERE, "..", "data", "7khz_tpr_test")
GROMACS_TPR = os.path.join(GROMACS_DIR, "md_complex_prod.tpr")


amber_skip = pytest.mark.skipif(
    not os.path.exists(AMBER_COMPLEX_PRMTOP), reason="Amber prmtop test fixture not present"
)
gromacs_skip = pytest.mark.skipif(
    not os.path.exists(GROMACS_TPR), reason="GROMACS tpr test fixture not present"
)


# --- _load_amber --------------------------------------------------------------

@amber_skip
def test_detect_mode_prmtop():
    from mmgbsa.inputs import InputManager, EngineMode
    assert InputManager.detect_mode(AMBER_COMPLEX_PRMTOP) == EngineMode.AMBER


@amber_skip
def test_load_amber_builds_system_and_topology():
    from mmgbsa.topology import TopologyLoader
    from mmgbsa.inputs import EngineMode

    system, topology, positions = TopologyLoader.load_system(
        AMBER_COMPLEX_PRMTOP, EngineMode.AMBER,
    )
    assert system.getNumParticles() > 0
    assert system.getNumParticles() == topology.getNumAtoms()
    # A bare .prmtop carries no coordinates of its own (unlike a
    # .prmtop+.inpcrd pair) -- _load_amber only ever reads the prmtop.
    assert positions is None


@amber_skip
def test_load_amber_respects_nonbonded_method_override():
    import openmm
    import openmm.app as app
    from mmgbsa.topology import TopologyLoader
    from mmgbsa.inputs import EngineMode

    system, _, _ = TopologyLoader.load_system(
        AMBER_COMPLEX_PRMTOP, EngineMode.AMBER, nonbondedMethod=app.NoCutoff,
    )
    forces = system.getForces()
    nb_forces = [f for f in forces if isinstance(f, openmm.NonbondedForce)]
    assert len(nb_forces) == 1
    assert nb_forces[0].getNonbondedMethod() == openmm.NonbondedForce.NoCutoff


# --- _load_gromacs (.tpr branch) -----------------------------------------------

@gromacs_skip
def test_detect_mode_tpr():
    from mmgbsa.inputs import InputManager, EngineMode
    assert InputManager.detect_mode(GROMACS_TPR) == EngineMode.GROMACS


@gromacs_skip
def test_load_gromacs_tpr_branch_raises_not_implemented():
    """`TopologyLoader._load_gromacs` never actually handles .tpr files --
    real native-TPR support lives in `mmgbsa/tpr_loader.py`'s
    `load_tpr_as_parmed` (built on the third-party `TprParser` library) and
    is invoked directly by `mmgbsa_core.py` before a .tpr path would ever
    reach `TopologyLoader.load_system`. This branch now raises
    `NotImplementedError` immediately instead of attempting (and failing) a
    `parmed.load_file()` call, since ParmEd has no .tpr parser."""
    from mmgbsa.topology import TopologyLoader
    from mmgbsa.inputs import EngineMode

    with pytest.raises(NotImplementedError):
        TopologyLoader.load_system(GROMACS_TPR, EngineMode.GROMACS)
