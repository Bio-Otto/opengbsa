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
def test_load_gromacs_tpr_branch_is_currently_broken():
    """KNOWN BUG, not a regression introduced by this test: `_load_gromacs`'s
    .tpr branch (topology.py:113-118) calls `parmed.load_file(path)`
    directly on the .tpr, but ParmEd's own format registry cannot identify
    raw GROMACS .tpr files (confirmed: raises `FormatNotFound`, not an
    ImportError or a missing-optional-dependency error -- ParmEd simply has
    no .tpr parser). This is a genuinely separate code path from the one
    `mmgbsa_core.py`'s real native-TPR analysis uses (`mmgbsa/tpr_loader.py`'s
    `load_tpr_as_parmed`, built on the third-party `TprParser` library) --
    `TopologyLoader._load_gromacs` appears to be dead/unreachable from the
    actual CLI/runner execution path (consistent with `mmgbsa/core/`'s
    broader "inert scaffolding from an in-progress refactor" state; see
    `mmgbsa/core/analysis.py`'s own docstring). This test documents the
    current behavior so a future fix (either wiring in real .tpr support
    here, or removing this dead branch) has a clear before/after signal --
    it is NOT a statement that this behavior is desired."""
    import parmed
    from mmgbsa.topology import TopologyLoader
    from mmgbsa.inputs import EngineMode

    with pytest.raises(parmed.exceptions.FormatNotFound):
        TopologyLoader.load_system(GROMACS_TPR, EngineMode.GROMACS)
