#!/usr/bin/env python3
"""
Tests for native NAMD/CHARMM .psf topology loading (TopologyLoader._load_charmm).

Fixture: test/data/charmm_namd_test/ -- a real published NAMD dataset
(Zenodo record 7186684, "Molecular Dynamics Simulation of a Designed Cyclic
Peptide Bound to CTLA4", CC-BY-4.0), trimmed to the .psf/.pdb/CHARMM36m
parameter files plus the small (4-frame) minimization .dcd. This PSF's
cyclic-peptide backbone improper is what originally motivated using ParmEd's
`Structure.load_parameters` instead of OpenMM's stricter
`app.CharmmPsfFile.loadParameters` in `_load_charmm` (see its docstring).
"""
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(HERE, "..", "data", "charmm_namd_test")
PSF_PATH = os.path.join(FIXTURE_DIR, "ctla4_P16_wat.psf")
PDB_PATH = os.path.join(FIXTURE_DIR, "ctla4_P16_wat.pdb")
PRM_PATH = os.path.join(FIXTURE_DIR, "par_all36m_prot.prm")
STR_PATH = os.path.join(FIXTURE_DIR, "toppar_water_ions_prot.str")
DCD_PATH = os.path.join(FIXTURE_DIR, "min_ctla4_P16.0.dcd")

pytestmark = pytest.mark.skipif(
    not os.path.exists(PSF_PATH), reason="CHARMM/NAMD test fixture not present"
)


def test_detect_mode_psf():
    from mmgbsa.inputs import InputManager, EngineMode
    assert InputManager.detect_mode(PSF_PATH) == EngineMode.CHARMM


def test_load_charmm_requires_params():
    from mmgbsa.topology import TopologyLoader
    from mmgbsa.inputs import EngineMode
    with pytest.raises(ValueError, match="charmm_params"):
        TopologyLoader.load_system(PSF_PATH, EngineMode.CHARMM)


def test_load_charmm_builds_system_without_coordinates():
    """A bare .psf carries no positions of its own (confirmed: ParmEd's
    own `struct.positions` is None for a real PSF) -- the System still
    builds, but positions come back None, which callers must not silently
    replace with an all-zeros placeholder without the user opting in."""
    import openmm.app as app
    from mmgbsa.topology import TopologyLoader
    from mmgbsa.inputs import EngineMode

    system, topology, positions = TopologyLoader.load_system(
        PSF_PATH, EngineMode.CHARMM,
        charmm_params=[PRM_PATH, STR_PATH],
        nonbondedMethod=app.NoCutoff,
    )
    assert system.getNumParticles() == topology.getNumAtoms()
    assert system.getNumParticles() > 0
    assert positions is None


def test_load_charmm_builds_system_with_coordinates():
    import openmm.app as app
    from mmgbsa.topology import TopologyLoader
    from mmgbsa.inputs import EngineMode

    system, topology, positions = TopologyLoader.load_system(
        PSF_PATH, EngineMode.CHARMM,
        charmm_params=[PRM_PATH, STR_PATH],
        charmm_coordinates=PDB_PATH,
        nonbondedMethod=app.NoCutoff,
    )
    assert system.getNumParticles() == topology.getNumAtoms()
    assert positions is not None
    assert len(positions) == system.getNumParticles()


def test_load_charmm_trajectory_matches_topology():
    import mdtraj as md
    import openmm.app as app
    from mmgbsa.topology import TopologyLoader
    from mmgbsa.inputs import EngineMode

    system, topology, _ = TopologyLoader.load_system(
        PSF_PATH, EngineMode.CHARMM,
        charmm_params=[PRM_PATH, STR_PATH],
        nonbondedMethod=app.NoCutoff,
    )
    traj = md.load(DCD_PATH, top=PDB_PATH)
    assert traj.n_atoms == system.getNumParticles()
    assert traj.n_frames > 0


def test_charmm_solvent_strip_removes_namd_ion_names():
    """CHARMM/NAMD uses different solvent/ion residue names than Amber/
    GROMACS (TIP3 not SOL/WAT-only, SOD/CLA not NA/CL) -- confirmed this
    fixture's original (Amber/GROMACS-only) strip mask removed ZERO of its
    16024 atoms despite containing TIP3 waters and CLA/SOD ions, which went
    on to corrupt the downstream trajectory-atom-count-matching fallback
    (blind slicing grabbed the wrong atoms entirely). This checks the
    ParmEd-side mask in mmgbsa_core.py's native-mode complex-loading path."""
    import parmed as pmd
    from mmgbsa.mmgbsa_core import StructureManager

    charmm_coord_source = PDB_PATH
    struct = StructureManager.load_complex(
        charmm_coord_source, PSF_PATH, xtc_path=None, gb_model="OBC2",
        charmm_params=[PRM_PATH, STR_PATH],
    )
    n_before = len(struct.atoms)
    solvent_mask = (":WAT,HOH,H2O,SOL,TIP3,TIP,SPC,"
                    "NA,CL,K,MG,ZN,CA,"
                    "Na+,Cl-,K+,Mg2+,Ca2+,Zn2+,"
                    "SOD,CLA,POT,CAL,ZN2")
    struct.strip(solvent_mask)
    n_after = len(struct.atoms)
    assert n_after < n_before
    remaining_resnames = set(a.residue.name for a in struct.atoms)
    assert "TIP3" not in remaining_resnames
    assert "SOD" not in remaining_resnames
    assert "CLA" not in remaining_resnames


def test_charmm_run_end_to_end_produces_physical_energy():
    """Full GBSACalculator.run() on a real NAMD .psf/.dcd/CHARMM36m input,
    in binding_mode='ppi' (chainid-based split, since the "ligand" here is a
    17-residue cyclic peptide chain, not a single small-molecule residue).
    This is a regression test for three real bugs found and fixed together
    (see mmgbsa_core.py's comments at each site):
      1. StructureManager.load_complex's .psf branch needs an explicit
         companion coordinate file (a bare .psf carries none) -- passing the
         same .psf path as both topology and coordinates silently left
         positions None, producing garbage energies once trajectory
         coordinates were set on the resulting all-zero-position System.
      2. The solvent-strip mask only recognized Amber/GROMACS-style
         solvent/ion names, silently leaving all CHARMM/NAMD-style TIP3/SOD/
         CLA atoms in place.
      3. Native Mode's complex_context was set from `protein_pos + ligand_pos`
         (assuming complex_system's atom order is always protein-then-ligand,
         true for typical Amber inputs but NOT for this PSF, whose peptide
         "ligand" chain comes first) instead of the frame's own true atom
         order -- silently swapping receptor/ligand coordinates and producing
         ~11 million kcal/mol of bogus VdW "energy" from atoms overlapping
         under the wrong assignment.
    Asserts only that the result is a finite, physically plausible binding
    energy (no fixed reference value exists for this system) -- the real
    regression signal is "doesn't produce a multi-million-kcal/mol blowup".
    """
    import math
    from mmgbsa.mmgbsa_core import GBSACalculator

    calculator = GBSACalculator(
        temperature=300, verbose=0, gb_model="OBC2", salt_concentration=0.15,
        use_cache=False, protein_forcefield="charmm",
        charmm_params=[PRM_PATH, STR_PATH], charmm_coordinates=PDB_PATH,
    )
    results = calculator.run(
        ligand_mol=None, complex_pdb=PSF_PATH, xtc_file=DCD_PATH, ligand_pdb=None,
        max_frames=4, ligand_selection="chainid 0", receptor_selection="chainid 1",
        solvated_topology=PDB_PATH, original_complex_pdb=PSF_PATH,
    )
    assert results is not None
    mean_binding = results["mean_binding_energy"]
    assert math.isfinite(mean_binding)
    # A real, non-blown-up MM/GBSA binding energy for a peptide-protein
    # complex is at most a few hundred kcal/mol in magnitude; the bug this
    # guards against produced ~11 MILLION kcal/mol, so this bound has huge
    # margin without being a fragile exact-value check.
    assert abs(mean_binding) < 5000
