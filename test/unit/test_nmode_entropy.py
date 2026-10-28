"""Normal-mode binding entropy (mmgbsa.nmode_entropy) against Amber's nmode program.

The reference is mmpbsa_py_nabnmode, the program MMPBSA.py runs for its nmode
entropy, on one frame of the baicalein ligand of the 3CLpro dataset
(test/data/nmode_baicalein/README.md): minimized energy -173.52915 kcal/mol,
entropies at 298.15 K of 42.641 (translational), 33.919 (rotational), 49.633
(vibrational) and 126.194 cal/(mol K) in total. OpenGBSA starts from an OBC2
system, as a run would build it, converts it to MMPBSA.py's nmode energy
function (HCT, rgbmax 15 A, no salt, no non-polar term), minimizes it and takes
the Hessian on the double-precision Reference platform.
"""
import re
from pathlib import Path

import numpy as np
import openmm
import pytest
from openmm import app, unit

from mmgbsa import nmode_entropy as nm
from mmgbsa.config import ConfigManager

DATA = Path(__file__).resolve().parents[1] / "data" / "nmode_baicalein"
AMBER_S = {'translational': 42.641, 'rotational': 33.919, 'vibrational': 49.633, 'total': 126.194}
AMBER_ENERGY = -173.52915
COULOMB_RATIO = 1.0000346        # OpenMM's Coulomb constant over Amber's


def species_system(gb_model=app.OBC2, constraints=None, prmtop=DATA / "ligand.prmtop"):
    """The ligand as an OpenGBSA run builds it: GB in group 1 (built-in) or 2 (custom)."""
    system = app.AmberPrmtopFile(str(prmtop)).createSystem(
        nonbondedMethod=app.NoCutoff, implicitSolvent=gb_model, gbsaModel=None,
        constraints=constraints, removeCMMotion=False)
    for force in system.getForces():
        force.setForceGroup(1 if isinstance(force, openmm.GBSAOBCForce)
                            else 2 if isinstance(force, openmm.CustomGBForce) else 0)
    return system


def positions():
    return app.PDBFile(str(DATA / "ligand_frame29701.pdb")).getPositions(asNumpy=True).value_in_unit(unit.nanometer)


def masses(system):
    return [system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(system.getNumParticles())]


def energy(system, pos):
    context = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("Reference"))
    context.setPositions(pos)
    return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)


def test_reproduces_amber_nmode():
    system = nm.amber_nmode_system(species_system())
    result = nm.species_entropy(nm.make_contexts(system, "CPU"), positions(), masses(system), 298.15)
    assert result['converged'] and result['n_imaginary'] == 0
    # same minimum: the energies differ only by the Coulomb constant
    assert result['energy'] == pytest.approx(AMBER_ENERGY, abs=0.02)
    for term, value in AMBER_S.items():
        assert result['entropy'][term] == pytest.approx(value, abs=0.003), term


def test_thermo_constants_reproduce_amber_from_its_frequencies():
    """Amber's own frequencies give its translational and vibrational entropy."""
    out = (DATA / "mmpbsa_py_nabnmode.out").read_text()
    freq = [float(m.group(1)) for m in
            re.finditer(r"^\s+\d+\s+(-?[\d.]+)\s+[\d.]+\s+[\d.]+\s+[\d.]+\s*$", out, re.M)]
    assert len(freq) == 3 * 29 - 6
    system = species_system()
    S = nm.thermo(np.r_[np.zeros(6), freq], masses(system), positions() * 10, 298.15)
    assert S['translational'] == pytest.approx(AMBER_S['translational'], abs=5e-4)
    assert S['vibrational'] == pytest.approx(AMBER_S['vibrational'], abs=5e-4)


@pytest.mark.parametrize("gb_model, radii", [(app.HCT, "mbondi"), (app.OBC1, "mbondi2"), (app.OBC2, "mbondi2"),
                                             (app.GBn, "bondi"), (app.GBn2, "mbondi3")])
def test_amber_nmode_system_from_any_gb_model(tmp_path, gb_model, radii):
    """HCT built from any model's radii and screens equals HCT built from the topology."""
    parmed = pytest.importorskip("parmed")
    from parmed.tools import changeRadii
    structure = parmed.load_file(str(DATA / "ligand.prmtop"))
    changeRadii(structure, radii).execute()
    prmtop = tmp_path / f"ligand_{radii}.prmtop"
    structure.save(str(prmtop))
    pos = positions()
    reference = species_system(app.HCT, prmtop=prmtop)
    hct = [f for f in reference.getForces() if isinstance(f, openmm.CustomGBForce)][0]
    nm.restrict_born_radius_integral(hct, nm.AMBER_NMODE_RGBMAX)
    topology = app.AmberPrmtopFile(str(prmtop)).topology
    converted = nm.amber_nmode_system(species_system(gb_model, prmtop=prmtop), topology)
    assert energy(converted, pos) == pytest.approx(energy(reference, pos), abs=1e-6)


def test_nmode_system_drops_non_polar_term_and_needs_all_bonds():
    system = species_system()
    sa = openmm.CustomExternalForce("0")
    sa.setForceGroup(4)
    system.addForce(sa)
    assert all(f.getForceGroup() != 4 for f in nm.nmode_system(system).getForces())
    with pytest.raises(ValueError, match="constraints"):
        nm.nmode_system(species_system(constraints=app.HBonds))


def test_binding_entropy_sign_units_and_unconverged_frames():
    def species(s_total, converged=True):
        terms = {'translational': 0.0, 'rotational': 0.0, 'vibrational': s_total, 'total': s_total}
        return {'converged': converged, 'rms_gradient': 1e-4 if converged else 0.1, 'n_imaginary': 0,
                'entropy': terms if converged else None}
    frames = [{'frame': 0, 'complex': species(1000.0), 'receptor': species(900.0), 'ligand': species(150.0)},
              {'frame': 1, 'complex': species(1010.0), 'receptor': species(900.0), 'ligand': species(150.0)},
              {'frame': 2, 'complex': species(1000.0, False), 'receptor': species(900.0), 'ligand': species(150.0)}]
    rows, summary = nm.binding_entropy(frames, 300.0)
    # dS = -50 and -40 cal/(mol K): -TdS = +15 and +12 kcal/mol at 300 K
    assert [r['minus_TdS'] for r in rows[:2]] == pytest.approx([15.0, 12.0])
    assert not rows[2]['converged'] and 'minus_TdS' not in rows[2]
    assert summary['n_converged'] == 2 and summary['n_frames'] == 3
    assert summary['value'] == pytest.approx(13.5)
    assert summary['sd'] == pytest.approx(np.std([15.0, 12.0], ddof=1))


def test_select_frames():
    assert nm.select_frames(5, None) == [0, 1, 2, 3, 4]
    assert nm.select_frames(5, 10) == [0, 1, 2, 3, 4]
    assert nm.select_frames(101, 5) == [0, 25, 50, 75, 100]


@pytest.mark.parametrize("setting, ok", [("nmode_frames: 5", True), ("nmode_frames: null", True),
                                         ("nmode_frames: 0", False), ("nmode_drms: -1", False),
                                         ("nmode_maxcyc: 2.5", False), ("nmode_precision: double", True),
                                         ("nmode_precision: single", False)])
def test_config_validation(tmp_path, setting, ok):
    for name in ("complex.prmtop", "traj.xtc"):
        (tmp_path / name).write_text("placeholder\n")
    path = tmp_path / "config.yaml"
    path.write_text(f"input_files:\n  complex_pdb: {tmp_path / 'complex.prmtop'}\n"
                    f"  trajectory: {tmp_path / 'traj.xtc'}\n"
                    f"analysis_settings:\n  entropy_method: nmode\n  {setting}\n")
    manager = ConfigManager(str(path))
    assert manager.validate_config() is ok, manager.get_validation_errors()
