#!/usr/bin/env python3
"""The per-residue decomposition decomposes the run's own energies exactly.

The decomposition reads every parameter from the OpenMM Systems the MM/GBSA
calculation evaluates (mmgbsa.solvation_peratom for the GB and surface-area
forces, mmgbsa.decomposition._NonbondedModel for the nonbonded forces). Up to
0.0.22 the polar term was recomputed with OBC2 and a fixed 25 A Born-radius
cutoff whatever gb_model and rgbmax the run used, and the non-polar term with
Shrake-Rupley SASA whatever sa_model it used.

Checked here, against OpenMM's Reference platform:
- per-atom polar energies sum to the GB force-group energy of the complex,
  receptor and ligand for HCT, OBC1, OBC2, GBn and GBn2, without and with
  rgbmax 25 A, at salt 0 and 0.15 M (and GBn2 with nucleic-acid parameters);
- per-atom non-polar energies sum to the surface-area group energy for LCPO
  (element and amber typing) and ACE;
- amber_compatible selects MMPBSA.py's decomposition convention for the
  per-residue non-polar term (numerical surface, as sander gbsa=2);
- the per-residue van der Waals and electrostatic halves add up to OpenMM's
  interaction energy, also when the run divides the Coulomb term by the
  solute dielectric, and residues plus the ligand's share give the binding
  energy; report_raw_energies and the parallel workers leave them unchanged.

Systems are built as GBSACalculator's Native Mode builds them
(StructureManager.create_openmm_system, GBSAForceManager.refine_gbsa_forces,
then run()'s force groups).

Fixtures (tracked in git): 6T1H with ligand 1151 (test/data/6t1h_1151_comp,
CHARMM36 chamber prmtops) at frame 0 of examples/quickstart/trajectory.xtc,
as a 10 A pocket for the GB matrix and whole for the binding-energy checks;
a DNA strand with a bound pentapeptide (test/data/gb_decomp_dna, tleap).
"""
import os

import numpy as np
import pytest
from openmm import app, openmm, unit

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
T1H = os.path.join(ROOT, "test", "data", "6t1h_1151_comp")
XTC = os.path.join(ROOT, "examples", "quickstart", "trajectory.xtc")
DNA = os.path.join(ROOT, "test", "data", "gb_decomp_dna")

MODELS = ["HCT", "OBC1", "OBC2", "GBn", "GBn2"]
TOL = 1e-3  # kcal/mol

pytestmark = pytest.mark.skipif(
    not (os.path.exists(os.path.join(T1H, "complex.prmtop")) and os.path.exists(XTC)
         and os.path.exists(os.path.join(DNA, "complex.prmtop"))),
    reason="decomposition fixtures not present")


# ---------------------------------------------------------------------------
# Systems built the way the MM/GBSA run builds them
# ---------------------------------------------------------------------------

def _force_groups(system):
    """GBSACalculator.run()'s assign_force_groups."""
    groups = {openmm.HarmonicBondForce: 10, openmm.HarmonicAngleForce: 11,
              openmm.PeriodicTorsionForce: 12, openmm.RBTorsionForce: 13,
              openmm.CMAPTorsionForce: 14, openmm.CustomTorsionForce: 15,
              openmm.NonbondedForce: 0, openmm.CustomGBForce: 2, openmm.CustomNonbondedForce: 3}
    for f in system.getForces():
        if type(f) in groups:
            f.setForceGroup(groups[type(f)])
        elif isinstance(f, openmm.CustomBondForce) and "sigma" in f.getEnergyFunction():
            f.setForceGroup(3)
    return system


def build(struct, model="OBC2", salt=0.15, rgbmax=None, sa="LCPO", typing="element",
          eps_in=1.0, radii=True, own_gb_builder=False):
    """own_gb_builder: replace the GB force by GBSAForceManager's own builder
    (_create_fallback_obc_force) applied to the OpenMM topology, as in
    Coordinate Mode: it gives DNA/RNA residues GBn2's nucleic-acid parameters
    (ParmEd's GBn2 force gives them the protein ones)."""
    from parmed.tools import changeRadii
    from mmgbsa.mmgbsa_core import GBSAForceManager, StructureManager

    manager = GBSAForceManager(gb_model=model, salt_concentration=salt, solute_dielectric=eps_in,
                               solvent_dielectric=78.5, sa_model=sa)
    manager.rgbmax = rgbmax
    manager.lcpo_typing = typing
    if radii:
        changeRadii(struct, StructureManager._GB_MODEL_RADII_SET[model]).execute()
    system = StructureManager.create_openmm_system(
        struct, implicitSolvent=manager.current_app_model,
        implicitSolventSaltConc=salt * unit.molar, nonbondedMethod=app.NoCutoff,
        nonbondedCutoff=None, solventDielectric=78.5)
    system = manager.refine_gbsa_forces(system, struct)
    if own_gb_builder:
        for k in reversed(range(system.getNumForces())):
            if isinstance(system.getForce(k), openmm.CustomGBForce):
                system.removeForce(k)
        gb = manager._create_fallback_obc_force(system, struct.topology)
        manager._apply_rgbmax(gb)
        system.addForce(gb)
    return _force_groups(system)


def reference_context(system, positions_nm):
    ctx = openmm.Context(system, openmm.VerletIntegrator(1e-3), openmm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(positions_nm)
    return ctx


def group_energy(ctx, groups):
    from mmgbsa.solvation_peratom import group_energy as ge
    return ge(ctx, groups)


@pytest.fixture(scope="module")
def t1h():
    """6T1H complex/receptor/ligand prmtops and frame 0 of the quickstart trajectory."""
    import mdtraj as md
    import parmed as pmd

    complex_top = os.path.join(T1H, "complex.prmtop")
    pos = md.load_frame(XTC, 0, top=complex_top).xyz[0].astype(np.float64)
    return {"complex": complex_top, "receptor": os.path.join(T1H, "receptor.prmtop"),
            "ligand": os.path.join(T1H, "ligand.prmtop"), "pos": pos,
            "lig": [a.idx for a in pmd.load_file(complex_top).atoms if a.residue.name == "LIG"]}


@pytest.fixture(scope="module")
def pocket(t1h):
    """The ligand and every residue within 10 A of it, plus the receptor and
    ligand parts. Slicing loses the chamber 1-4 LJ, which the GB and
    surface-area checks below do not involve."""
    import parmed as pmd

    full = pmd.load_file(t1h["complex"])
    full.coordinates = t1h["pos"] * 10.0
    cplx = full[":LIG<:10"]
    is_lig = [a.residue.name == "LIG" for a in cplx.atoms]
    pos = cplx.coordinates / 10.0
    return {"complex": (lambda: cplx.copy(pmd.Structure), pos),
            "receptor": (lambda: cplx[[not x for x in is_lig]], pos[[not x for x in is_lig]]),
            "ligand": (lambda: cplx[is_lig], pos[is_lig])}


@pytest.fixture(scope="module")
def dna():
    import parmed as pmd

    c = pmd.load_file(os.path.join(DNA, "complex.prmtop"), os.path.join(DNA, "complex.inpcrd"))
    pos = c.coordinates / 10.0
    is_lig = np.array([not a.residue.name.startswith("D") for a in c.atoms])
    return {"complex": (lambda: pmd.load_file(os.path.join(DNA, "complex.prmtop")), pos),
            "receptor": (lambda: pmd.load_file(os.path.join(DNA, "receptor.prmtop")), pos[~is_lig]),
            "ligand": (lambda: pmd.load_file(os.path.join(DNA, "ligand.prmtop")), pos[is_lig]),
            "lig": np.nonzero(is_lig)[0]}


def _check_sums(parts, **kw):
    """Per-atom sums vs OpenMM for complex, receptor, ligand; returns the complex GB energy."""
    from mmgbsa.solvation_peratom import POLAR_GROUPS, NONPOLAR_GROUPS, SolvationEnergies

    out = {}
    for name in ("complex", "receptor", "ligand"):
        make, pos = parts[name]
        system = build(make(), **kw)
        polar, nonpolar = SolvationEnergies(system).per_atom(pos)
        ctx = reference_context(system, pos)
        ref_polar, ref_nonpolar = group_energy(ctx, POLAR_GROUPS), group_energy(ctx, NONPOLAR_GROUPS)
        assert polar.sum() == pytest.approx(ref_polar, abs=TOL), f"{name} polar"
        assert nonpolar.sum() == pytest.approx(ref_nonpolar, abs=TOL), f"{name} non-polar"
        assert abs(ref_polar) > 1.0
        out[name] = (ref_polar, ref_nonpolar, system)
    return out


# ---------------------------------------------------------------------------
# Polar term: every GB model, rgbmax, salt
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("salt", [0.0, 0.15])
@pytest.mark.parametrize("model", MODELS)
def test_polar_term_sums_to_openmm_gb_energy(pocket, model, salt):
    plain = _check_sums(pocket, model=model, salt=salt, rgbmax=None, sa="LCPO", typing="amber")
    cut = _check_sums(pocket, model=model, salt=salt, rgbmax=25.0, sa="LCPO", typing="amber")
    # the pocket is wider than 25 A, so the cutoff is exercised
    assert abs(plain["complex"][0] - cut["complex"][0]) > 0.01


@pytest.mark.parametrize("rgbmax", [None, 25.0])
@pytest.mark.parametrize("salt", [0.0, 0.15])
def test_gbn2_nucleic_acid_parameters(dna, salt, rgbmax):
    out = _check_sums(dna, model="GBn2", salt=salt, rgbmax=rgbmax, sa="LCPO", typing="element",
                      own_gb_builder=True)
    gb = next(f for f in out["receptor"][2].getForces() if isinstance(f, openmm.CustomGBForce))
    names = [gb.getPerParticleParameterName(k) for k in range(gb.getNumPerParticleParameters())]
    alpha = np.array([gb.getParticleParameters(i)[names.index("alpha")] for i in range(gb.getNumParticles())])
    from openmm.app.internal.customgbforces import GBSAGBn2Force
    nucleic = {p[0] for p in GBSAGBn2Force._atom_params_nucleic.values()}
    assert np.isin(np.round(alpha, 6), np.round(list(nucleic), 6)).sum() > 50  # DNA atoms use them


# ---------------------------------------------------------------------------
# Non-polar term: LCPO (element / amber typing) and ACE
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("fixture,sa,typing", [("dna", "LCPO", "element"), ("dna", "LCPO", "amber"),
                                                ("dna", "ACE", "element"), ("pocket", "LCPO", "amber"),
                                                ("pocket", "ACE", "element")])
def test_nonpolar_term_sums_to_openmm_sa_energy(request, fixture, sa, typing):
    out = _check_sums(request.getfixturevalue(fixture), model="OBC2", salt=0.15, sa=sa, typing=typing)
    assert out["complex"][1] > 1.0


def test_gbsaobcforce_with_ace(pocket):
    """OpenMM's own GBSAOBCForce (OBC2 + ACE), engulfed-atom term included."""
    from mmgbsa.solvation_peratom import SolvationEnergies

    make, pos = pocket["complex"]
    custom = build(make(), model="OBC2", salt=0.0, sa="LCPO", typing="amber")
    gb = next(f for f in custom.getForces() if isinstance(f, openmm.CustomGBForce))
    force = openmm.GBSAOBCForce()
    force.setSoluteDielectric(1.0)
    force.setSolventDielectric(78.5)
    force.setSurfaceAreaEnergy(2.25936)
    for i in range(gb.getNumParticles()):
        q, orad, srad = gb.getParticleParameters(i)
        force.addParticle(q, orad + 0.009, srad / orad)
    force.setForceGroup(1)
    system = openmm.System()
    for _ in range(gb.getNumParticles()):
        system.addParticle(1.0)
    system.addForce(force)
    polar, _ = SolvationEnergies(system).per_atom(pos)
    assert polar.sum() == pytest.approx(group_energy(reference_context(system, pos), (1,)), abs=TOL)


# ---------------------------------------------------------------------------
# Nonbonded halves, binding-energy sum, serial/parallel/raw consistency
# ---------------------------------------------------------------------------

def _systems(parts, **kw):
    return [build(parts[k][0](), **kw) for k in ("complex", "receptor", "ligand")]


def _residue_map(struct, lig):
    lig = set(int(i) for i in lig)
    out = {}
    for a in struct.atoms:
        if a.idx not in lig:
            out.setdefault(f"{a.residue.name}_{a.residue.number}_{a.residue.chain or 'A'}", []).append(a.idx)
    return out


def _decomposer(systems, residue_map, lig, raw=False, platform="Reference"):
    from mmgbsa.decomposition import _FrameDecomposer
    plat = openmm.Platform.getPlatformByName(platform)
    ctx = {k: openmm.Context(s, openmm.VerletIntegrator(1e-3), plat)
           for k, s in zip(("complex", "receptor", "ligand"), systems)}
    return _FrameDecomposer(*systems, residue_map, lig, report_raw_energies=raw, contexts=ctx)


def _openmm_delta(systems, pos, lig, groups):
    rec = np.setdiff1d(np.arange(len(pos)), lig)
    e = [group_energy(reference_context(s, p), groups) for s, p in zip(systems, (pos, pos[rec], pos[lig]))]
    return e[0] - e[1] - e[2]


def test_residues_and_ligand_share_sum_to_binding_energy(t1h):
    """Whole 6T1H complex with its own receptor and ligand prmtops (CHARMM36
    chamber: NBFIX LJ table, 1-4 LJ CustomBondForce)."""
    import parmed as pmd

    structs = [pmd.load_file(t1h[k]) for k in ("complex", "receptor", "ligand")]
    systems = [build(s, model="OBC2", salt=0.15, sa="LCPO", typing="element") for s in structs]
    residue_map = _residue_map(structs[0], t1h["lig"])
    pos, lig = t1h["pos"], np.array(t1h["lig"])
    result, check = _decomposer(systems, residue_map, lig, platform="CPU").decompose(pos)
    elec = sum(r["electrostatic"] for r in result.values())
    vdw = sum(r["vdw"] for r in result.values())
    assert 2 * elec == pytest.approx(_openmm_delta(systems, pos, lig, (0,)), abs=TOL)
    assert 2 * vdw == pytest.approx(_openmm_delta(systems, pos, lig, (3,)), abs=TOL)
    assert abs(vdw) > 10.0
    # the polar/non-polar sums against the CPU platform's (single-precision) energies
    assert check["max_solvation_deviation"] < 0.01
    assert check["decomposed_binding_energy"] == pytest.approx(check["openmm_binding_energy"], abs=0.01)
    assert sum(r["solvation"] for r in result.values()) == pytest.approx(check["receptor_solvation"])


def test_report_raw_energies_and_parallel_workers_give_identical_values(dna):
    from mmgbsa import decomposition as d

    systems = _systems(dna, model="GBn2", salt=0.15, rgbmax=25.0, sa="LCPO", typing="element")
    pos, lig = dna["complex"][1].astype(np.float32), dna["lig"]
    residue_map = _residue_map(dna["complex"][0](), lig)
    plain, _ = _decomposer(systems, residue_map, lig, raw=False, platform="CPU").decompose(pos)
    raw, _ = _decomposer(systems, residue_map, lig, raw=True, platform="CPU").decompose(pos)
    for res_id, e in plain.items():
        for comp in d.FRAME_COMPONENTS:
            assert raw[res_id][comp] == e[comp]
        # receptor value = complex value - binding contribution
        assert raw[res_id]["receptor_solvation"] == pytest.approx(raw[res_id]["complex_solvation"] - e["solvation"])
    # the complex view counts every pair, intra-residue 1-2/1-3/1-4 exceptions included
    complex_elec = sum(r["complex_electrostatic"] for r in raw.values())
    rec = np.setdiff1d(np.arange(len(pos)), lig)
    receptor_ref = group_energy(reference_context(systems[1], pos[rec].astype(np.float64)), (0,))
    assert complex_elec - sum(r["electrostatic"] for r in raw.values()) == pytest.approx(receptor_ref, abs=TOL)

    # the parallel path: forked workers on the parent's decomposer, OpenMM check in the parent
    dec = _decomposer(systems, residue_map, lig, raw=True, platform="CPU")
    frames = [pos, pos + np.float32(0.001)]
    parallel = d._decompose_parallel(dec, frames, 2)
    for xyz, (result, check) in zip(frames, parallel):
        serial, serial_check = dec.decompose(xyz)
        assert result == serial                   # bit for bit
        assert d._frame_row(result, 7, 1, 1) == d._frame_row(serial, 7, 1, 1)
        assert "_pending" not in check
        assert check["max_solvation_deviation"] < 0.01
        # OpenMM CPU-platform energies vary by ~1e-4 between evaluations
        assert check["openmm_binding_energy"] == pytest.approx(serial_check["openmm_binding_energy"], abs=TOL)


def _scaled_coulomb_systems(dna, eps_in):
    """The run's Systems at solute dielectric eps_in. A build that does not yet
    divide the Coulomb term by eps_in (NonbondedForce charges equal to the GB
    charges) gets that scaling here: charges times 1/sqrt(eps_in), exception
    charge products times 1/eps_in. The GB force keeps the unscaled charges."""
    systems = _systems(dna, model="OBC2", salt=0.15, sa="LCPO", typing="element", eps_in=eps_in)
    for system in systems:
        nb = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
        gb = next(f for f in system.getForces() if isinstance(f, openmm.CustomGBForce))
        i0 = next(i for i in range(nb.getNumParticles())
                  if abs(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)) > 0.05)
        q_nb = nb.getParticleParameters(i0)[0].value_in_unit(unit.elementary_charge)
        if eps_in != 1.0 and q_nb == pytest.approx(gb.getParticleParameters(i0)[0]):
            for i in range(nb.getNumParticles()):
                q, s, e = nb.getParticleParameters(i)
                nb.setParticleParameters(i, q / np.sqrt(eps_in), s, e)
            for k in range(nb.getNumExceptions()):
                i, j, qq, s, e = nb.getExceptionParameters(k)
                nb.setExceptionParameters(k, i, j, qq / eps_in, s, e)
    return systems


def test_electrostatics_follow_the_runs_coulomb_scaling(dna):
    """With solute dielectric 2 the per-residue electrostatics sum to the
    scaled Coulomb interaction (half of the unscaled one); the polar term,
    from the GB force's unscaled charges with eps_in in its prefactor, still
    sums to OpenMM's GB energy."""
    from mmgbsa.solvation_peratom import POLAR_GROUPS

    pos, lig = dna["complex"][1], dna["lig"]
    residue_map = _residue_map(dna["complex"][0](), lig)
    elec = {}
    for eps_in in (1.0, 2.0):
        systems = _scaled_coulomb_systems(dna, eps_in)
        result, check = _decomposer(systems, residue_map, lig).decompose(pos)
        elec[eps_in] = sum(r["electrostatic"] for r in result.values())
        assert 2 * elec[eps_in] == pytest.approx(_openmm_delta(systems, pos, lig, (0,)), abs=TOL)
        assert check["max_solvation_deviation"] < TOL
        assert check["receptor_solvation"] + check["ligand_solvation"] == pytest.approx(
            _openmm_delta(systems, pos, lig, POLAR_GROUPS), abs=TOL)
    assert abs(elec[1.0]) > 1.0
    assert elec[2.0] == pytest.approx(elec[1.0] / 2.0, rel=1e-9)


@pytest.mark.parametrize("amber_compatible", [False, True])
def test_amber_compatible_selects_the_numerical_surface_nonpolar_term(dna, amber_compatible):
    """Default: the per-residue non-polar term is the run's own LCPO term per
    atom and sums to it. amber_compatible: MMPBSA.py's decomposition
    convention (sander gbsa=2, a numerical surface), computed as up to 0.0.22:
    Shrake-Rupley areas of the complex minus the receptor alone, x 0.0072."""
    import mdtraj as md
    from mmgbsa.decomposition import (NONPOLAR_FORCE, NONPOLAR_NUMERICAL_SURFACE,
                                      PerResidueDecomposition)

    class Calculator:
        decomposition_platform = None

    Calculator.amber_compatible = amber_compatible
    systems = _systems(dna, model="OBC2", salt=0.15, sa="LCPO", typing="element")
    pos, lig = dna["complex"][1], dna["lig"]
    residue_map = _residue_map(dna["complex"][0](), lig)
    rec = np.setdiff1d(np.arange(len(pos)), lig)
    plat = openmm.Platform.getPlatformByName("Reference")
    sys_dict = {"complex_system": systems[0], "protein_system": systems[1], "ligand_system": systems[2]}
    for key, system in zip(("complex", "protein", "ligand"), systems):
        sys_dict[f"{key}_context"] = openmm.Context(system, openmm.VerletIntegrator(1e-3), plat)
    top = md.load_prmtop(os.path.join(DNA, "complex.prmtop"))
    decomposer = PerResidueDecomposition(Calculator())._make_decomposer(sys_dict, residue_map, lig, md_topology=top)
    result, check = decomposer.decompose(pos)
    nonpolar = np.array([result[r]["nonpolar_solvation"] for r in residue_map])

    # polar, vdW and electrostatic follow the run's forces in both modes
    assert check["decomposed_binding_energy_without_nonpolar"] == pytest.approx(
        check["openmm_binding_energy_without_nonpolar"], abs=TOL)
    if amber_compatible:
        assert decomposer.nonpolar == NONPOLAR_NUMERICAL_SURFACE
        frame = md.Trajectory(pos.astype(np.float32)[None], top)
        # b61ec93 _sa_decomposition, line by line
        complex_a2 = md.shrake_rupley(frame, mode="atom")[0] * 100.0
        receptor_a2 = md.shrake_rupley(frame.atom_slice(rec), mode="atom")[0] * 100.0
        correction = np.zeros(frame.n_atoms)
        correction[rec] = complex_a2[rec] - receptor_a2
        expected = np.array([np.sum(correction[atoms]) * 0.0072 for atoms in residue_map.values()])
        np.testing.assert_allclose(nonpolar, expected, rtol=0, atol=1e-12)
        assert np.abs(expected).max() > 0.05
        # not the LCPO term, as in Amber
        assert abs(check["decomposed_nonpolar"] - check["openmm_nonpolar"]) > 0.01
    else:
        assert decomposer.nonpolar == NONPOLAR_FORCE
        assert check["decomposed_nonpolar"] == pytest.approx(check["openmm_nonpolar"], abs=TOL)
        assert check["decomposed_binding_energy"] == pytest.approx(check["openmm_binding_energy"], abs=TOL)


# ---------------------------------------------------------------------------
# Output format and failures
# ---------------------------------------------------------------------------

def test_frame_rows_keep_every_component_and_the_chain():
    from mmgbsa.decomposition import FRAME_COMPONENTS, _frame_row

    energies = {c: 1.0 for c in FRAME_COMPONENTS}
    row = _frame_row({"ARG_10_0": energies, "ARG_10_1": dict(energies, total=2.0)}, 3, 1, 5)
    assert row["ARG10_0_total"] == 1.0 and row["ARG10_1_total"] == 2.0
    assert {f"ARG10_1_{c}" for c in FRAME_COMPONENTS} <= set(row)
    assert (row["frame_index"], row["frame_number"], row["total_frames"]) == (3, 1, 5)


def test_failed_decomposition_raises():
    from mmgbsa.decomposition import DecompositionError, PerResidueDecomposition

    class Calc:
        systems = {"complex_system": openmm.System(), "complex_topology": None}
        salt_concentration = 0.15

    decomp = PerResidueDecomposition(Calc())
    with pytest.raises(DecompositionError):
        decomp._prepare_decomposition_systems(None, None, None)

    def fail(*a, **k):
        raise DecompositionError("boom")

    decomp._perform_per_residue_decomposition = fail
    with pytest.raises(DecompositionError):
        decomp.run_per_residue_analysis(None, None, None, None, skip_baseline=True,
                                         baseline_results={"mean_binding_energy": 0.0})
