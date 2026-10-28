"""
Normal-mode entropy of binding (-TdS), as in Amber MMPBSA.py's nmode calculation.

For each selected frame the complex, the receptor and the ligand are each
minimized with their MM/GB energy function until the root-mean-square
gradient falls below `drms`. The non-polar surface term is left out, as in
MMPBSA.py. The normal modes come from a finite-difference Hessian of the
OpenMM forces in double precision. The rigid-rotor harmonic-oscillator model
of an ideal gas at 1 atm then gives the translational, rotational and
vibrational entropy of each species, and

    -TdS = -T (S_complex - S_receptor - S_ligand).

Constants and conventions follow Amber's nmode program (mmpbsa_py_nabnmode;
the same thermo routine is DataSet_Modes::Thermo in cpptraj): frequencies use
108.587 cm^-1 per sqrt(kcal/mol/A^2/amu), the six lowest modes are skipped,
vibrational temperatures use c = 3.0e10 cm/s, the rotational symmetry number
is 1, and a mode with 1 - exp(-hv/kT) <= 1e-7, an imaginary mode included,
contributes no entropy. With these, the entropies printed by
mmpbsa_py_nabnmode are reproduced from its own frequencies to 0.001 cal/(mol K).

MMPBSA.py minimizes in HCT GB (nmode_igb=1) without salt (nmode_istrng=0),
with a Born-radius cutoff of 15 A (NAB's default rgbmax) and a solvent
dielectric of 78.5. `amber_nmode_system` builds that energy function from a
run's own system; by default OpenGBSA keeps the run's GB model instead.
MMPBSA.py also fixes the temperature of its normal-mode entropy at 298.15 K
(AMBER_THERMO_TEMPERATURE); by default OpenGBSA uses the run temperature.
"""
import copy
import re

import numpy as np
import openmm
from openmm import unit

AMBER_FREQ_CONV = 108.587       # cm^-1 per sqrt(kcal/mol/A^2/amu), Amber's nmode
AMBER_NMODE_RGBMAX = 15.0       # A, NAB's default Born-radius cutoff
AMBER_NMODE_SOLVENT_DIELECTRIC = 78.5
AMBER_THERMO_TEMPERATURE = 298.15   # K, fixed by MMPBSA.py and NAB's thermo output

# Physical constants of Amber's nmode thermo routine (SI)
_BOLTZ = 1.380622e-23           # J/K
_PLANCK = 6.626196e-34          # J s
_AVOG = 6.022169e23             # 1/mol
_JPCAL = 4.18674                # J/cal
_AMU = 1.660531e-27             # kg
_PSTD = 1.01325e5               # Pa (1 atm)
_C_VIB = 3.0e10                 # cm/s, as in Amber's vibrational temperatures
_GAS = _AVOG * _BOLTZ           # J/(mol K)

SA_GROUPS = (4, 16)             # non-polar solvation (mmgbsa_core._frame_energy_terms)
GB_GROUPS = (1, 2)              # polar GB
_KJ_NM_PER_KCAL_A = 41.84       # 1 kcal/mol/A in kJ/mol/nm


def thermo(frequencies, masses, coordinates, temperature):
    """Ideal-gas RRHO entropies (cal/(mol K)) of one minimized molecule.

    `frequencies` are all 3N normal-mode frequencies in cm^-1 in ascending
    order (imaginary ones negative), `masses` in amu, `coordinates` in A.
    Returns the translational, rotational, vibrational and total entropy.
    """
    m = np.asarray(masses, dtype=float)
    T = float(temperature)
    weight = m.sum() * _AMU
    kT = _BOLTZ * T
    arg = (kT ** 1.5 / _PLANCK) / _PSTD * (kT / _PLANCK)
    arg *= (2.0 * np.pi) ** 1.5 * (weight / _PLANCK) * np.sqrt(weight) * np.exp(2.5)
    s_trans = _GAS * np.log(arg)
    s_rot = s_vib = 0.0
    if len(m) > 1:
        x = np.asarray(coordinates, dtype=float)
        x = x - (m[:, None] * x).sum(0) / m.sum()
        inertia = (m * (x ** 2).sum(1)).sum() * np.eye(3) - (m[:, None, None] * x[:, :, None] * x[:, None, :]).sum(0)
        moments = np.sort(np.linalg.eigvalsh(inertia))            # amu A^2, ascending
        con = _PLANCK / (_BOLTZ * 8.0 * np.pi ** 2) / _AMU * (_PLANCK / 1e-20)
        linear = len(m) == 2
        if linear:
            sigma = 2.0 if m[0] == m[1] else 1.0
            s_rot = _GAS * np.log(T / (con / moments[2]) * np.e / sigma)
        else:
            rtemp = con / moments
            s_rot = _GAS * np.log(np.sqrt(np.pi * np.e ** 3) * np.sqrt(T ** 3 / np.prod(rtemp)))
        tovt = np.asarray(frequencies, dtype=float)[5 if linear else 6:] * (_PLANCK / _BOLTZ * _C_VIB) / T
        argd = 1.0 - np.exp(-tovt)
        ok = argd > 1e-7
        s = np.zeros_like(tovt)
        s[ok] = tovt[ok] / np.expm1(tovt[ok]) - np.log(argd[ok])
        s_vib = _GAS * s.sum()
    S = {'translational': s_trans / _JPCAL, 'rotational': s_rot / _JPCAL, 'vibrational': s_vib / _JPCAL}
    S['total'] = S['translational'] + S['rotational'] + S['vibrational']
    return S


def gb_offset(force):
    """Dielectric offset (nm) of an Amber-style CustomGBForce (0.009; GBn2 0.0195141)."""
    for k in range(force.getNumComputedValues()):
        found = re.search(r"offset=([0-9.eE+-]+)", force.getComputedValueParameters(k)[1])
        if found:
            return float(found.group(1))
    for k in range(force.getNumEnergyTerms()):
        found = re.search(r"offset=([0-9.eE+-]+)", force.getEnergyTermParameters(k)[0])
        if found:
            return float(found.group(1))
    return 0.009


def gb_particle_parameters(force):
    """(charge e, intrinsic radius nm, screening factor) of every atom of a GB force."""
    out = []
    if isinstance(force, openmm.GBSAOBCForce):
        for i in range(force.getNumParticles()):
            q, r, s = force.getParticleParameters(i)
            out.append((q.value_in_unit(unit.elementary_charge), r.value_in_unit(unit.nanometer), float(s)))
        return out
    offset = gb_offset(force)
    for i in range(force.getNumParticles()):
        q, rho, sr = force.getParticleParameters(i)[:3]
        out.append((q, rho + offset, sr / rho if rho else 0.0))
    return out


def restrict_born_radius_integral(gb_force, rgbmax):
    """Limit the effective-Born-radius integral to pairs within rgbmax (A), as Amber does.

    Same rule as mmgbsa_core's rgbmax option: atoms beyond rgbmax + s_j are
    skipped and atoms straddling the cutoff are integrated up to rgbmax.
    """
    for k in range(gb_force.getNumComputedValues()):
        name, expr, kind = gb_force.getComputedValueParameters(k)
        if kind != openmm.CustomGBForce.ParticlePairNoExclusions or "U=r+sr2;" not in expr:
            continue
        new = expr.replace("select(step(r+sr2-or1), ", "select(step(r+sr2-or1)*step(rgbmax-L), ", 1)
        new = new.replace("U=r+sr2;", "U=min(r+sr2, rgbmax);", 1) + f";rgbmax={rgbmax * 0.1:.10g}"
        if new.count("rgbmax") != 3:
            raise ValueError(f"Unexpected Born-radius integral for rgbmax: {expr}")
        gb_force.setComputedValueParameters(k, name, new, kind)
        return
    raise ValueError("No Born-radius integral found to apply rgbmax to.")


def nmode_system(system):
    """Copy of a species System without its non-polar solvation forces.

    Normal modes need every bond term: a System with bond constraints has
    dropped those terms, and a massless particle has no vibration.
    """
    if system.getNumConstraints():
        raise ValueError("Normal-mode entropy needs a System without bond constraints "
                         f"({system.getNumConstraints()} constraints found).")
    if any(system.getParticleMass(i).value_in_unit(unit.dalton) <= 0 for i in range(system.getNumParticles())):
        raise ValueError("Normal-mode entropy needs a mass on every particle.")
    new = openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(system))
    for i in reversed(range(new.getNumForces())):
        if new.getForce(i).getForceGroup() in SA_GROUPS:
            new.removeForce(i)
    return new


def amber_nmode_system(system, topology=None, solute_dielectric=1.0):
    """MMPBSA.py's nmode energy function built from a species System.

    The GB force is replaced by HCT (Amber igb=1) with the same charges and
    radii, a solvent dielectric of 78.5, no salt and a Born-radius cutoff of
    15 A; the non-polar term is removed. HCT and OBC forces keep their
    screening factors, which are the topology's. GBn and GBn2 carry their own
    screening factors, so for them the HCT factors come from Amber's element
    table, which needs the species `topology`.
    """
    from openmm.app.internal.customgbforces import GBSAHCTForce, _screen_parameter

    new = nmode_system(system)
    gb = [i for i in range(new.getNumForces()) if new.getForce(i).getForceGroup() in GB_GROUPS
          and isinstance(new.getForce(i), (openmm.GBSAOBCForce, openmm.CustomGBForce))]
    if len(gb) != 1:
        raise ValueError(f"Expected one GB force in the species system, found {len(gb)}.")
    force = new.getForce(gb[0])
    params = gb_particle_parameters(force)
    neck = isinstance(force, openmm.CustomGBForce) and any(
        force.getPerParticleParameterName(k) == 'radindex' for k in range(force.getNumPerParticleParameters()))
    if neck:
        if topology is None or topology.getNumAtoms() != len(params):
            raise ValueError("A GBn/GBn2 system needs its own topology for the HCT screening factors.")
        params =[(q, r, _screen_parameter(atom)[0]) for (q, r, _), atom in zip(params, topology.atoms())]
    new.removeForce(gb[0])
    hct = GBSAHCTForce(solventDielectric=AMBER_NMODE_SOLVENT_DIELECTRIC,
                       soluteDielectric=solute_dielectric, SA=None, cutoff=None, kappa=0.0)
    for q, r, s in params:
        hct.addParticle([q, r, s])
    hct.finalize()
    restrict_born_radius_integral(hct, AMBER_NMODE_RGBMAX)
    hct.setForceGroup(2)
    new.addForce(hct)
    return new


def rms_gradient(context):
    """Root-mean-square force component (kcal/mol/A), Amber's RMSG."""
    f = context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(
        unit.kilocalorie_per_mole / unit.angstrom)
    return float(np.sqrt(np.mean(f ** 2)))


def minimize(context, drms=0.001, maxcyc=10000, rounds=5):
    """L-BFGS until the RMS gradient is at most drms (kcal/mol/A).

    OpenMM's tolerance is the same RMS force criterion. The minimizer is
    restarted up to `rounds` times if it stops early; `maxcyc` bounds the
    iterations of each round. Returns the final RMS gradient.
    """
    rmsg = rms_gradient(context)
    for _ in range(rounds):
        if rmsg <= drms:
            break
        openmm.LocalEnergyMinimizer.minimize(context, drms * _KJ_NM_PER_KCAL_A, maxcyc)
        rmsg = rms_gradient(context)
    return rmsg


def hessian(context, positions, step=1e-5):
    """Finite-difference Hessian (kJ/mol/nm^2) from central differences of the forces.

    `positions` (nm, N x 3) is the minimized structure; `step` is the
    displacement in nm. The result is symmetrized.
    """
    pos = np.array(positions, dtype=float)
    n3 = pos.size
    H = np.empty((n3, n3))
    kj_nm = unit.kilojoule_per_mole / unit.nanometer
    for a in range(len(pos)):
        for d in range(3):
            p = pos.copy()
            p[a, d] += step
            context.setPositions(p)
            f_plus = context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(kj_nm).ravel()
            p[a, d] -= 2.0 * step
            context.setPositions(p)
            f_minus = context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(kj_nm).ravel()
            H[3 * a + d] = (f_minus - f_plus) / (2.0 * step)
    context.setPositions(pos)
    H += H.T
    H *= 0.5
    return H


def frequencies(H, masses):
    """Normal-mode frequencies (cm^-1, ascending, imaginary negative) of a Hessian in kJ/mol/nm^2."""
    w = 1.0 / np.sqrt(np.repeat(np.asarray(masses, dtype=float), 3))
    H *= w[:, None]
    H *= w[None, :]
    lam = np.linalg.eigvalsh(H) / 418.4                        # kcal/mol/A^2/amu
    return np.sign(lam) * np.sqrt(np.abs(lam)) * AMBER_FREQ_CONV


def species_entropy(contexts, positions, masses, temperature, drms=0.001, maxcyc=10000):
    """Minimize one species and return its normal-mode entropy terms.

    `contexts` comes from `make_contexts`. With a mixed-precision context the
    structure is first minimized there, which is fast but stops short of drms,
    then finished in double precision, and the Hessian uses mixed precision.
    `positions` in nm. Returns a dict with the RMS gradient, the minimized
    energy (kcal/mol), the number of imaginary modes beyond the six lowest and,
    when the minimization converged, the entropies of `thermo`.
    """
    exact, fast = contexts['double'], contexts.get('mixed')
    pos = np.asarray(positions, dtype=float)
    if fast is not None:
        fast.setPositions(pos)
        openmm.LocalEnergyMinimizer.minimize(fast, drms * _KJ_NM_PER_KCAL_A, maxcyc)
        pos = fast.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    exact.setPositions(pos)
    rmsg = minimize(exact, drms, maxcyc)
    state = exact.getState(getEnergy=True, getPositions=True)
    result = {'rms_gradient': rmsg, 'converged': rmsg <= drms,
              'energy': state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)}
    if not result['converged']:
        return result
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    H = hessian(fast, pos, FD_STEP['mixed']) if fast is not None else hessian(exact, pos, FD_STEP['double'])
    freq = frequencies(H, masses)
    result['n_imaginary'] = int((freq[6:] < 0).sum())
    result['entropy'] = thermo(freq, masses, pos * 10.0, temperature)
    return result


# Finite-difference displacement (nm) for each force precision
FD_STEP = {'double': 1e-5, 'mixed': 3e-4}


def make_contexts(system, platform_name, precision='mixed', device_properties=None):
    """Contexts for one species: {'double': ..., 'mixed': ... or None}.

    The final minimization always runs in double precision. On CUDA and
    OpenCL, precision 'mixed' adds a mixed-precision context for the first
    minimization and the Hessian. Other platforms use the double-precision
    Reference platform, since the CPU platform computes forces in single
    precision.
    """
    def context(platform, props):
        return openmm.Context(system, openmm.VerletIntegrator(0.001), platform, props)

    if platform_name in ('CUDA', 'OpenCL'):
        platform = openmm.Platform.getPlatformByName(platform_name)
        props = dict(device_properties or {})
        return {'double': context(platform, {**props, 'Precision': 'double'}),
                'mixed': context(platform, {**props, 'Precision': 'mixed'}) if precision == 'mixed' else None}
    return {'double': context(openmm.Platform.getPlatformByName('Reference'), {}), 'mixed': None}


def select_frames(n_frames, n_select):
    """Indices of `n_select` frames spread evenly over `n_frames` (None: all)."""
    if n_select is None or n_select >= n_frames:
        return list(range(n_frames))
    return sorted({int(round(x)) for x in np.linspace(0, n_frames - 1, n_select)})


def binding_entropy(frames, temperature):
    """-TdS (kcal/mol) per frame and its mean from `species_entropy` results.

    `frames` is a list of dicts with 'complex', 'receptor' and 'ligand'
    results. Frames in which any species did not converge are left out.
    """
    T = float(temperature)
    rows = []
    for k, fr in enumerate(frames):
        row = {'frame': fr.get('frame', k)}
        ok = all(fr[s]['converged'] for s in ('complex', 'receptor', 'ligand'))
        row['converged'] = ok
        for s in ('complex', 'receptor', 'ligand'):
            row[f'{s}_rms_gradient'] = fr[s]['rms_gradient']
            row[f'{s}_n_imaginary'] = fr[s].get('n_imaginary')
        if ok:
            for term in ('translational', 'rotational', 'vibrational', 'total'):
                dS = (fr['complex']['entropy'][term] - fr['receptor']['entropy'][term]
                      - fr['ligand']['entropy'][term])
                row[f'TdS_{term}'] = T * dS / 1000.0
            row['minus_TdS'] = -row['TdS_total']
        rows.append(row)
    values = np.array([r['minus_TdS'] for r in rows if r['converged']])
    summary = {'n_frames': len(rows), 'n_converged': int(values.size), 'temperature': T}
    if values.size:
        summary.update(value=float(values.mean()),
                       sd=float(values.std(ddof=1)) if values.size > 1 else 0.0,
                       sem=float(values.std(ddof=1) / np.sqrt(values.size)) if values.size > 1 else 0.0)
        for term in ('translational', 'rotational', 'vibrational'):
            summary[f'TdS_{term}'] = float(np.mean([r[f'TdS_{term}'] for r in rows if r['converged']]))
    return rows, summary
