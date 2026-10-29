#!/usr/bin/env python3
"""Per-atom polar (GB) and non-polar solvation energies of an OpenMM System.

The per-residue decomposition attributes the run's own solvation energy to
atoms. Every parameter is read from the force objects of the System that the
MM/GBSA calculation evaluates: per-particle charges, radii, screening factors
(and GBn2's alpha/beta/gamma), the dielectric constants, kappa, the offset and
the Born-radius cutoff (rgbmax) are all literals or parameters of those forces.
The per-atom energies therefore sum to the force-group energies OpenMM computes
for the same coordinates, for every GB model and setting the run used.

Supported forces
----------------
openmm.CustomGBForce
    The Amber GB models of openmm.app.internal.customgbforces (HCT, OBC1, OBC2,
    GBn, GBn2), including the rgbmax edit of the Born-radius integral
    (GBSAForceManager._apply_rgbmax), salt screening and an ACE term when one is
    present. Computed values and energy terms are evaluated from the force's own
    expression strings (getComputedValueParameters / getEnergyTermParameters)
    with numpy, and the GBn/GBn2 neck tables are the force's own tabulated
    functions (OpenMM's d0/m0 tables interpolated to the radii present).
openmm.GBSAOBCForce
    OBC2 (and its ACE term) as OpenMM's reference platform implements it
    (platforms/reference/src/SimTKReference/ReferenceObc.cpp), including the
    term for an atom engulfed by a neighbour's descreening sphere.
openmm.LCPOForce
    LCPO as in ReferenceLCPOIxn.cpp / ReferenceKernels.cpp: each atom's
    surface tension x (P1 S_i + P2 sum_j A_ij + P3 sum_jk A_jk + P4 sum_j A_ij sum_k A_jk).

Attribution
-----------
A single-particle term belongs to its particle. A pair term is split equally
between the two particles, as Amber's `decpair` does (sander decomp.F90:
hfval = 0.5*fval for every pair). An LCPO term belongs to the atom whose
surface it approximates.

Force groups follow GBSACalculator's convention (`_frame_energy_terms`):
groups 1 and 2 hold the polar GB force, groups 4 and 16 the non-polar one.
"""

import math
import os
import re
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from openmm import openmm, unit

KJ_TO_KCAL = 1.0 / 4.184

POLAR_GROUPS = (1, 2)
NONPOLAR_GROUPS = (4, 16)

# OpenMM's Coulomb constant 1/(4 pi eps0) in kJ nm/(mol e^2), as defined in
# SimTKOpenMMRealType.h and used by GBSAOBCForce. The CustomGBForce models carry
# their own constant (138.935485) inside their expressions.
_E_CHARGE = 1.602176634e-19
_AVOGADRO = 6.02214076e23
_EPSILON0 = 1e-6 * 8.8541878128e-12 / (_E_CHARGE * _E_CHARGE * _AVOGADRO)
ONE_4PI_EPS0 = 1.0 / (4.0 * math.pi * _EPSILON0)

# Elements per block in the pair loops (rows x columns). Each numpy temporary of
# a block is CHUNK_ELEMENTS doubles (2 MB), so memory stays bounded for any
# system size.
CHUNK_ELEMENTS = 1 << 18


class SolvationDecompositionError(RuntimeError):
    """The solvation energy of a System cannot be decomposed exactly."""


def _n_threads():
    try:
        n = int(os.environ.get('OPENGBSA_DECOMP_THREADS', '0'))
    except ValueError:
        n = 0
    if n <= 0:
        n = min(8, os.cpu_count() or 1)
    return max(1, n)


_EXECUTOR = None


def _executor():
    global _EXECUTOR
    if _EXECUTOR is None:
        _EXECUTOR = ThreadPoolExecutor(max_workers=_n_threads())
    return _EXECUTOR


def _map_blocks(func, blocks):
    """Run func over the row blocks; numpy releases the GIL inside its loops."""
    if len(blocks) <= 1 or _n_threads() == 1:
        return [func(*b) for b in blocks]
    return list(_executor().map(lambda b: func(*b), blocks))


def _value(x, u=None):
    """Strip a unit (if any) and return a float."""
    if unit.is_quantity(x):
        x = x.value_in_unit(u) if u is not None else x._value
    return float(x)


# ---------------------------------------------------------------------------
# Lepton expressions -> numpy
# ---------------------------------------------------------------------------

def _step(x):
    return np.where(x >= 0, 1.0, 0.0)


def _delta(x):
    return np.where(x == 0, 1.0, 0.0)


def _select(x, a, b):
    return np.where(x != 0, a, b)


def _erf(x):
    from scipy.special import erf
    return erf(x)


def _erfc(x):
    from scipy.special import erfc
    return erfc(x)


# Lepton's built-in functions (lepton/src/Operation.h), elementwise in numpy.
_BUILTINS = {
    'sqrt': np.sqrt, 'exp': np.exp, 'log': np.log,
    'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
    'sec': lambda x: 1.0 / np.cos(x), 'csc': lambda x: 1.0 / np.sin(x), 'cot': lambda x: 1.0 / np.tan(x),
    'asin': np.arcsin, 'acos': np.arccos, 'atan': np.arctan, 'atan2': np.arctan2,
    'sinh': np.sinh, 'cosh': np.cosh, 'tanh': np.tanh,
    'erf': _erf, 'erfc': _erfc,
    'step': _step, 'delta': _delta, 'select': _select,
    'square': np.square, 'cube': lambda x: x * x * x, 'recip': lambda x: 1.0 / x,
    'min': np.minimum, 'max': np.maximum, 'abs': np.abs,
    'floor': np.floor, 'ceil': np.ceil,
}

_TOKEN = re.compile(r"\s*(?:(?P<num>(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)"
                    r"|(?P<name>[A-Za-z_][A-Za-z0-9_]*)|(?P<op>[-+*/^(),]))")


def _discrete_function(fn):
    """A CustomGBForce tabulated function as a numpy callable (Discrete only)."""
    if isinstance(fn, openmm.Discrete1DFunction):
        values = np.asarray(fn.getFunctionParameters(), dtype=float)
        sizes = (values.size,)
    elif isinstance(fn, openmm.Discrete2DFunction):
        xs, ys, values = fn.getFunctionParameters()
        values, sizes = np.asarray(values, dtype=float), (int(xs), int(ys))
    elif isinstance(fn, openmm.Discrete3DFunction):
        xs, ys, zs, values = fn.getFunctionParameters()
        values, sizes = np.asarray(values, dtype=float), (int(xs), int(ys), int(zs))
    else:
        raise SolvationDecompositionError(
            f"Tabulated function of type {type(fn).__name__} is not supported by the decomposition")

    def lookup(*args):
        if len(args) != len(sizes):
            raise SolvationDecompositionError("Tabulated function called with the wrong number of arguments")
        index, stride = 0, 1
        for a, size in zip(args, sizes):
            k = np.rint(a).astype(np.int64)
            if np.any(k < 0) or np.any(k >= size):
                raise SolvationDecompositionError("Tabulated function argument out of range")
            index = index + stride * k
            stride *= size
        return values[index]
    return lookup


class _Expression:
    """One Lepton expression ("main; name=def; ...") compiled for numpy.

    Lepton parses the definitions after the main expression from the last one
    to the first, and each may use the ones that follow it; they are evaluated
    in that order here. Definitions shadow variables of the same name.
    """

    def __init__(self, text, functions):
        parts = text.split(';')
        self.text = text
        self.constants = {}
        self.functions = functions
        definitions = []
        for part in parts[1:]:
            if not part.strip():
                continue
            name, sep, rhs = part.partition('=')
            name = name.strip()
            if not sep or not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", name):
                raise SolvationDecompositionError(f"Cannot parse definition '{part}' in '{text}'")
            definitions.append((name, rhs))
        self.definitions = [(name, self._compile(rhs)) for name, rhs in reversed(definitions)]
        self.main = self._compile(parts[0])

    def _compile(self, src):
        out, pos, prev_name = [], 0, None
        src = src.rstrip()
        while pos < len(src):
            m = _TOKEN.match(src, pos)
            if not m or m.end() == pos:
                raise SolvationDecompositionError(f"Cannot parse '{src[pos:]}' in expression '{self.text}'")
            pos = m.end()
            if m.group('num') is not None:
                key = f"_K{len(self.constants)}"
                self.constants[key] = np.float64(float(m.group('num')))
                out.append(key)
                prev_name = None
            elif m.group('name') is not None:
                out.append(m.group('name'))
                prev_name = len(out) - 1
            else:
                op = m.group('op')
                if op == '(' and prev_name is not None:
                    fname = out[prev_name]
                    if fname not in self.functions:
                        raise SolvationDecompositionError(f"Unknown function '{fname}' in '{self.text}'")
                    out[prev_name] = f"_f_{fname}"
                elif prev_name is not None:
                    out[prev_name] = f"_v_{out[prev_name]}"
                out.append('**' if op == '^' else op)
                prev_name = None
        if prev_name is not None:
            out[prev_name] = f"_v_{out[prev_name]}"
        code = ' '.join(out)
        try:
            return compile(code, '<lepton>', 'eval')
        except SyntaxError as e:
            raise SolvationDecompositionError(f"Cannot compile expression '{self.text}': {e}")

    def evaluate(self, variables):
        """variables: {name: array or float}; returns the value (broadcast shapes)."""
        ns = {'__builtins__': {}}
        ns.update(self.constants)
        for name, f in self.functions.items():
            ns[f"_f_{name}"] = f
        for name, v in variables.items():
            ns[f"_v_{name}"] = v
        with np.errstate(all='ignore'):
            try:
                for name, code in self.definitions:
                    ns[f"_v_{name}"] = eval(code, ns)
                return eval(self.main, ns)
            except NameError as e:
                raise SolvationDecompositionError(f"Undefined variable in GB expression '{self.text}': {e}")


def _full(v, shape):
    v = np.asarray(v, dtype=float)
    if v.shape != tuple(shape):
        v = np.broadcast_to(v, shape)
    return v


def _distances(a, b):
    d = a[:, None, :] - b[None, :, :]
    return np.sqrt(np.einsum('ijk,ijk->ij', d, d))


# ---------------------------------------------------------------------------
# GB models as expression sets (CustomGBForce, or GBSAOBCForce translated)
# ---------------------------------------------------------------------------

_SINGLE = openmm.CustomGBForce.SingleParticle
_PAIR = openmm.CustomGBForce.ParticlePair
_PAIR_NO_EXCL = openmm.CustomGBForce.ParticlePairNoExclusions


class _ExpressionGBModel:
    """Per-atom energies of a CustomGBForce-like model, evaluated in numpy.

    Mirrors ReferenceCustomGBIxn: computed values in order (single-particle, or a
    sum over every other particle j with particle 1 = i); single-particle energy
    terms per particle; pair energy terms over i < j with particle 1 = i, split
    half to each particle. With a cutoff, pairs at r >= cutoff are skipped
    (`cutoff_inclusive` False: GBSAOBCForce skips r > cutoff only).
    """

    def __init__(self, name, param_names, params, globals_, computed, energy_terms,
                 functions, exclusions=None, cutoff=None, cutoff_inclusive=False):
        self.name = name
        self.n = params.shape[0]
        self.param_names = list(param_names)
        self.params = {p: np.ascontiguousarray(params[:, k], dtype=float) for k, p in enumerate(param_names)}
        self.globals = dict(globals_)
        funcs = dict(_BUILTINS)
        funcs.update(functions)
        self.computed = [(n, _Expression(e, funcs), k) for n, e, k in computed]
        self.energy_terms = [(_Expression(e, funcs), k) for e, k in energy_terms]
        self.cutoff = cutoff
        self.cutoff_inclusive = cutoff_inclusive
        self.exclusions = {}
        for i, j in (exclusions or []):
            self.exclusions.setdefault(i, set()).add(j)
            self.exclusions.setdefault(j, set()).add(i)

    def _single_vars(self, pos, values):
        v = dict(self.globals)
        v.update(self.params)
        v.update(values)
        v.update({'x': pos[:, 0], 'y': pos[:, 1], 'z': pos[:, 2]})
        return v

    def _pair_vars(self, a, b, c0, c1, values, r):
        v = dict(self.globals)
        for src in (self.params, values):
            for name, arr in src.items():
                v[name + '1'] = arr[a:b, None]
                v[name + '2'] = arr[None, c0:c1]
        v['r'] = r
        return v

    def _mask(self, a, b, c0, r, exclude, upper):
        rows = np.arange(a, b)[:, None]
        cols = np.arange(c0, c0 + r.shape[1])[None, :]
        valid = (cols > rows) if upper else (cols != rows)
        if self.cutoff is not None:
            valid &= (r <= self.cutoff) if self.cutoff_inclusive else (r < self.cutoff)
        if exclude and self.exclusions:
            for i in range(a, b):
                for j in self.exclusions.get(i, ()):
                    if c0 <= j < c0 + r.shape[1]:
                        valid[i - a, j - c0] = False
        return valid

    def _row_blocks(self, upper):
        n = self.n
        blocks, a = [], 0
        while a < n:
            cols = (n - a) if upper else n
            m = max(1, CHUNK_ELEMENTS // max(cols, 1))
            b = min(n, a + m)
            blocks.append((a, b))
            a = b
        return blocks

    def _pair_value(self, expr, pos, values, exclude):
        n = self.n

        def work(a, b):
            r = _distances(pos[a:b], pos)
            v = _full(expr.evaluate(self._pair_vars(a, b, 0, n, values, r)), r.shape)
            return np.where(self._mask(a, b, 0, r, exclude, upper=False), v, 0.0).sum(axis=1)

        return np.concatenate(_map_blocks(work, self._row_blocks(upper=False))) if n else np.zeros(0)

    def _pair_energy(self, expr, pos, values, exclude):
        n = self.n
        energy = np.zeros(n)

        def work(a, b):
            r = _distances(pos[a:b], pos[a:])
            v = _full(expr.evaluate(self._pair_vars(a, b, a, n, values, r)), r.shape)
            e = np.where(self._mask(a, b, a, r, exclude, upper=True), v, 0.0)
            return a, b, 0.5 * e.sum(axis=1), 0.5 * e.sum(axis=0)

        for a, b, rows, cols in _map_blocks(work, self._row_blocks(upper=True)):
            energy[a:b] += rows
            energy[a:] += cols
        return energy

    def per_atom(self, pos):
        """Per-particle energies (kJ/mol) at positions pos (n x 3, nm)."""
        values = {}
        for name, expr, kind in self.computed:
            if kind == _SINGLE:
                values[name] = np.array(_full(expr.evaluate(self._single_vars(pos, values)), (self.n,)))
            else:
                values[name] = self._pair_value(expr, pos, values, exclude=(kind == _PAIR))
        energy = np.zeros(self.n)
        for expr, kind in self.energy_terms:
            if kind == _SINGLE:
                energy += _full(expr.evaluate(self._single_vars(pos, values)), (self.n,))
            else:
                energy += self._pair_energy(expr, pos, values, exclude=(kind == _PAIR))
        self.last_values = values
        return energy


def _custom_gb_model(force):
    if force.getNonbondedMethod() == openmm.CustomGBForce.CutoffPeriodic:
        raise SolvationDecompositionError("A periodic GB force cannot be decomposed per atom")
    cutoff = None
    if force.getNonbondedMethod() == openmm.CustomGBForce.CutoffNonPeriodic:
        cutoff = _value(force.getCutoffDistance(), unit.nanometer)
    names = [force.getPerParticleParameterName(k) for k in range(force.getNumPerParticleParameters())]
    n = force.getNumParticles()
    params = np.zeros((n, len(names)))
    for i in range(n):
        p = force.getParticleParameters(i)
        params[i, :] = [_value(x) for x in p]
    globals_ = {force.getGlobalParameterName(k): _value(force.getGlobalParameterDefaultValue(k))
                for k in range(force.getNumGlobalParameters())}
    computed = [tuple(force.getComputedValueParameters(k)) for k in range(force.getNumComputedValues())]
    energy = [tuple(force.getEnergyTermParameters(k)) for k in range(force.getNumEnergyTerms())]
    functions = {force.getTabulatedFunctionName(k): _discrete_function(force.getTabulatedFunction(k))
                 for k in range(force.getNumTabulatedFunctions())}
    exclusions = [tuple(force.getExclusionParticles(k)) for k in range(force.getNumExclusions())]
    return _ExpressionGBModel('CustomGBForce', names, params, globals_, computed, energy,
                              functions, exclusions=exclusions, cutoff=cutoff)


def _gbsa_obc_model(force):
    """GBSAOBCForce as the equivalent expression set (ReferenceObc.cpp)."""
    if force.getNonbondedMethod() == openmm.GBSAOBCForce.CutoffPeriodic:
        raise SolvationDecompositionError("A periodic GB force cannot be decomposed per atom")
    cutoff = None
    if force.getNonbondedMethod() != openmm.GBSAOBCForce.NoCutoff:
        cutoff = _value(force.getCutoffDistance(), unit.nanometer)
    n = force.getNumParticles()
    params = np.zeros((n, 3))
    for i in range(n):
        q, rad, scale = force.getParticleParameters(i)
        params[i] = (_value(q, unit.elementary_charge), _value(rad, unit.nanometer), _value(scale))
    eps_in = _value(force.getSoluteDielectric())
    eps_out = _value(force.getSolventDielectric())
    sa = _value(force.getSurfaceAreaEnergy(), unit.kilojoule_per_mole / unit.nanometer ** 2)
    pref = -ONE_4PI_EPS0 * (1.0 / eps_in - 1.0 / eps_out) if (eps_in != 0 and eps_out != 0) else 0.0
    offset = 0.009
    computed = [
        ('orad', f"radius-{offset!r}", _SINGLE),
        ('srad', f"scale*(radius-{offset!r})", _SINGLE),
        # ReferenceObc::computeBornRadii: the sum there omits the 0.5 that it
        # applies afterwards; the second select is the engulfed-atom term
        # 2*(1/or_i - 1/L) of that loop, halved likewise.
        ('I', "select(step(r+srad2-orad1), 0.5*(1/L-1/U+0.25*(r-srad2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r)"
              "+select(step(srad2-r-orad1), 1/orad1-1/L, 0), 0);"
              "U=r+srad2; L=max(orad1, D); D=abs(r-srad2)", _PAIR_NO_EXCL),
        ('B', "1/(1/orad-tanh(psi-0.8*psi^2+4.85*psi^3)/radius); psi=I*orad", _SINGLE),
    ]
    energy = [(f"0.5*{pref!r}*charge^2/B", _SINGLE)]
    if cutoff is None:
        energy.append((f"{pref!r}*charge1*charge2/f; f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))", _PAIR_NO_EXCL))
    else:
        energy.append((f"{pref!r}*charge1*charge2*(1/f-1/{cutoff!r}); f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))",
                       _PAIR_NO_EXCL))
    if sa != 0.0:
        energy.append((f"select(step(B), {4.0 * math.pi * sa!r}*(radius+0.14)^2*(radius/B)^6, 0)", _SINGLE))
    # Negative constants come out as "-1.0"; wrap them so "a*-1.0" parses.
    energy = [(e.replace(f"{pref!r}", f"({pref!r})"), k) for e, k in energy]
    return _ExpressionGBModel('GBSAOBCForce', ['charge', 'radius', 'scale'], params, {}, computed, energy,
                              {}, cutoff=cutoff, cutoff_inclusive=True)


class _LCPOModel:
    """Per-atom LCPO energies of an openmm.LCPOForce (ReferenceLCPOIxn.cpp)."""

    def __init__(self, force):
        if force.usesPeriodicBoundaryConditions():
            raise SolvationDecompositionError("A periodic LCPO force cannot be decomposed per atom")
        self.name = 'LCPOForce'
        n = force.getNumParticles()
        self.n = n
        p = np.zeros((n, 5))
        for i in range(n):
            p[i] = [_value(x) for x in force.getParticleParameters(i)]
        gamma = _value(force.getSurfaceTension(), unit.kilojoule_per_mole / unit.nanometer ** 2)
        self.radius = p[:, 0]
        self.p = p[:, 1:] * gamma
        self.one_body = 4.0 * math.pi * self.p[:, 0] * self.radius ** 2
        self.active = np.nonzero(self.radius != 0.0)[0]
        self.cutoff = 2.0 * self.radius.max() if n else 0.0

    def per_atom(self, pos):
        from scipy import sparse
        from scipy.spatial import cKDTree

        energy = self.one_body.copy()
        act = self.active
        if act.size < 2:
            return energy
        R = self.radius[act]
        x = pos[act]
        pairs = cKDTree(x).query_pairs(self.cutoff, output_type='ndarray')
        if len(pairs) == 0:
            return energy
        i, j = pairs[:, 0], pairs[:, 1]
        d = x[j] - x[i]
        r = np.sqrt(np.einsum('ij,ij->i', d, d))
        keep = r < R[i] + R[j]
        i, j, r = i[keep], j[keep], r[keep]
        dr = (R[i] ** 2 - R[j] ** 2) / r
        a_ij = math.pi * R[i] * (2.0 * R[i] - r - dr)
        a_ji = math.pi * R[j] * (2.0 * R[j] - r + dr)
        m = act.size
        rows = np.concatenate([i, j])
        cols = np.concatenate([j, i])
        A = sparse.csr_matrix((np.concatenate([a_ij, a_ji]), (rows, cols)), shape=(m, m))
        M = sparse.csr_matrix((np.ones(rows.size), (rows, cols)), shape=(m, m))
        term2 = np.asarray(A.sum(axis=1)).ravel()
        # sum over j in N(i), k in N(i) & N(j) of A_jk, and of A_ij*A_jk
        term3 = np.asarray((M @ A).multiply(M).sum(axis=1)).ravel()
        term4 = np.asarray((A @ A).multiply(M).sum(axis=1)).ravel()
        p = self.p[act]
        energy[act] += p[:, 1] * term2 + p[:, 2] * term3 + p[:, 3] * term4
        return energy


def _force_model(force):
    if isinstance(force, openmm.CustomGBForce):
        return _custom_gb_model(force)
    if isinstance(force, openmm.GBSAOBCForce):
        return _gbsa_obc_model(force)
    if isinstance(force, openmm.LCPOForce):
        return _LCPOModel(force)
    raise SolvationDecompositionError(
        f"{type(force).__name__} in force group {force.getForceGroup()} cannot be decomposed per atom")


class SolvationEnergies:
    """Per-atom polar and non-polar solvation energies of one System.

    Parameters are read once from the System's solvation forces (groups
    POLAR_GROUPS and NONPOLAR_GROUPS); `per_atom(positions)` then returns
    (polar, nonpolar) arrays in kcal/mol whose sums are those groups' energies.
    """

    def __init__(self, system, polar_groups=POLAR_GROUPS, nonpolar_groups=NONPOLAR_GROUPS):
        self.n = system.getNumParticles()
        self.polar, self.nonpolar = [], []
        self.polar_groups, self.nonpolar_groups = tuple(polar_groups), tuple(nonpolar_groups)
        for force in system.getForces():
            group = force.getForceGroup()
            if group in self.polar_groups:
                self.polar.append(_force_model(force))
            elif group in self.nonpolar_groups:
                self.nonpolar.append(_force_model(force))
        if not self.polar:
            raise SolvationDecompositionError(
                f"No GB force in force groups {self.polar_groups} of the System")

    def describe(self):
        return {'polar': [m.name for m in self.polar], 'nonpolar': [m.name for m in self.nonpolar]}

    def per_atom(self, positions):
        if unit.is_quantity(positions):
            positions = positions.value_in_unit(unit.nanometer)
        pos = np.ascontiguousarray(np.asarray(positions, dtype=np.float64))
        if pos.shape != (self.n, 3):
            raise SolvationDecompositionError(
                f"Positions have shape {pos.shape}, the System has {self.n} particles")
        polar = np.zeros(self.n)
        for m in self.polar:
            polar += m.per_atom(pos)
        nonpolar = np.zeros(self.n)
        for m in self.nonpolar:
            nonpolar += m.per_atom(pos)
        return polar * KJ_TO_KCAL, nonpolar * KJ_TO_KCAL


def group_energy(context, groups):
    """OpenMM's energy (kcal/mol) of the given force groups of a Context."""
    mask = 0
    for g in groups:
        mask |= 1 << g
    return context.getState(getEnergy=True, groups=mask).getPotentialEnergy().value_in_unit(
        unit.kilocalories_per_mole)
