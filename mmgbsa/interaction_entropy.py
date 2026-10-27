"""
Entropy estimates (-TdS) from the fluctuation of the receptor-ligand
interaction energy.

Both estimates use the gas-phase interaction energy of every frame,
dE_int = dE_vdW + dE_elec (the per-frame `delta_vdw` and `delta_elec`
columns), not the binding enthalpy with its solvation terms:

- Interaction entropy (IE), Duan, Liu & Zhang, J Am Chem Soc 2016, 138, 5722:
  -TdS = kT ln < exp(beta (dE_int - <dE_int>)) >
- Second-order cumulant (C2), Ekberg & Ryde, J Chem Theory Comput 2021,
  17, 5379:  -TdS = sigma^2 / (2 kT)

Ekberg & Ryde showed that IE cannot be converged in practice when the
standard deviation sigma of dE_int exceeds 15 kJ/mol, and that C2 gives
unrealistically large entropies above 25 kJ/mol. Both limits are returned
with the estimates so that callers can flag unreliable values.
"""
import numpy as np
from scipy.special import logsumexp

GAS_CONSTANT = 0.0019872041  # kcal/(mol*K)
IE_SIGMA_LIMIT = 15.0 / 4.184  # kcal/mol, 15 kJ/mol
C2_SIGMA_LIMIT = 25.0 / 4.184  # kcal/mol, 25 kJ/mol

METHOD_LABELS = {'interaction': 'Interaction Entropy', 'c2': 'C2 Entropy'}


def interaction_energies(frames):
    """Per-frame gas-phase interaction energy dE_vdW + dE_elec (kcal/mol).

    `frames` is a frame table (DataFrame or dict of arrays) with the
    `delta_vdw` and `delta_elec` columns that OpenGBSA writes per frame.
    """
    missing = [c for c in ('delta_vdw', 'delta_elec') if c not in frames]
    if missing:
        raise ValueError(
            f"Entropy estimates need the per-frame interaction energy terms; "
            f"missing column(s): {', '.join(missing)}."
        )
    return (np.asarray(frames['delta_vdw'], dtype=float)
            + np.asarray(frames['delta_elec'], dtype=float))


def _checked(energies):
    E = np.asarray(energies, dtype=float).ravel()
    if E.size < 2:
        raise ValueError(f"Entropy estimates need at least 2 frames, got {E.size}.")
    if not np.all(np.isfinite(E)):
        raise ValueError("Entropy estimates need finite interaction energies.")
    return E


def interaction_entropy(energies, temperature):
    """-TdS by interaction entropy (Duan et al. 2016), in kcal/mol.

    The exponential average is evaluated with log-sum-exp, so a frame with
    a large positive fluctuation cannot overflow it.
    """
    E = _checked(energies)
    kT = GAS_CONSTANT * float(temperature)
    dE = E - E.mean()
    return float(kT * (logsumexp(dE / kT) - np.log(E.size)))


def c2_entropy(energies, temperature):
    """-TdS by the second-order cumulant (Ekberg & Ryde 2021), in kcal/mol."""
    E = _checked(energies)
    return float(E.var() / (2.0 * GAS_CONSTANT * float(temperature)))


def entropy_estimates(energies, temperature):
    """Both estimates with sigma(dE_int) and their reliability flags."""
    E = _checked(energies)
    sigma = float(E.std())
    return {
        'temperature': float(temperature),
        'n_frames': int(E.size),
        'sigma_interaction_energy': sigma,
        'interaction_entropy': interaction_entropy(E, temperature),
        'c2_entropy': c2_entropy(E, temperature),
        'interaction_entropy_reliable': sigma <= IE_SIGMA_LIMIT,
        'c2_entropy_reliable': sigma <= C2_SIGMA_LIMIT,
    }


def entropy_result(energies, method, temperature):
    """Entropy result for `method` ('interaction' or 'c2').

    `energies` are per-frame gas-phase interaction energies (see
    `interaction_energies`). Returns the estimates of `entropy_estimates`
    plus `method`, `label` and `value`, the -TdS of the requested method.
    """
    method = method.lower()
    if method not in METHOD_LABELS:
        raise ValueError(f"Unknown interaction-energy entropy method '{method}'; "
                         f"use one of {sorted(METHOD_LABELS)}.")
    result = entropy_estimates(energies, temperature)
    key = 'interaction_entropy' if method == 'interaction' else 'c2_entropy'
    result.update(method=method, label=METHOD_LABELS[method], value=result[key],
                  reliable=result[f'{key}_reliable'])
    return result


def describe(result):
    """One line per estimate for logs and text summaries."""
    lines = [f"sigma(dE_int) = {result['sigma_interaction_energy']:.2f} kcal/mol "
             f"over {result['n_frames']} frames at {result['temperature']:.1f} K"]
    for key, limit, name in (('interaction_entropy', IE_SIGMA_LIMIT, 'Interaction entropy'),
                             ('c2_entropy', C2_SIGMA_LIMIT, 'C2 entropy')):
        flag = '' if result[f'{key}_reliable'] else \
            f"  (unreliable: sigma > {limit:.1f} kcal/mol, Ekberg & Ryde 2021)"
        lines.append(f"{name} (-TdS) = {result[key]:.2f} kcal/mol{flag}")
    return lines
