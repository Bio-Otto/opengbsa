Entropy Methods
===============

The entropic contribution (:math:`-T\Delta S`) is the most challenging component of free energy calculations. It represents the loss of translational, rotational, and vibrational degrees of freedom upon binding. OpenGBSA offers four methods to estimate this term.

1. Interaction Entropy (IE)
---------------------------

The Interaction Entropy method approximates the entropy change directly from the fluctuations of the gas-phase receptor--ligand interaction energy along the MD trajectory, :math:`\Delta E_{int} = \Delta E_{vdW} + \Delta E_{elec}` (the per-frame ``delta_vdw`` and ``delta_elec`` columns). The solvation terms of the binding enthalpy are not part of :math:`\Delta E_{int}`, as in Duan et al. and gmx_MMPBSA.

**Theory**:
Derived from the rigorous Zwanzig relation, it avoids the assumption of Gaussian distribution of conformations.

.. math::

   -T\Delta S_{IE} = kT \ln \langle e^{\beta \Delta E_{int}^{fluct}} \rangle

Where :math:`\Delta E_{int}^{fluct} = \Delta E_{int} - \langle \Delta E_{int} \rangle`.

*   **Pros**:
    *   **Computationally Efficient**: Requires only the interaction energies (already calculated for enthalpy) and simple post-processing.
    *   **No Diagonalization**: Does not require matrix diagonalization.
*   **Cons**:
    *   **Convergence**: The exponential average cannot be converged in practice when the standard deviation :math:`\sigma` of :math:`\Delta E_{int}` exceeds 15 kJ/mol (3.6 kcal/mol). OpenGBSA reports :math:`\sigma` and flags larger values.
*   **Citation**: Duan, L., et al. (2016). Interaction entropy: A new paradigm for computing binding free energy. *Journal of the American Chemical Society*, 138(17), 5722-5728.

Second-order cumulant (C2)
~~~~~~~~~~~~~~~~~~~~~~~~~~

Expanding the exponential average to second order gives

.. math::

   -T\Delta S_{C2} = \frac{\sigma^2}{2kT}

which converges much faster than IE but becomes unrealistically large when :math:`\sigma` exceeds 25 kJ/mol (6.0 kcal/mol). Select it with ``entropy_method: c2``. OpenGBSA reports both estimates with :math:`\sigma` whichever method is selected.

*   **Citation**: Ekberg, V., & Ryde, U. (2021). On the use of interaction entropy and related methods to estimate binding entropies. *Journal of Chemical Theory and Computation*, 17(8), 5379-5391.

2. Quasi-Harmonic Analysis (QHA)
--------------------------------

QHA (or Principal Component Analysis - PCA) estimates the configurational entropy by fitting a multi-variate Gaussian to the distribution of conformations sampled in the trajectory.

**Theory**:
It diagonalizes the mass-weighted covariance matrix of atomic fluctuations. The entropy is calculated from the eigenvalues (:math:`\lambda_i`) of this matrix.

.. math::

   S_{config} \approx k_B \sum_{i} \left[ \frac{1}{2}\ln(1 + \frac{k_B T e^2}{\hbar^2 \omega_i^2}) \right]

*   **Pros**:
    *   **Anharmonicity**: Captures some anharmonic motions and correlations between atoms.
*   **Cons**:
    *   **Sampling**: Requires that the simulation is long enough to sample the relevant phase space. If the trajectory is too short (fewer frames than degrees of freedom), the entropy will be underestimated.
*   **Citation**: Karplus, M., & Kushick, J. N. (1981). Method for estimating the configurational entropy of macromolecules. *Macromolecules*, 14(2), 325-332.

3. Normal Mode Analysis (NMode)
-------------------------------

NMode estimates vibrational entropy by treating the system as a collection of harmonic oscillators.

**Theory**:
The system is first minimized to a local energy minimum. The Hessian matrix (second derivatives of the potential energy) is calculated and diagonalized to obtain vibrational frequencies (:math:`\nu_i`).

.. math::

   S_{vib} = R \sum_{i} \left[ \frac{h\nu_i/kT}{e^{h\nu_i/kT} - 1} - \ln(1 - e^{-h\nu_i/kT}) \right]

*   **Pros**:
    *   **Standard Reference**: The most established method in MM/PBSA literature.
*   **Cons**:
    *   **Expensive**: Minimization and Hessian diagonalization for large systems is extremely slow.
    *   **Approximation**: The harmonic approximation may effectively fail for highly flexible biomolecules.
