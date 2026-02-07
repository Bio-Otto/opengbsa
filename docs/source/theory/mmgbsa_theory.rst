MM/GBSA Theory
==============

Molecular Mechanics / Generalized Born Surface Area (MM/GBSA) is a powerful and widely used method for estimating the relative binding free energies of biomolecular complexes. It balances computational efficiency with accuracy by combining molecular mechanics energies, continuum solvation models, and entropy estimates.

OpenGBSA leverages the high-performance OpenMM engine to calculate these terms rapidly on modern hardware.

Thermodynamic Cycle
-------------------

The binding free energy (:math:`\Delta G_{bind}`) is defined as the difference between the free energy of the complex (:math:`G_{complex}`) and the sum of the free energies of the unbound receptor (:math:`G_{receptor}`) and ligand (:math:`G_{ligand}`):

.. math::
   :label: binding_energy

   \Delta G_{bind} = G_{complex} - (G_{receptor} + G_{ligand})

Since calculating absolute free energies is computationally prohibitive, MM/GBSA employs a thermodynamic cycle. The calculations are typically performed on an ensemble of snapshots extracted from a Molecular Dynamics (MD) trajectory.

The free energy of each state (:math:`G`) is estimated as:

.. math::
   :label: free_energy_decomp

   G = E_{MM} + G_{solv} - T\Delta S

Where:
   - :math:`E_{MM}` is the gas-phase molecular mechanics energy.
   - :math:`G_{solv}` is the solvation free energy.
   - :math:`T\Delta S` is the conformational entropy contribution.

Combining these, the binding free energy is approximated as:

.. math::
   :label: mmgbsa_equation

   \Delta G_{bind} \approx \Delta E_{MM} + \Delta G_{solv} - T\Delta S

Energy Components
-----------------

Molecular Mechanics Energy (:math:`E_{MM}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The molecular mechanics energy represents the enthalpy of the system in the gas phase. It is the sum of bonded and non-bonded interactions:

.. math::
   :label: emm_components

   E_{MM} = E_{bonded} + E_{elec} + E_{vdw}

1. **Bonded Interactions** (:math:`E_{bonded}`):
   Includes bond stretching, angle bending, and torsional (dihedral) terms. In the "single trajectory" approach (where receptor and ligand snapshots are extracted from the complex trajectory), :math:`\Delta E_{bonded}` is often assumed to be zero, as the internal geometries are identical. However, OpenGBSA calculates this explicitly to support multiple-trajectory protocols.

2. **Electrostatic Interactions** (:math:`E_{elec}`):
   Calculated using Coulomb's law. In MM/GBSA, this represents the gas-phase electrostatic energy (usually with no cutoff or a very large cutoff, as PME is not used in the implicit solvent post-processing).

   .. math::

      E_{elec} = \sum_{i<j} \frac{q_i q_j}{4\pi \epsilon_0 r_{ij}}

3. **Van der Waals Interactions** (:math:`E_{vdw}`):
   Modeled using the Lennard-Jones 12-6 potential, capturing steric repulsion and dispersion attraction.

   .. math::

      E_{vdw} = \sum_{i<j} 4\epsilon_{ij} \left[ \left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^{6} \right]

Solvation Free Energy (:math:`G_{solv}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The solvation energy is the energy required to transfer the solute from the gas phase to the solvent. It is decomposed into polar and non-polar parts:

.. math::
   :label: g_solv

   G_{solv} = G_{polar} + G_{nonpolar}

1. **Polar Solvation** (:math:`G_{polar}`):
   This term accounts for the electrostatic interaction between the solute and the implicit solvent. In MM/GBSA, this is solved using the **Generalized Born (GB)** equation, which is an analytical approximation to the Poisson-Boltzmann (PB) equation.

   .. math::

      G_{GB} = - \frac{1}{2} \left( \frac{1}{\epsilon_{in}} - \frac{1}{\epsilon_{out}} \right) \sum_{i,j} \frac{q_i q_j}{f_{GB}(r_{ij})}

   Where :math:`\epsilon_{in}` is the solute dielectric (typically 1.0 or 2.0), :math:`\epsilon_{out}` is the solvent dielectric (80.0 for water), and :math:`f_{GB}` depends on the effective Born radii of the atoms.

2. **Non-Polar Solvation** (:math:`G_{nonpolar}`):
   This term accounts for the cost of cavity formation and van der Waals interactions between the solute and the solvent. It is standardly approximated as linearly dependent on the Solvent Accessible Surface Area (SASA):

   .. math::

      G_{nonpolar} = \gamma \cdot \text{SASA} + b

   Standard values (from AMBER): :math:`\gamma = 0.00542` kcal/mol/Å\ :sup:`2`, :math:`b = 0.92` kcal/mol.

Entropic Contribution (:math:`-T\Delta S`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Binding restricts the conformational freedom of the ligand and receptor, resulting in an entropy penalty (:math:`\Delta S < 0`). This term is crucial for obtaining absolute binding free energies but is computationally expensive and has high statistical uncertainty.

Common approximations include:
- **Interaction Entropy**: Approximated from potential energy fluctuations.
- **Normal Mode Analysis (NMA)**: Based on vibrational frequencies of the minimized structure.
- **Quasi-Harmonic Analysis (QHA)**: Derived from the covariance of atomic fluctuations.

See :doc:`entropy_methods` for detailed methodologies.

References
----------

1.  **Massova, I., & Kollman, P. A. (2000).** Combined molecular mechanical and continuum solvent approach (MM-PBSA/GBSA) to predict ligand binding. *Perspectives in Drug Discovery and Design*, 18(1), 113-135.
2.  **Genheden, S., & Ryde, U. (2015).** The MM/PBSA and MM/GBSA methods to estimate ligand-binding affinities. *Expert Opinion on Drug Discovery*, 10(5), 449-461.
