# Theoretical Background

This software performs **End-Point Free Energy Calculations** using the **MM/GBSA** (Molecular Mechanics / Generalized Born Surface Area) method. It simulates the thermodynamic cycle of ligand binding to estimate the **Binding Free Energy** ($\Delta G_{bind}$).

## Thermodynamic Cycle

The binding free energy is calculated as the difference between the free energy of the complex and the sum of the free energies of the unbound receptor and ligand:

$$
\Delta G_{bind} = G_{complex} - (G_{receptor} + G_{ligand})
$$

Where each free energy term ($G$) is estimated as:

$$
G = H - TS \approx E_{MM} + G_{solv} - T\Delta S
$$

## Energy Components

The total approximate free energy is decomposed into the following physical components:

1. **Molecular Mechanics Energy ($E_{MM}$)**:
   Describes the internal and interaction energies in vacuum.

   $$
   E_{MM} = E_{bonded} + E_{electrostatic} + E_{vdw}
   $$

   *   **$E_{bonded}$**: Bond stretching, angle bending, and torsional rotation energies (Force Group 10-14).
   *   **$E_{electrostatic}$**: Coulombic interactions between atomic partial charges (Group 0).
   *   **$E_{vdw}$**: Van der Waals interactions (Lennard-Jones potential) (Group 0).

2. **Solvation Free Energy ($G_{solv}$)**:
   Describes the free energy cost of transferring the molecule from vacuum to solvent.

   $$
   G_{solv} = G_{polar} + G_{nonpolar}
   $$

   *   **Polar Solvation ($G_{polar}$)**: Calculated using the **GBn2** (Onufriev-Bashford-Case GBSA II) Implicit Solvent model. This approximates the Poisson-Boltzmann equation.
   *   **Non-Polar Solvation ($G_{nonpolar}$)**: Correlates with the Solvent Accessible Surface Area (SASA).

   $$
   G_{nonpolar} = \gamma \cdot SASA + b
   $$

   Where $\gamma$ is the surface tension coefficient (typically 0.00542 kcal/mol/Å²).

## Interaction Entropy ($-T\Delta S$)

Standard MM/GBSA typically neglects the conformational entropy change or uses expensive Normal Mode Analysis (NMA). This implementation uses the **Interaction Entropy (IE)** method [Duan et al., JACS 2016], which is rigorously derived from statistical mechanics to estimate the entropy penalty from the fluctuation of binding energies.

The entropy contribution is calculated as:

$$
-T\Delta S = kT \ln \left\langle e^{\beta \Delta E_{int}} \right\rangle
$$

Where:
*   $\beta = 1/kT$
*   $\Delta E_{int} = E_{int} - \langle E_{int} \rangle$ is the fluctuation of the interaction energy.
*   The bracket $\langle \dots \rangle$ denotes the ensemble average.

This method avoids the harmonic approximation of NMA and captures anharmonic entropy contributions, providing a more rigorous estimate of the entropic penalty upon binding.

## Force Field Specifications

*   **Protein Force Field**: AMBER14 (`amber14-all.xml` usually paired with `amber14/tip3p.xml`).
*   **Ligand Force Field**: OpenFF/SMIRNOFF (via `openff-toolkit`) or GAFF (via `antechamber`).
*   **Charges**: AM1-BCC (Standard) or Gasteiger (Fast/Robust approximation).
