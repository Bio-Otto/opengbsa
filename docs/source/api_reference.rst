API Reference
=============

Core Analysis Engine
--------------------

.. module:: mmgbsa.core

.. class:: GBSACalculator(system_generator, complex_topology, complex_positions, temperature=300*unit.kelvin, ...)

   The main engine for performing fixed-trajectory MM/GBSA calculations with enhanced force group handling.

   :param system_generator: OpenMM SystemGenerator object containing forcefield definitions.
   :param complex_topology: OpenMM Topology of the receptor-ligand complex.
   :param complex_positions: Atomic positions of the complex.
   :param temperature: Simulation temperature (default: 300 K).
   :param gb_model: Generalized Born model ('OBC2', 'GBn', 'GBn2').

   .. method:: run_comprehensive(self, trajectory_file, frame_indices=None)

      Executes the MM/GBSA energy analysis on the provided trajectory.

      1.  **Decomposition**: Splits the complex into Receptor and Ligand subsets.
      2.  **Force Grouping**: Assigns specific Force Groups to internal (10-14), nonbonded (0), and GBSA (1-4) forces.
      3.  **Energy Calculation**: Evaluates $E_{MM} + G_{solv}$ for Complex, Receptor, and Ligand for each frame.
      4.  **Binding Energy**: $\Delta G_{bind} = E_{complex} - E_{receptor} - E_{ligand}$.

      :param trajectory_file: Path to the .xtc/.dcd trajectory.
      :param frame_indices: List of specific frame indices to analyze; if None, analyzes all.
      :return: Dictionary containing mean binding energies and standard deviations.

   .. method:: calculate_interaction_entropy(self, binding_energies, temperature=300.0)

      Calculates the Entropic Penalty using the **Interaction Entropy (IE)** method.
      
      .. math:: -T\Delta S = \frac{1}{\beta} \ln \langle e^{\beta(E - \langle E \rangle)} \rangle

      :param binding_energies: Array of binding energies ($\Delta H$) from the trajectory frames.
      :param temperature: Temperature in Kelvin.
      :return: Entropy penalty value (positive kcal/mol).

   .. method:: assign_force_groups(sys_obj)

      **Static Helper**. Iterates through the OpenMM System and assigns Force Groups to ensure correct energy component separation.

      *   **Group 0**: Standard Nonbonded (Coulomb + VdW).
      *   **Group 1**: Standard GBSAOBCForce.
      *   **Group 2**: CustomGBForce (GBn2 part 1).
      *   **Group 3**: CustomNonbondedForce (Salt Screening).
      *   **Group 4**: CustomBondForce (Surface Area / Non-polar).
      *   **Group 10-14**: Harmonic Bond, Angle, Torsion, CMPA (Internal Energies).

Running the Analysis
--------------------

.. module:: mmgbsa.runner

.. function:: run_analysis(config_file)

   Main entry point. Parses the YAML configuration, builds the system using `GBSACalculator`, and executes the specific analysis protocol (e.g., Single Trajectory Protocol).

   :param config_file: Path to the YAML configuration file.
