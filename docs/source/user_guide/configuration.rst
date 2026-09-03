Configuration Reference
=======================

The OpenGBSA configuration is a standard YAML file. This reference guide details every parameter, its units, default values, and scientific context.

Input Files
-----------

Defines the topology and trajectory data for the system.

.. list-table::
   :widths: 20 15 65
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``complex_pdb``
     - File Path
     - **Required**. The reference structure. Can be a PDB file or an AMBER ``.prmtop`` file. If using a trajectory, this file provides the topology.
   * - ``trajectory``
     - File Path
     - **Required**. The MD trajectory containing snapshots to analyze. Supports ``.xtc``, ``.dcd``, ``.netcdf``.
   * - ``ligand_mol``
     - File Path
     - **Conditional**. Path to the ligand file (SDF or Mol2). Required if using a PDB ``complex_pdb`` to correctly parameterize the ligand.
   * - ``receptor_topology``
     - File Path
     - **Optional**. Specific topology for the receptor (e.g., ``receptor.prmtop``). Used in Amber-native mode.
   * - ``ligand_topology``
     - File Path
     - **Optional**. Specific topology for the ligand (e.g., ``ligand.prmtop``). Used in Amber-native mode.

Analysis Settings
-----------------

Controls the physics model and sampling protocols.

``gb_model``
~~~~~~~~~~~~
*   **Type**: String
*   **Options**: ``OBC1``, ``OBC2``, ``GBn``, ``GBn2``, ``HCT``
*   **Default**: ``OBC2``
*   **Description**: Selects the Generalized Born model for polar solvation. ``OBC2`` is generally recommended for protein-ligand systems due to its optimized rescaling parameters. See :doc:`../theory/gb_models` for details.

``salt_concentration``
~~~~~~~~~~~~~~~~~~~~~~
*   **Type**: Float
*   **Units**: Molar (M)
*   **Default**: ``0.15`` (Physiological salt)
*   **Description**: Defines the ionic strength of the solvent. This value is used to calculate the Debye screening length (:math:`\kappa`) which attenuates electrostatic interactions.

``entropy_method``
~~~~~~~~~~~~~~~~~~
*   **Type**: String
*   **Options**: ``none``, ``interaction``, ``quasiharmonic``, ``normal_mode``
*   **Default**: ``none``
*   **Description**: Specifies the method to estimate the entropic penalty (:math:`-T\Delta S`). 
    *   ``interaction``: Fast, good for rigid binding.
    *   ``normal_mode``: Slow, standard reference.
    *   ``none``: Returns :math:`\Delta G` without entropy (often sufficient for ranking).

``start_frame``, ``end_frame``, ``interval``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*   **Type**: Integer
*   **Description**: Controls the subset of trajectory frames to analyze.
    *   ``start_frame``: Index of the first frame (0-based).
    *   ``end_frame``: Index of the last frame.
    *   ``interval``: Stride. Analyze every Nth frame. Increasing this (e.g., to 10 or 100) significantly speeds up analysis with minimal loss of statistical accuracy if the trajectory is converged.

``run_per_residue_decomposition``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*   **Type**: Boolean
*   **Default**: ``true``
*   **Description**: If true, decomposes the binding free energy into contributions from each receptor residue. Essential for identifying "hotspot" residues critical for binding.

``nonbonded_cutoff``
~~~~~~~~~~~~~~~~~~~~
*   **Type**: Float or ``null``
*   **Units**: Angstrom
*   **Default**: ``null`` (no cutoff -- exact ``O(N^2)`` nonbonded/GB Born-radius evaluation)
*   **Description**: Optional distance cutoff for both the nonbonded and GB
    Born-radius calculations (OpenMM's ``CutoffNonPeriodic`` nonbonded
    method, analogous to Amber ``sander``'s ``rgbmax``). Leaving this unset
    is the most accurate option and matches the values used throughout
    this project's own validation work (see the Tutorial section), but is
    the single most expensive part of the calculation for larger systems --
    setting a cutoff (e.g. ``16``-``25``) trades a small amount of accuracy
    for a significant speedup on systems with more than a few thousand
    atoms. Changing this value changes the computed energies (it is not a
    performance-only knob), so a config's `nonbonded_cutoff` should stay
    fixed across any comparison you intend to make (e.g. against an Amber
    reference computed with a specific ``rgbmax``).
