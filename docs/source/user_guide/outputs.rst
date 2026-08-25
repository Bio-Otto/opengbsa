Outputs and Results
===================

OpenGBSA generates a comprehensive set of results in the specified output directory. This section details the file structure and provides guidelines on interpreting the data.

Directory Structure
-------------------

.. code-block:: text

    results/
    ├── analysis_YYYYMMDD_HHMMSS/
    │   ├── final_report.txt               # Executive summary text file
    │   ├── mmgbsa_results.csv             # Time-series global energies (Total, VdW, Elec, etc.)
    │   ├── interactive_report.html        # ⚡ The Main Interactive Dashboard
    │   ├── mmgbsa.log                     # Detailed execution log
    │   ├── figures/                       # Static PNG plots (High Resolution)
    │   │   ├── binding_energy_plot.png
    │   │   ├── energy_components.png
    │   │   └── per_residue_decomposition.png
    │   └── decompositions/                # (If enabled)
    │       ├── per_residue_ligand.csv     # Ligand residue contributions (Energy vs Residue)
    │       ├── per_residue_receptor.csv   # Receptor residue contributions
    │       └── decomposition_summary.csv  # Mean/Std for each residue

Interactive Report
------------------

The ``interactive_report.html`` is the primary interface for exploring your results. It aggregates all data into a single, portable file.

1.  **Summary Table**: Key statistics (Mean Binding Energy, Standard Deviation).
2.  **Dynamic Plots**: Zoomable, pannable Plotly charts for time-series data.
3.  **3D Structure**: An embedded NGLView/3Dmol viewer showing the complex with mapped binding hotspots.
4.  **Residue Table**: A sortable, searchable table of per-residue energy components.

Interpreting Results
--------------------

Global Binding Free Energy
~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``mmgbsa_results.csv`` file contains the time series of the binding free energy components.

*   **Total Binding Energy (:math:`\Delta G_{bind}`)**:
    *   **Negative Value**: Indicates favorable binding. The more negative, the stronger the affinity.
    *   **Positive Value**: Indicates unfavorable binding (no affinity).
*   **VDW vs Electrostatics**:
    *   Hydrophobic binders are typically driven by :math:`\Delta E_{vdw}` and non-polar solvation (:math:`\Delta G_{SA}`).
    *   Polar/Salt-bridge driven binders show strong favorable electrostatic energy (:math:`\Delta E_{elec}`), often partially offset by the penalty of desolvation (:math:`\Delta G_{GB}`).

Per-Residue Decomposition
~~~~~~~~~~~~~~~~~~~~~~~~~
Decomposition analysis identifies which residues contribute most to binding.

*   **Hotspots**: Residues with highly negative interaction energies (e.g., < -1.0 kcal/mol) are key binding determinants.
*   **Interpret with Caution**: Decomposition energies are effective pairwise interactions and do not strictly sum to the total binding free energy due to the non-additive nature of the GB solvation term. They are best used qualitatively to identify important contacts.

CSV Formats
-----------

**mmgbsa_results.csv**

Columns correspond to energy terms per frame:

   - ``delta_total``: Final binding free energy (~ :math:`\Delta H - T\Delta S`).
   - ``delta_vdw``: Van der Waals contribution.
   - ``delta_elec``: Electrostatic contribution.
   - ``delta_gb``: Polar solvation (Generalized Born).
   - ``delta_sa``: Non-polar solvation (Surface Area).

**per_residue_*.csv**

Matrix format where:

   - **Rows**: Trajectory Frames.
   - **Columns**: Residue IDs (e.g., ``ALA:12``).
   - **Values**: Total interaction energy of that residue with the binding partner.
