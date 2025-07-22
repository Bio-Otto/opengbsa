Output Files and Directories
===========================

OpenGBSA generates a structured output directory containing all results, logs, and intermediate files. This page describes each output type and how to interpret the results.

Sample Output Directory Structure
---------------------------------
.. code-block:: text

   mmgbsa_results/
     analysis_YYYYMMDD_HHMMSS/   # Main analysis directory (timestamped)
       plots/                    # All generated plots (PNG, PDF, SVG)
       data/                     # Main result tables, per-residue CSVs, entropy results
       logs/                     # Log files for debugging and reproducibility
       reports/                  # Final reports and summaries (TXT, YAML, PDF)
       cache/                    # Cached intermediate systems and data
       debug/                    # Optional debug files (if enabled)

Output Subdirectories
---------------------

plots/
~~~~~~~
- Contains all visualizations generated during analysis.
- Typical files: ``energy_profile.png``, ``per_residue_decomposition.png``, ``convergence_plot.pdf``
- Formats: PNG, PDF, SVG

.. tip::
   Review these plots to assess convergence, per-residue contributions, and overall trends.

data/
~~~~~~
- Contains main numerical results and tables.
- Typical files:
  - ``summary.csv``: Overall MM/GBSA results (total, complex, receptor, ligand energies)
  - ``per_residue_decomposition.csv``: Per-residue energy breakdown
  - ``frame_by_frame.csv``: Frame-by-frame results (if enabled)
  - ``entropy_results.csv``: Entropy analysis results (if enabled)
- Formats: CSV, HDF5, JSON (depending on config)

.. hint::
   Import these files into Excel, pandas, or other tools for further analysis. ``summary.csv`` is the main result table.

logs/
~~~~~~
- Contains log files for the analysis run.
- Typical files: ``run.log``, ``warnings.log``, ``validation.log``
- Formats: TXT

.. important::
   Check for errors, warnings, and validation results here. Useful for troubleshooting and reproducibility.

reports/
~~~~~~~~
- Contains final reports and summaries.
- Typical files: ``final_report.txt``, ``config_summary.yaml``, ``convergence_report.txt``
- Formats: TXT, YAML, PDF

.. note::
   ``final_report.txt`` provides a human-readable summary of the analysis. ``config_summary.yaml`` records the config used.

cache/
~~~~~~~
- Contains cached intermediate systems, forcefields, and data to speed up repeated runs.
- Typical files: ``system_cache.pkl``, ``forcefield_cache.pkl``
- Formats: PKL, HDF5

.. caution::
   Usually not needed by end users, but can be useful for debugging or advanced workflows.

debug/
~~~~~~~
- Contains optional debug files if debug mode is enabled.
- Typical files: ``intermediate_structures.pdb``, ``energy_components.csv``, ``debug.log``
- Formats: PDB, CSV, TXT

.. warning::
   For advanced troubleshooting and development only. Debug files may be large or contain sensitive intermediate data.

Interpreting Key Results
------------------------
.. tip::
   ``summary.csv``: Main results table. Columns include total energy, complex, receptor, ligand, and optionally entropy.

.. hint::
   ``per_residue_decomposition.csv``: Shows energy contribution of each residue. Useful for identifying hot spots.

.. note::
   ``frame_by_frame.csv``: Shows how results change per frame. Useful for convergence and stability analysis.

.. note::
   ``entropy_results.csv``: Entropy contributions (if enabled).

.. tip::
   Plots: Visualize trends, convergence, and per-residue effects.

Best Practices
--------------
.. tip::
   Always check ``logs/`` for warnings or errors after a run.

.. note::
   Use ``reports/`` for a summary and for reproducibility records.

.. hint::
   For publication, use plots and summary tables from ``plots/`` and ``data/``. 