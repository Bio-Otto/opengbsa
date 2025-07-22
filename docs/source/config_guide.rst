Configuration Guide
===================

This page provides a comprehensive guide to all YAML configuration options in OpenGBSA.

YAML Structure Overview
-----------------------
A typical configuration file consists of several main sections:

- ``input_files``: Paths to all required input files
- ``analysis_settings``: Main analysis parameters
- ``output_settings``: Output directory and file options
- ``advanced_settings``: Optional advanced analysis parameters
- ``platform_settings``: Hardware and platform preferences
- ``validation_settings``: Input/result validation options
- ``reporting_settings``: Report and summary options
- ``debug_settings``: Debug and logging options
- ``reproducibility_settings``: Reproducibility and environment tracking
- ``performance_settings``: Performance and resource control

Example: Minimal Complete YAML
-----------------------------
.. code-block:: yaml

   input_files:
     ligand_mol: "test/ligand.sdf"
     complex_pdb: "test/complex.pdb"
     ligand_pdb: "test/ligand.pdb"
     trajectory: "test/complex.xtc"

   analysis_settings:
     temperature: 300
     gb_model: "OBC2"
     max_frames: 50
     run_per_residue_decomposition: false
     run_entropy_analysis: false

   output_settings:
     output_directory: "mmgbsa_results"
     save_plots: true

input_files
-----------
- ``ligand_mol``: Path to ligand structure (SDF, MOL2, PDB)
- ``complex_pdb``: Path to protein-ligand complex PDB
- ``ligand_pdb``: Path to isolated ligand PDB
- ``trajectory``: Path to MD trajectory (XTC, DCD, TRR)

analysis_settings
-----------------
- ``temperature``: Simulation temperature (K)
- ``gb_model``: GB model (OBC2, OBC1, HCT, GBn, GBn2)
- ``max_frames``: Max number of frames to analyze (optional, auto-calculated if not set)
- ``frame_start``, ``frame_end``, ``frame_stride``: Frame selection (see :doc:`advanced`)
- ``frame_selection``: 'sequential', 'random', or 'equidistant'
- ``run_per_residue_decomposition``: Enable per-residue decomposition (true/false)
- ``decomp_frames``: Number of frames for decomposition (optional)
- ``run_entropy_analysis``: Enable entropy analysis (true/false)
- ``save_frame_by_frame_csv``: Save frame-by-frame decomposition (true/false)
- ``frame_output_format``: 'csv', 'json', or 'hdf5'
- ``random_seed``: Random seed for reproducibility

output_settings
---------------
- ``output_directory``: Where to save results
- ``save_plots``: Save plots (true/false)
- ``save_trajectories``: Save processed trajectories (true/false)
- ``output_format``: Output file format (csv, txt, yaml)

advanced_settings
-----------------
- ``minimization_tolerance``: Energy minimization tolerance
- ``nma_quality_threshold``: Quality threshold for normal mode analysis
- ``hot_spot_threshold``: Threshold for hot spot detection
- ``bootstrap_confidence``: Confidence for bootstrapping
- ``bootstrap_samples``: Number of bootstrap samples
- ``convergence_window``: Window size for convergence analysis

platform_settings
-----------------
- ``preferred_platform``: 'CUDA', 'CPU', or 'OpenCL'
- ``cuda_device``: CUDA device index
- ``cuda_precision``: 'single' or 'double'
- ``deterministic_forces``: Use deterministic forces (true/false)

validation_settings
-------------------
- ``validate_inputs``: Validate input files (true/false)
- ``validate_results``: Validate results (true/false)
- ``max_std_dev``: Max allowed standard deviation
- ``min_successful_frames``: Minimum number of successful frames
- ``check_convergence``: Check convergence (true/false)
- ``generate_validation_report``: Generate validation report (true/false)

reporting_settings
------------------
- ``generate_final_report``: Generate final report (true/false)
- ``include_config_summary``: Include config summary (true/false)
- ``include_convergence``: Include convergence analysis (true/false)
- ``include_hot_spots``: Include hot spot analysis (true/false)
- ``report_format``: 'txt', 'yaml', etc.

debug_settings
--------------
- ``debug_mode``: Enable debug mode (true/false)
- ``save_energy_components``: Save energy components (true/false)
- ``save_forcefield_params``: Save forcefield parameters (true/false)
- ``save_intermediate_systems``: Save intermediate systems (true/false)
- ``save_snapshots``: Save snapshots (true/false)
- ``verbose_errors``: Verbose error output (true/false)

reproducibility_settings
------------------------
- ``global_random_seed``: Global random seed
- ``save_configuration``: Save config file (true/false)
- ``save_environment``: Save environment info (true/false)
- ``save_system_hashes``: Save system hashes (true/false)
- ``save_versions``: Save package versions (true/false)
- ``strict_reproducibility``: Enforce strict reproducibility (true/false)

performance_settings
--------------------
- ``batch_size``: Batch size for analysis
- ``cpu_affinity``: CPU affinity (optional)
- ``gpu_memory_fraction``: Fraction of GPU memory to use
- ``memory_limit``: Memory limit (GB)
- ``optimize_memory``: Optimize memory usage (true/false)
- ``track_progress``: Show progress bars (true/false)

Best Practices
--------------
- Use only the parameters you need; defaults are set for most options.
- For advanced use, see the example configs in the ``config/`` directory.
- Always check the output directory for logs and validation reports.

.. note::
   For the most up-to-date and detailed config examples, see ``config/COMPLETE_CONFIG_GUIDE.md`` in the repository. 