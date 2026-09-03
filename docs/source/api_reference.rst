API Reference
=============

This section documents the public API of OpenGBSA.

Core Analysis
-------------

.. module:: mmgbsa.mmgbsa_core

.. class:: GBSACalculator

   The main engine for MM/GBSA calculations.

   .. method:: __init__(temperature=300, verbose=1, gb_model='OBC2', salt_concentration=0.15, use_cache=True, parallel_processing=False, max_workers=None, protein_forcefield='amber', charge_method='am1bcc', solute_dielectric=1.0, solvent_dielectric=78.5, entropy_method='none', decomposition_method='full', visualization_settings=None, platform=None, reporting_settings=None, sa_model='ACE', cache_dir=None, nonbonded_cutoff=None, reimage_trajectory=True, charmm_params=None, charmm_coordinates=None)

      Initialize the calculator with specific physics and analysis parameters.

   .. method:: run_comprehensive(ligand_mol, complex_pdb, xtc_file, ligand_pdb, max_frames=50, energy_decomposition=False, output_dir=None)

      **Main Entry Point**. Performs the complete MM/GBSA analysis pipeline:
      
      1.  System parameterization (Protein + Ligand)
      2.  Trajectory processing
      3.  Energy calculation (GB + SA + Entropy)
      4.  Results compilation

      :param ligand_mol: Ligand molecule object (RDKit/OpenFF)
      :param complex_pdb: Path to complex structure/topology
      :param xtc_file: Path to trajectory
      :param output_dir: Directory to save results

   .. method:: set_ligand_forcefield(forcefield_name)
   
      Set the forcefield for ligand parameterization ('gaff' or 'openff').

   .. method:: calculate_interaction_entropy(binding_energies, temperature=300.0)
   
      Calculate Interaction Entropy (IE) from binding energy fluctuations.

Execution Management
--------------------

.. module:: mmgbsa.runner

.. class:: MMGBSARunner(config_file, output_dir=None)

   Handles the execution flow based on a configuration dictionary.

   .. method:: run_analysis()
   
      Executes the analysis defined in the configuration. Use this for programmatic access.

Configuration
-------------

.. module:: mmgbsa.config

.. class:: ConfigManager(config_path=None)

   Manages validation and loading of YAML configuration files.

   .. method:: validate_config()

      Checks if the configuration output meets all schema requirements.

Reporting
---------

.. module:: mmgbsa.reporting

.. class:: HTMLReportGenerator(output_dir, config=None)

   Generates interactive HTML reports with Plotly charts.

   .. method:: generate_report(analysis_results, frame_data, global_results=None, complex_pdb_path=None, ligand_resname=None)

      Compiles all analysis data into a single, shareable HTML file.
