API Reference
=============

This section documents the public API of OpenGBSA.

Core Analysis
-------------

.. module:: mmgbsa.core

.. class:: GBSACalculator

   The main engine for MM/GBSA calculations.

   .. method:: __init__(temperature=300.0, verbose=1, gb_model='OBC2', salt_concentration=0.15, charge_method='am1bcc', solute_dielectric=1.0, solvent_dielectric=78.5, entropy_method='interaction', decomposition_method='full', protein_forcefield='amber', use_cache=True, visualization_settings=None, platform=None, reporting_settings=None, sa_model='ACE')

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

.. class:: MMGBSARunner(config, output_dir)

   Handles the execution flow based on a configuration dictionary.

   .. method:: run_analysis()
   
      Executes the analysis defined in the configuration. Use this for programmatic access.

Configuration
-------------

.. module:: mmgbsa.config

.. class:: ConfigManager(config_file)

   Manages validation and loading of YAML configuration files.

   .. method:: validate_config()
   
      Checks if the configuration output meets all schema requirements.

Reporting
---------

.. module:: mmgbsa.reporting

.. class:: HTMLReportGenerator(output_dir)

   Generates interactive HTML reports with Plotly charts.

   .. method:: generate_report(results, plots, output_filename='report.html')
   
      Compiles all analysis data into a single, shareable HTML file.
