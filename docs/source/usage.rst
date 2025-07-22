Basic Usage
===========

This page provides a quick start for running OpenGBSA analyses.

Quick Start
-----------

1. Prepare your input files (see :doc:`config_guide`)
2. Create a YAML configuration file
3. Run the analysis:

.. code-block:: bash

   opengbsa my_config.yaml

Minimal YAML Configuration Example
----------------------------------

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