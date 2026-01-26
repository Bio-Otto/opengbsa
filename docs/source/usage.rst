Usage Guide
===========

This guide explains how to prepare and execute a new MM/GBSA analysis.

1. Required Input Files
-----------------------

To start a new run, you need the following four files:

1.  **Complex Topology**: A PDB file containing both the protein and the ligand.
    *   *Example*: `complex.pdb` or `system.pdb`
    *   *Note*: The ligand residue name must match your config (default: `LIG` or `UNL`).

2.  **Ligand Structure**: An SDF or MOL2 file defining the ligand's bond order and connectivity.
    *   *Example*: `ligand.sdf`
    *   *Importance*: Necessary for accurate parameterization (OpenFF/GAFF).

3.  **Trajectory**: A Simulation trajectory file.
    *   *Example*: `trajectory.xtc` or `trajectory.dcd`
    *   *Requirement*: Must assume the complex topology corresponds to this trajectory. It is recommended to center/image the trajectory first.

4.  **Configuration File**: A YAML file to control the analysis.

2. Configuration (config.yaml)
------------------------------

Create a `config.yaml` file with the following structure. This uses the robust settings (GBn2 + Gasteiger) established in this project.

.. code-block:: yaml

    input_files:
      complex_pdb: "path/to/complex.pdb"      # Protein + Ligand Topology
      trajectory: "path/to/trajectory.xtc"    # MD Simulation Trajectory
      ligand_mol: "path/to/ligand.sdf"        # Ligand with correct bond orders

    output_settings:
      output_directory: "mmgbsa_results"
      analysis_name: "my_analysis_name"
      save_plots: true

    analysis_settings:
      gb_model: "GBn2"              # Recommended: Accurate and Efficient
      charge_method: "gasteiger"    # Recommended: Robust and Fast
      salt_concentration: 0.15      # Molar (M)
      temperature: 300.0            # Kelvin
      
      # Frame Selection
      max_frames: 100               # Number of frames to analyze
      frame_stride: 10              # Skip every N frames (optional)
      frame_selection: "sequential" # or "random"

3. Running the Analysis
-----------------------

Once files are ready, activate your environment and run:

.. code-block:: bash

    # 1. Activate Environment
    conda activate mmgbsa

    # 2. Run Analysis
    mmgbsa config.yaml

4. Output Interpretation
------------------------

Results will be saved in `mmgbsa_results/my_analysis_name/`. Look for:

*   **final_report.txt**: Summary of Binding Energy ($\Delta G$) and Interaction Entropy.
*   **fixed_enhanced_mmgbsa_results_gbn2.csv**: Frame-by-frame energy breakdown.
*   **energy_analysis.png**: Plots of energy convergence.