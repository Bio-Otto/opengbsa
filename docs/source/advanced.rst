Advanced Usage
==============

This page describes advanced analysis options in OpenGBSA, including frame selection, per-residue decomposition, frame-by-frame output, and entropy analysis.

Frame Selection
---------------
You can analyze a specific range of frames from your trajectory using the following parameters in your YAML config:

.. code-block:: yaml

   analysis_settings:
     frame_start: 200      # Start from frame 200 (0-indexed)
     frame_end: 1000       # End at frame 1000 (inclusive)
     frame_stride: 1       # Analyze every frame (set to 5 for every 5th frame)

- If ``max_frames`` is not set, the number of frames is automatically determined from these parameters.
- You can combine with ``frame_selection: 'sequential'``, ``'random'``, or ``'equidistant'`` for different sampling strategies.

Per-Residue Decomposition
-------------------------
Enable per-residue energy decomposition to identify key binding residues:

.. code-block:: yaml

   analysis_settings:
     run_per_residue_decomposition: true
     decomp_frames: 100   # (optional) Number of frames for decomposition (default: all selected frames)

- The output file ``per_residue_detailed.csv`` contains mean and std for vdw, electrostatic, solvation, and total energy per residue.

Frame-by-Frame Decomposition Output
-----------------------------------
To save per-residue decomposition for each frame:

.. code-block:: yaml

   analysis_settings:
     run_per_residue_decomposition: true
     save_frame_by_frame_csv: true
     frame_by_frame_csv_name: "frame_by_frame_decomposition"
     frame_output_components: ["vdw", "electrostatic", "solvation", "total"]
     frame_output_format: "csv"  # or "json", "hdf5"
     include_residue_summary: true

- Output files: ``frame_by_frame_decomposition.csv``, ``frame_by_frame_decomposition_residue_summary.csv``

Entropy Analysis (Normal Mode Analysis)
---------------------------------------
To include vibrational entropy in your binding free energy calculation:

.. code-block:: yaml

   analysis_settings:
     run_entropy_analysis: true

- Entropy results are included in the final report and results summary.

Combining Advanced Options
--------------------------
You can combine all advanced options in a single YAML config:

.. code-block:: yaml

   analysis_settings:
     frame_start: 200
     frame_end: 1000
     frame_stride: 5
     run_per_residue_decomposition: true
     save_frame_by_frame_csv: true
     run_entropy_analysis: true

.. note::
   For a full list of configuration options, see the :doc:`config_guide` section. 