Output files
============

Every run writes its results to a new folder,
``<output_directory>/analysis_<date>_<time>/`` (or
``<DIR>/mmgbsa_analysis_<date>_<time>/`` with ``--output-dir DIR``), so runs
never overwrite each other. The files below are those of the
:doc:`../quickstart`, which uses the OBC2 model, the per-residue decomposition
and the entropy estimate.

.. code-block:: text

    analysis_<date>_<time>/
    ├── fixed_enhanced_mmgbsa_results_obc2.csv          energies of every frame
    ├── results_summary.yaml                            all results, machine-readable
    ├── final_report.txt                                summary of the run
    ├── mmgbsa_report.txt                               statistics of the binding enthalpy
    ├── analysis_config.yaml                            the configuration as run
    ├── energy_analysis.png                             binding enthalpy plots
    │   (with run_per_residue_decomposition: true)
    ├── per_residue_detailed.csv                        mean contribution of every receptor residue
    ├── frame_by_frame_decomposition.csv                contribution of every residue in every frame
    ├── frame_by_frame_decomposition_residue_summary.csv  column index of the file above
    ├── binding_hot_spots.csv                           residues that contribute most
    ├── decomposition_summary.csv                       sum and top residue
    ├── decomposition_energy_check.csv                  decomposed vs OpenMM energies, per frame
    ├── per_residue_decomposition.png                   hot-spot plots
    ├── per_residue_heatmap.png                         top residues frame by frame
    ├── interactive_report.html                         all results in one web page
    ├── hotspots.pdb, view_binding.pml                  PyMOL view coloured by contribution
    │   (with the optional package PandaMap)
    ├── pandamap_interactions.png                       2D map of the ligand's contacts
    ├── structure_3d.html                               3D view of the binding site, shown in the report
    ├── cache/                                          parameterised systems reused within the run
    └── temp_*.pdb                                      intermediate structures

The binding enthalpy of every frame
-----------------------------------

``fixed_enhanced_mmgbsa_results_<gb model>.csv`` has one row per frame. All
energies are in kcal/mol.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Column
     - Meaning
   * - ``frame``
     - Index of the frame among the evaluated frames (from 0).
   * - ``complex_energy``, ``protein_energy``, ``ligand_energy``
     - Total energy of the complex, the receptor and the ligand.
   * - ``binding_energy``
     - :math:`\Delta H_\mathrm{bind}` = complex - receptor - ligand.
   * - ``delta_vdw``, ``delta_elec``
     - Van der Waals and electrostatic interaction energy in the gas phase.
   * - ``delta_nb``
     - Their sum.
   * - ``delta_gb``
     - Polar solvation (Generalized Born) part of the binding enthalpy.
   * - ``delta_sa``
     - Non-polar solvation (surface-area) part.
   * - ``complex_nb``, ``complex_gb``, ``complex_sa``, ``complex_bond``, ``complex_angle``, ``complex_torsion``
     - The terms of the complex.

``binding_energy`` equals ``delta_vdw + delta_elec + delta_gb + delta_sa``. The
bonded terms cancel because the receptor and ligand coordinates are taken from
the complex (single-trajectory approach).

Summaries
---------

``results_summary.yaml`` holds every result in machine-readable form. Under
``mmgbsa``: ``mean_binding_energy``, ``std_dev`` (frame-to-frame standard
deviation), ``std_error`` (standard error of the mean),
``median_binding_energy``, ``min_binding_energy`` and ``max_binding_energy``,
``n_frames``, ``bootstrap_results`` and ``convergence_check``, the
platform (``platform_requested``, ``platform_used``, ``platform_fallback``)
and, when an entropy estimate was made, ``entropy`` and ``delta_g``. Under
``entropy``: the estimate of each method with the standard deviation of the
interaction energy and the reliability flags. Read it in Python with:

.. code-block:: python

    import yaml
    summary = yaml.safe_load(open("results_summary.yaml"))
    print(summary["mmgbsa"]["mean_binding_energy"], summary["mmgbsa"]["std_error"])

``final_report.txt`` repeats the configuration, the binding enthalpy and the
entropy estimates in plain text. ``mmgbsa_report.txt`` adds the median, the
range, a bootstrap 95% confidence interval of the mean and a convergence check
that compares the first and second halves of the frames.
``analysis_config.yaml`` is the configuration as it was run, with the defaults
filled in.

Per-residue decomposition
-------------------------

``per_residue_detailed.csv`` has one row per receptor residue with the mean
over the decomposed frames and its standard deviation (``_std`` columns):

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Column
     - Meaning
   * - ``residue_id``
     - ``NAME_NUMBER_CHAIN``, for example ``PHE_78_0`` (chain as index).
   * - ``vdw``, ``electrostatic``
     - Gas-phase interaction of the residue with the ligand.
   * - ``solvation``
     - Polar (GB) solvation part of the residue.
   * - ``nonpolar_solvation``
     - Non-polar (surface-area) part of the residue.
   * - ``total``
     - ``vdw + electrostatic + solvation``.
   * - ``n_frames``
     - Number of decomposed frames.

``total`` is the value shown in the plots and the PyMOL view; adding
``nonpolar_solvation`` gives the residue's full share including the surface
term. Negative values favour binding.

Every term comes from the OpenMM systems of the MM/GBSA calculation itself, so
it follows the run's GB model, ``rgbmax``, salt concentration, dielectric
constants, ``sa_model`` and ``lcpo_typing``. The van der Waals and electrostatic
terms are half of each pair energy between the residue and the ligand. The
polar and non-polar terms are the change of the residue's atoms between the
complex and the receptor; an atom's GB energy is its self term plus half of
each of its pair terms. Amber MMPBSA.py splits the terms the same way, but it
decomposes the non-polar term with a numerical surface (sander ``gbsa=2``)
instead of LCPO, so that column differs from Amber's.

The sum over the receptor residues is not the binding enthalpy: the ligand
keeps the other half of each pair term and its own solvation change. With the
ligand's share added, the decomposition gives the binding energy of every
decomposed frame term by term. ``decomposition_energy_check.csv`` lists, per
frame, the receptor and ligand sums, the decomposed binding energy, OpenMM's
binding energy for the same coordinates and the largest difference between
the per-atom solvation energies and OpenMM's GB and surface-area energies. A
failed decomposition is reported as an error and writes no per-residue files.
``frame_by_frame_decomposition.csv`` gives the same terms for every frame, in
columns named ``<NAME><NUMBER>_<CHAIN>_<term>``.
``binding_hot_spots.csv`` lists the residues with the largest contributions,
and ``decomposition_summary.csv`` gives their sum, the top residue and the mean
decomposed and OpenMM binding energies.

Plots and reports
-----------------

The plots, the interactive HTML report and the PyMOL view are shown in
:doc:`visualization`.

Intermediate files
------------------

``temp_from_traj.pdb`` (first frame), ``temp_from_prmtop.pdb`` and
``temp_fixed_for_panda.pdb`` are structures that the run writes for the
viewers; they can be deleted. A converted GROMACS topology is written next to
the ``.top`` file, and a converted OpenMM System to ``openmm_system/`` in the
results folder.
