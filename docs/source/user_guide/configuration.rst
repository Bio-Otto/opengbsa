Configuration reference
=======================

This page lists every option of the configuration file with its default and
its allowed values. Each name links to the section of the user guide that
explains what the option does, how it is applied and when to change it.
``opengbsa --create-complete-config config.yaml`` writes the same options into
a commented file. An option that is left out takes its default. The sections
``input_files``, with ``complex_pdb`` and ``trajectory``, and
``analysis_settings`` are required; ``analysis_settings`` may be empty, so
that every option of it takes its default.

input_files
-----------

.. list-table::
   :header-rows: 1
   :widths: 28 16 56

   * - Option
     - Default
     - Values
   * - :ref:`complex_pdb <opt-complex_pdb>`
     - required
     - ``.prmtop``, ``.parm7``, ``.top``, ``.tpr``, ``.psf``, ``.xml``, ``.pdb``, ``.gro``
   * - :ref:`trajectory <opt-trajectory>`
     - required
     - any format MDTraj reads
   * - :ref:`solvated_topology <opt-solvated_topology>`
     - none
     - ``.gro``, ``.pdb``, ``.prmtop``
   * - :ref:`receptor_topology, ligand_topology <opt-receptor_topology>`
     - none
     - Amber ``.prmtop`` or ``.parm7``
   * - :ref:`ligand_mol <opt-ligand_mol>`
     - none
     - ``.sdf`` or ``.mol2``
   * - ``ligand_pdb``, ``ligand_resname``
     - none
     - a ``.pdb`` file; a residue name (:doc:`input_formats`)

output_settings
---------------

.. list-table::
   :header-rows: 1
   :widths: 28 16 56

   * - Option
     - Default
     - Values
   * - :ref:`output_directory <opt-output_directory>`
     - ``mmgbsa_results``
     - a folder; ``--output-dir`` on the command line overrides it

analysis_settings
-----------------

.. list-table::
   :header-rows: 1
   :widths: 28 16 56

   * - Option
     - Default
     - Values
   * - :ref:`binding_mode <opt-binding_mode>`
     - ``standard``
     - ``standard``, ``ppi``
   * - :ref:`ligand_selection, receptor_selection <opt-ligand_selection>`
     - none
     - MDTraj selections
   * - :ref:`gb_model <opt-gb_model>`
     - ``OBC2``
     - ``HCT``, ``OBC1``, ``OBC2``, ``GBn``, ``GBn2``
   * - :ref:`gb_radii_set <opt-gb_radii_set>`
     - ``auto``
     - ``auto``, ``mbondi``, ``mbondi2``, ``mbondi3``, ``bondi``
   * - :ref:`salt_concentration <opt-salt_concentration>`
     - ``0.15``
     - 0 to 2 mol/L
   * - :ref:`temperature <opt-temperature>`
     - ``300.0``
     - above 0, up to 1000 K
   * - :ref:`solute_dielectric <opt-solute_dielectric>`
     - ``1.0``
     - above 0
   * - :ref:`solvent_dielectric <opt-solvent_dielectric>`
     - ``78.5``
     - above 0
   * - :ref:`sa_model <opt-sa_model>`
     - ``LCPO``
     - ``LCPO``, ``ACE``
   * - :ref:`lcpo_typing <opt-lcpo_typing>`
     - ``element``
     - ``element``, ``amber``
   * - :ref:`rgbmax <opt-rgbmax>`
     - ``null``
     - ``null`` or above 0 Å
   * - :ref:`amber_compatible <opt-amber_compatible>`
     - ``false``
     - ``true``, ``false``
   * - ``charge_method``
     - ``am1bcc``
     - ``am1bcc``, ``mmff94``, ``gasteiger`` (:ref:`opt-ligand_mol`)
   * - :ref:`frame_start, frame_end <opt-frame_start>`
     - whole trajectory
     - frame indices from 0; ``frame_end`` excluded
   * - :ref:`frame_stride <opt-frame_stride>`
     - ``1``
     - a positive whole number
   * - :ref:`max_frames <opt-max_frames>`
     - ``50``
     - 1 to 10,000, or ``null`` for all
   * - :ref:`frame_selection <opt-frame_selection>`
     - ``sequential``
     - ``sequential``, ``equidistant``, ``random``
   * - :ref:`random_seed <opt-random_seed>`
     - ``42``
     - a whole number
   * - :ref:`reimage_trajectory <opt-reimage_trajectory>`
     - ``true``
     - ``true``, ``false``
   * - :ref:`run_per_residue_decomposition <opt-run_per_residue_decomposition>`
     - ``false``
     - ``true``, ``false``
   * - :ref:`decomp_frames <opt-decomp_frames>`
     - ``10``
     - a positive whole number, at most ``max_frames``
   * - :ref:`report_raw_energies <opt-report_raw_energies>`
     - ``false``
     - ``true``, ``false``
   * - :ref:`plot_top_residues <opt-plot_top_residues>`
     - ``10``
     - a positive whole number
   * - :ref:`energy_decomposition, decomposition_method <opt-energy_decomposition>`
     - ``false``, ``full``
     - ``true``, ``false``; ``full``, ``fast``
   * - :ref:`run_entropy_analysis <opt-run_entropy_analysis>`
     - ``false``
     - ``true``, ``false``
   * - :ref:`entropy_method <opt-entropy_method>`
     - ``interaction``
     - ``interaction``, ``c2``, ``quasiharmonic``, ``nmode``, ``none``
   * - :ref:`qha_analyze_complex <opt-qha_analyze_complex>`
     - ``false``
     - ``true``, ``false``
   * - :ref:`nmode_frames <opt-nmode_frames>`
     - ``10``
     - a positive whole number, or ``null`` for all analysed frames
   * - :ref:`nmode_drms <opt-nmode_drms>`
     - ``0.001``
     - above 0 kcal/mol/Å
   * - :ref:`nmode_maxcyc <opt-nmode_maxcyc>`
     - ``10000``
     - a positive whole number
   * - :ref:`nmode_precision <opt-nmode_precision>`
     - ``mixed``
     - ``mixed``, ``double``
   * - :ref:`max_workers <opt-max_workers>`
     - one less than the cores, at most 4
     - a positive whole number
   * - :ref:`print_interval <opt-print_interval>`
     - ``10``
     - a positive whole number
   * - :ref:`verbose <opt-verbose>`
     - ``true``
     - ``true``, ``false``

forcefield_settings
-------------------

.. list-table::
   :header-rows: 1
   :widths: 28 16 56

   * - Option
     - Default
     - Values
   * - ``protein_forcefield``
     - ``amber``
     - ``amber`` (ff14SB), ``amber14``, ``charmm``, ``charmm_gromacs``; Coordinate Mode only (:ref:`coordinate-mode`)
   * - ``charmm_params``
     - none
     - CHARMM parameter files for ``.psf`` input (:ref:`charmm-namd`)
   * - ``charmm_coordinates``
     - none
     - coordinates for ``.psf`` input (:ref:`charmm-namd`)

platform_settings, advanced_settings
------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 28 16 56

   * - Option
     - Default
     - Values
   * - :ref:`preferred_platform <opt-preferred_platform>`
     - automatic
     - ``CUDA``, ``OpenCL``, ``CPU``, ``Reference``
   * - :ref:`decomposition_platform <opt-decomposition_platform>`
     - as the calculation
     - as ``preferred_platform``
   * - :ref:`strict_platform <opt-strict_platform>`
     - ``false``
     - ``true``, ``false``
   * - :ref:`n_jobs <opt-n_jobs>` (``advanced_settings``)
     - ``-1``
     - ``-1`` for every core, or a positive whole number

reporting_settings
------------------

.. list-table::
   :header-rows: 1
   :widths: 28 16 56

   * - Option
     - Default
     - Values
   * - :ref:`include_all_plots <opt-reporting>`
     - ``true``
     - ``true``, ``false``
   * - :ref:`export_format <opt-reporting>`
     - ``png``
     - a Plotly image format
   * - :ref:`export_width, export_height <opt-reporting>`
     - ``4200``, ``2800``
     - pixels

Equivalents in Amber MMPBSA.py
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - MMPBSA.py
     - OpenGBSA
   * - ``-cp``, ``-rp``, ``-lp``
     - ``complex_pdb``, ``receptor_topology``, ``ligand_topology``
   * - ``-sp``
     - ``solvated_topology``
   * - ``-y``
     - ``trajectory``
   * - ``receptor_mask``, ``ligand_mask``
     - ``receptor_selection``, ``ligand_selection``, in MDTraj syntax
   * - ``startframe``, ``endframe``
     - ``frame_start`` = ``startframe`` - 1, ``frame_end`` = ``endframe``
   * - ``interval``
     - ``frame_stride``
   * - ``igb`` = 1, 2, 5, 7, 8
     - ``gb_model``: ``HCT``, ``OBC1``, ``OBC2``, ``GBn``, ``GBn2``
   * - ``saltcon``
     - ``salt_concentration``
   * - ``surften`` 0.0072, ``surfoff`` 0
     - ``sa_model: LCPO``, which uses these values
   * - fixed ``intdiel`` 1, ``extdiel`` 80.0 or 78.3, ``rgbmax`` 25 Å
     - ``solute_dielectric``, ``solvent_dielectric``, ``rgbmax``; ``amber_compatible`` sets MMPBSA.py's values
   * - ``&decomp``
     - ``run_per_residue_decomposition`` (:doc:`decomposition`)
   * - ``&nmode`` ``nmstartframe``, ``nmendframe``, ``nminterval``
     - ``nmode_frames``
   * - ``&nmode`` ``drms``, ``maxcyc``
     - ``nmode_drms``, ``nmode_maxcyc``
   * - ``entropy = 1`` (quasi-harmonic)
     - ``entropy_method: quasiharmonic``

.. TODO decomposition idecomp equivalence after the decomposition validation

Unknown, renamed and removed options
------------------------------------

``opengbsa`` checks the configuration before the run. It stops with a message
when a required file is missing or when an option has a value that is not
allowed, and the message says which values are allowed. An option or section
that OpenGBSA does not know gives a warning with the closest known name, for
example ``unknown option 'gb_modle' in analysis_settings; did you mean
'gb_model'?``, and is otherwise ignored. An option written in the wrong section
gives a warning that names the right one. The run writes the configuration it
used, with every default filled in, to ``analysis_config.yaml``
(:doc:`outputs`), so a run can be repeated from that file.

Older spellings that still work, with a warning:

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Old
     - Use instead
   * - top-level ``params: output_directory``
     - ``output_settings: output_directory``
   * - top-level ``input:`` section, with ``topology``
     - ``input_files:``, with ``complex_pdb``
   * - ``analysis_settings: per_residue_decomposition``
     - ``analysis_settings: run_per_residue_decomposition``
   * - ``advanced_settings: run_per_residue_decomposition``
     - ``analysis_settings: run_per_residue_decomposition``
   * - ``binding_mode: dimer_ligand``
     - ``binding_mode: standard``, which it equals

Options removed in version 0.0.23. They had no effect, or did not work; each
now gives a warning that says why, except ``nonbonded_cutoff``, which stops the
run when it is set to a distance:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Option
     - Why
   * - ``analysis_settings: nonbonded_cutoff``
     - It cut some interactions and not others and stopped GPU runs. MM/GBSA
       evaluates every atom pair, as MMPBSA.py does; ``rgbmax`` limits the
       Born-radius calculation.
   * - ``analysis_settings: use_cache``
     - The cache never reused anything within or between runs.
   * - ``analysis_settings: parallel_processing``
     - CPU runs always divide the frames among worker processes;
       ``max_workers`` sets how many.
   * - ``analysis_settings: dimer_mode``, ``receptor_subunits``
     - The interface between two units is computed by a second run with
       ``receptor_selection`` and ``ligand_selection``.
   * - ``output_settings: analysis_name``; ``analysis_settings: compound_name``,
       ``save_frame_by_frame_csv``, ``frame_by_frame_csv_name``,
       ``include_residue_summary``, ``frame_output_components``,
       ``frame_output_format``; ``reporting_settings: generate_final_report``,
       ``include_plots``, ``include_statistics``, ``include_convergence``;
       section ``reproducibility_settings``
     - They were read only by a second program that ``opengbsa`` never ran.
       The files they named are always written.
   * - ``reporting_settings: entropy_approximation``
     - The report shows the result of ``entropy_method``.
   * - section ``visualization_settings``
     - It was never read. The PyMOL script and the heat maps are always
       written.
   * - ``input_files: receptor_resname``, ``ligand_charge``
     - They were never read.
