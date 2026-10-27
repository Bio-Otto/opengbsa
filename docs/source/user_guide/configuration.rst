Configuration reference
=======================

This page lists every option of the configuration file, with its default and
its meaning. ``opengbsa --create-complete-config`` writes the same options into
a commented file. An option that is left out takes the default shown here.
Only ``input_files`` and ``analysis_settings`` are required.

.. contents:: Sections
   :local:
   :depth: 1

input_files
-----------

.. list-table::
   :header-rows: 1
   :widths: 24 14 62

   * - Option
     - Default
     - Meaning
   * - ``complex_pdb``
     - required
     - The system. A topology for Native Mode (``.prmtop``, ``.top``,
       ``.tpr``, ``.psf``, or an OpenMM ``.xml`` with a ``.pdb`` of the same
       name) or a structure for Coordinate Mode (``.pdb``, ``.gro``). See
       :doc:`input_formats`.
   * - ``trajectory``
     - required
     - The MD trajectory, any format MDTraj reads (``.xtc``, ``.trr``,
       ``.dcd``, ``.nc``, ``.mdcrd``/``.trj``, multi-model ``.pdb``).
   * - ``solvated_topology``
     - none
     - Topology of the atoms in the trajectory when it still contains water
       and ions (``.gro``, ``.pdb`` or solvated ``.prmtop``). For ``.psf``
       input it also gives the coordinates.
   * - ``receptor_topology``
     - none
     - Amber topology of the receptor. With ``ligand_topology`` it defines
       receptor and ligand by atom count.
   * - ``ligand_topology``
     - none
     - Amber topology of the ligand.
   * - ``ligand_mol``
     - none
     - Coordinate Mode: the small-molecule ligand as ``.sdf`` (or ``.mol2``)
       with bond orders and partial charges, parameterised with OpenFF.
   * - ``ligand_pdb``
     - none
     - Optional ligand-only ``.pdb``, used only to check that its residue
       names occur in the complex.
   * - ``ligand_resname``
     - none
     - Residue name shown in the ligand box of the interactive 3D viewer. It
       does not choose the ligand; use ``ligand_selection`` for that.

Every other key in this section is taken as a file path that must exist.

output_settings
---------------

.. list-table::
   :header-rows: 1
   :widths: 24 18 58

   * - Option
     - Default
     - Meaning
   * - ``output_directory``
     - ``mmgbsa_results``
     - Each run writes ``<output_directory>/analysis_<date>_<time>/``. The
       command-line option ``--output-dir DIR`` overrides it.

analysis_settings
-----------------

Receptor and ligand
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 14 62

   * - Option
     - Default
     - Meaning
   * - ``binding_mode``
     - ``standard``
     - ``standard``: the ligand is found automatically (see
       :doc:`input_formats`) and the receptor is the rest. ``ppi``: receptor =
       chain 0, ligand = chain 1. ``dimer_ligand`` behaves like ``standard``.
   * - ``ligand_selection``
     - none
     - MDTraj selection of the ligand, for example ``"chainid 1"`` or
       ``"resname LIG"``. Overrides ``binding_mode``.
   * - ``receptor_selection``
     - none
     - MDTraj selection of the receptor. Overrides ``binding_mode``.

Solvation model
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 14 62

   * - Option
     - Default
     - Meaning
   * - ``gb_model``
     - ``OBC2``
     - Generalized Born model: ``HCT``, ``OBC1``, ``OBC2``, ``GBn`` or
       ``GBn2`` (Amber ``igb`` 1, 2, 5, 7 and 8). See :doc:`../theory/gb_models`.
   * - ``gb_radii_set``
     - ``auto``
     - Intrinsic GB radii. ``auto`` keeps the radii of an Amber ``.prmtop``
       and assigns the model's own set for other input (``mbondi`` for HCT and
       GBn, ``mbondi2`` for OBC1 and OBC2, ``mbondi3`` for GBn2). Or force one
       of ``mbondi``, ``mbondi2``, ``mbondi3``, ``bondi``.
   * - ``salt_concentration``
     - ``0.15``
     - Ionic strength in M (0 to 2) for Debye-Hückel screening of the GB term.
   * - ``temperature``
     - ``300.0``
     - Temperature in K (0 to 1000), used by salt screening and the entropy
       estimates.
   * - ``solute_dielectric``
     - ``1.0``
     - Dielectric constant inside the solute.
   * - ``solvent_dielectric``
     - ``78.5``
     - Dielectric constant of the solvent. With ``amber_compatible`` the value
       of MMPBSA.py is used (80.0 for CHARMM and chamber topologies, 78.3
       otherwise).
   * - ``sa_model``
     - ``LCPO``
     - Non-polar term: ``LCPO`` (surface tension 0.0072 kcal/mol/Å², offset 0,
       as in Amber's GB calculations) or ``ACE`` (OpenMM's analytical
       continuum electrostatics approximation, 0.0054 kcal/mol/Å²).
   * - ``lcpo_typing``
     - ``element``
     - How LCPO parameters are assigned: ``element`` (by element and bonding)
       or ``amber`` (by atom-type name, as Amber does). They differ only for
       non-Amber atom types.
   * - ``rgbmax``
     - none
     - Born-radius cutoff in Å: only atoms within this distance enter an
       atom's Born radius. None means every atom pair.
   * - ``nonbonded_cutoff``
     - none
     - Cutoff in Å for the non-bonded and GB interactions. None (or 999 and
       above) means no cutoff, the usual choice for MM/GBSA.
   * - ``amber_compatible``
     - ``false``
     - ``true`` applies Amber MMPBSA.py's conventions at once: ``rgbmax`` 25
       Å, MMPBSA.py's solvent dielectric and ``lcpo_typing: amber``. Use it to
       compare with MMPBSA.py results. See :doc:`../theory/mmgbsa_theory`.

Frames
~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 14 62

   * - Option
     - Default
     - Meaning
   * - ``frame_start``
     - first frame
     - Index of the first frame (counting from 0).
   * - ``frame_end``
     - last frame
     - Index after the last frame (exclusive). Must be larger than
       ``frame_start``.
   * - ``frame_stride``
     - ``1``
     - Use every Nth frame of the range.
   * - ``max_frames``
     - ``50``
     - Number of frames to evaluate (1 to 10,000); ``null`` for every frame
       of the range.
   * - ``frame_selection``
     - ``sequential``
     - ``sequential``: the first ``max_frames`` frames of the range after the
       stride. ``equidistant``: ``max_frames`` frames spread evenly over the
       range. ``random``: ``max_frames`` random frames.
   * - ``random_seed``
     - ``42``
     - Seed of the ``random`` frame selection.
   * - ``reimage_trajectory``
     - ``true``
     - Make every molecule whole and keep receptor and ligand in one periodic
       image before computing energies. Switch it off only for a trajectory
       that is already imaged and has no box.

Per-residue decomposition
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 12 58

   * - Option
     - Default
     - Meaning
   * - ``run_per_residue_decomposition``
     - ``false``
     - Decompose the binding enthalpy into the contributions of the receptor
       residues (van der Waals, electrostatic, polar and non-polar
       solvation), per frame and on average, with plots, the HTML report and
       a PyMOL script.
   * - ``decomp_frames``
     - ``10``
     - Number of frames for the decomposition, from the same selection; must
       not exceed ``max_frames``.
   * - ``report_raw_energies``
     - ``false``
     - Also write the complex, receptor and ligand energies of each residue.
   * - ``plot_top_residues``
     - ``10``
     - Number of residues in the hot-spot plots.
   * - ``energy_decomposition``
     - ``false``
     - A second, faster decomposition inside the main calculation
       (``decomposition_method: fast`` or ``full``), mainly for development.

Entropy
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Option
     - Default
     - Meaning
   * - ``run_entropy_analysis``
     - ``false``
     - Estimate the entropy term :math:`-T\Delta S` and report
       :math:`\Delta G = \Delta H - T\Delta S`. Setting ``entropy_method``
       also turns the estimate on.
   * - ``entropy_method``
     - ``interaction``
     - ``interaction`` (interaction entropy) or ``c2`` (its second-order
       cumulant), both from the fluctuation of the interaction energy and
       reported with its standard deviation and reliability.
       ``quasiharmonic``: quasi-harmonic entropy of the ligand
       (``qha_analyze_complex: true`` adds receptor and complex to give a
       binding entropy). ``normal_mode``: normal-mode entropy of the ligand
       alone (needs ``ligand_mol``), reported but not added to
       :math:`\Delta G`.
       See :doc:`../theory/entropy_methods`.
   * - ``qha_analyze_complex``
     - ``false``
     - With ``quasiharmonic``: also analyse the receptor and the complex.
   * - ``nmode_cutoff_frequency``
     - ``10.0``
     - With ``normal_mode``: modes below this frequency (cm⁻¹) are left out.

Execution
~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 22 54

   * - Option
     - Default
     - Meaning
   * - ``max_workers``
     - min(cores - 1, 4)
     - Worker processes on the CPU platform. Each evaluates its own frames.
   * - ``print_interval``
     - ``10``
     - Print progress every N frames.
   * - ``use_cache``
     - ``true``
     - Reuse the parameterised systems within a run.

forcefield_settings
-------------------

.. list-table::
   :header-rows: 1
   :widths: 24 14 62

   * - Option
     - Default
     - Meaning
   * - ``protein_forcefield``
     - ``amber``
     - Coordinate Mode only. ``amber``: Amber ff14SB. ``amber14``: OpenMM's
       ``amber14-all`` (proteins, DNA, RNA). ``charmm`` or ``charmm_gromacs``:
       CHARMM36.
   * - ``charmm_params``
     - none
     - ``.psf`` input: the CHARMM parameter files (``.prm``, ``.str``,
       ``.rtf``), one path or a list.
   * - ``charmm_coordinates``
     - none
     - ``.psf`` input: coordinates (``.pdb``, ``.coor``, ``.crd``) when no
       ``solvated_topology`` is given.

platform_settings
-----------------

.. list-table::
   :header-rows: 1
   :widths: 24 14 62

   * - Option
     - Default
     - Meaning
   * - ``preferred_platform``
     - automatic
     - ``CUDA``, ``OpenCL``, ``CPU`` or ``Reference``. Automatic tries CUDA,
       then OpenCL, then CPU, unless the environment variable
       ``OPENMM_DEFAULT_PLATFORM`` names one.
   * - ``decomposition_platform``
     - as above
     - Platform of the per-residue decomposition. With ``CPU`` it runs on
       several processes (``advanced_settings.n_jobs``).
   * - ``strict_platform``
     - ``false``
     - ``true``: stop with an error when the requested platform cannot be
       used, instead of falling back to the next one.

GPU runs use mixed precision on the first visible device; choose the device
with ``CUDA_VISIBLE_DEVICES``.

advanced_settings
-----------------

.. list-table::
   :header-rows: 1
   :widths: 24 14 62

   * - Option
     - Default
     - Meaning
   * - ``n_jobs``
     - ``-1``
     - Worker processes of a per-residue decomposition on the CPU platform;
       -1 uses every core.

reporting_settings
------------------

These options shape the interactive HTML report, which the per-residue
decomposition writes.

.. list-table::
   :header-rows: 1
   :widths: 24 14 62

   * - Option
     - Default
     - Meaning
   * - ``include_all_plots``
     - ``true``
     - Include the time series and distribution plots.
   * - ``export_format``
     - ``png``
     - Format of the plot export button (any Plotly format).
   * - ``export_width``, ``export_height``
     - ``4200``, ``2800``
     - Size of exported plots in pixels.

Validation
----------

``opengbsa`` checks the configuration before the run. It stops with a message
when a required file is missing, when ``gb_model``, ``gb_radii_set``,
``sa_model``, ``lcpo_typing``, ``binding_mode``, ``frame_selection``,
``entropy_method`` or ``decomposition_method`` has an unknown value, when
``temperature``, ``salt_concentration`` or ``max_frames`` is out of range, when
``frame_end`` is not after ``frame_start``, or when ``decomp_frames`` exceeds
``max_frames``. Unknown options in other sections are ignored, so check the
spelling of an option that seems to have no effect.
