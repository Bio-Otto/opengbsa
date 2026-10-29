Per-residue decomposition
=========================

The per-residue decomposition splits the binding enthalpy into the
contributions of the receptor residues. It shows which residues hold the ligand
and is used to find hot spots, residues whose mutation would weaken binding
most. It belongs to ``analysis_settings``.

.. code-block:: yaml

    analysis_settings:
      run_per_residue_decomposition: true
      decomp_frames: 50            # a subset of the analysed frames
      plot_top_residues: 10        # residues in the plots, tables and report
      report_raw_energies: false   # also write complex, receptor and ligand energies
    platform_settings:
      decomposition_platform: CUDA

How the contributions are computed
----------------------------------

For every decomposed frame the energies of the complex, the receptor and the
ligand are divided among their atoms, and a residue's contribution is the
change of its atoms' energies between the complex and the free receptor:

* **Van der Waals and electrostatic:** half of each pair energy between an atom
  of the residue and an atom of the ligand, the other half going to the ligand,
  as in the scheme of Gohlke et al. (2003) that Amber uses. The parameters are
  those of the run's own OpenMM system, so CHARMM's pair-specific
  Lennard-Jones terms and ``solute_dielectric`` are included.
* **Polar solvation:** each atom's GB self term plus half of each of its GB pair
  terms, computed from the run's own GB force for every GB model, with the
  run's Born radii, ``rgbmax``, salt screening and dielectric constants.
* **Non-polar solvation:** by default each atom's part of the run's own LCPO
  (or ACE) energy. With ``amber_compatible: true`` OpenGBSA follows MMPBSA.py,
  whose decomposition uses a numerical surface instead of LCPO: the change of
  each atom's Shrake-Rupley surface area times 0.0072 kcal/mol/Å². As in Amber,
  these terms then do not add up to the LCPO term of the binding enthalpy.

The receptor residues together with the ligand's share give the binding
enthalpy of each frame term by term. OpenGBSA checks this on every frame and
writes the result to ``decomposition_energy_check.csv`` (:doc:`outputs`). The
sum over the receptor residues alone is smaller, because the ligand keeps half
of every pair term and its own change of solvation.

*Validation.* For all five GB models, with and without ``rgbmax`` and salt, and
for LCPO and ACE, the per-atom terms add up to OpenMM's energies of the complex,
the receptor and the ligand on the Reference platform to within 10⁻⁹ kcal/mol.
``test/unit/test_decomposition_exactness.py`` checks this with a tolerance of
10⁻³ kcal/mol.
Against Amber MMPBSA.py's decomposition (``idecomp=2``) on 3CLpro with
baicalein, 306 receptor residues over three frames, the van der Waals,
electrostatic and polar terms of every model agree with a Pearson
:math:`r` of 1.000 and differences of at most 0.008 kcal/mol, and with
``amber_compatible`` the totals with :math:`r \geq 0.999998`. On 18 runs of nine
OXA-48 complexes, 4,320 receptor residues averaged over 200 frames each,
OpenGBSA with ``amber_compatible`` and MMPBSA.py's ``idecomp=2`` agree with
:math:`r > 0.9999` for every term. The per-residue totals differ by a median of
3.9 × 10⁻⁵ kcal/mol and by at most 0.0034 kcal/mol.

.. _opt-run_per_residue_decomposition:

``run_per_residue_decomposition``
---------------------------------

:Default: ``false``
:Values: ``true``, ``false``
:Amber: MMPBSA.py's ``&decomp`` namelist

Switches the decomposition on. It writes the per-residue tables, the hot-spot
plots, the interactive HTML report and the PyMOL view (:doc:`outputs`,
:doc:`visualization`). A failed decomposition is reported as failed in the log,
``results_summary.yaml`` and ``final_report.txt``; the binding enthalpy of the
run is kept. With ``strict_platform: true`` the run then stops with an error.

.. _opt-decomp_frames:

``decomp_frames``
-----------------

:Default: ``10``
:Values: a positive whole number, at most ``max_frames``, or ``null`` for every analysed frame

The number of frames decomposed. They are always a subset of the analysed
frames, chosen like ``frame_selection``: the first ones with ``sequential``,
evenly spread with ``equidistant`` and a random set with ``random``
(:doc:`frames`). The decomposition takes longer per frame than the binding
enthalpy, because every receptor residue is evaluated, so fewer frames than
``max_frames`` are often used.

.. _opt-plot_top_residues:

``plot_top_residues``
---------------------

:Default: ``10``
:Values: a positive whole number

The number of residues in the hot-spot plots, ``binding_hot_spots.csv`` and the
tables and heat map of the HTML report.

.. _opt-report_raw_energies:

``report_raw_energies``
-----------------------

:Default: ``false``
:Values: ``true``, ``false``

Also write the energies of each residue in the complex, the receptor and the
ligand (``per_residue_complex.csv``, ``per_residue_receptor.csv`` and
``per_residue_ligand.csv``), from which the binding contributions are the
differences. The contributions themselves do not change.

.. _opt-energy_decomposition:

``energy_decomposition`` and ``decomposition_method``
-----------------------------------------------------

:Default: ``false`` and ``full``
:Values: ``true``, ``false``; ``full``, ``fast``

A second decomposition inside the main calculation, mainly for development.
``fast`` gives only the gas-phase van der Waals and electrostatic terms per
residue (``final_decomposition_fast.csv``), with a PyMOL script and plots;
``full`` writes ``final_decomposition.csv`` over all analysed frames. For
results use ``run_per_residue_decomposition``.

Platforms
---------

The per-atom energies are computed with NumPy from the parameters of the run's
forces. ``platform_settings.decomposition_platform``
(:ref:`opt-decomposition_platform`) sets the OpenMM platform that computes the
complex, receptor and ligand energies of every decomposed frame for the check
in ``decomposition_energy_check.csv``. With ``CPU`` the frames are divided among
``advanced_settings.n_jobs`` processes (:ref:`opt-n_jobs`); with a GPU platform
they are decomposed one after another. Serial and parallel runs give the same
contributions.

On CUDA and OpenCL the check computes in double precision, whatever precision
the run uses. In mixed precision a GPU's GB energy of a large complex differs
from the exact sum by up to about 10⁻⁵ of its value: 0.046 kcal/mol for the
3,400 kcal/mol of an antibody Fv bound to PD-1, against 4 × 10⁻⁴ kcal/mol in
double precision. A frame fails the check when a per-atom sum differs from
OpenMM's energy by more than 0.01 kcal/mol plus 10⁻⁵ of the energy, or 10⁻⁴ on
a GPU without double precision.
