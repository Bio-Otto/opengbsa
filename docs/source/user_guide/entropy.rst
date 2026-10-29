Entropy
=======

OpenGBSA reports the binding enthalpy :math:`\Delta H` by default. An estimate
of the entropy term :math:`-T\Delta S` is added on request, and the binding free
energy is then :math:`\Delta G = \Delta H - T\Delta S`. The entropy term is much
less precise than the enthalpy, so it is off by default, as in Amber MMPBSA.py.
This page explains the options. The methods and their validation are described
in :doc:`../theory/entropy_methods`.

.. code-block:: yaml

    analysis_settings:
      run_entropy_analysis: true
      entropy_method: interaction   # interaction, c2, quasiharmonic or nmode
      temperature: 300.0

.. _opt-run_entropy_analysis:

``run_entropy_analysis``
------------------------

:Default: ``false``
:Values: ``true``, ``false``

Switches the entropy estimate on. With ``false`` no entropy is computed and
``entropy_method`` is ignored.

.. _opt-entropy_method:

``entropy_method``
------------------

:Default: ``interaction``
:Values: ``interaction``, ``c2``, ``quasiharmonic``, ``nmode`` (also written ``normal_mode``), ``none``

.. list-table::
   :header-rows: 1
   :widths: 18 30 26 26

   * - Method
     - What it uses
     - Cost
     - Use it when
   * - ``interaction``
     - the fluctuation of the receptor-ligand interaction energy over the
       analysed frames
     - none beyond the enthalpy
     - the standard deviation :math:`\sigma` of the interaction energy is
       below 3.6 kcal/mol
   * - ``c2``
     - the same fluctuation, to second order
     - none beyond the enthalpy
     - :math:`\sigma` is below 6.0 kcal/mol
   * - ``quasiharmonic``
     - the covariance of the atomic fluctuations over the analysed frames
     - seconds to minutes
     - the trajectory samples the motions of the ligand well
   * - ``nmode``
     - minimized structures of the complex, the receptor and the ligand
     - minutes per frame on a GPU, hours for large complexes
     - a comparison with MMPBSA.py or gmx_MMPBSA normal-mode results is needed

*How it is applied.*

* ``interaction`` and ``c2`` use the gas-phase interaction energy
  :math:`\Delta E_\mathrm{vdW} + \Delta E_\mathrm{elec}` of every analysed frame.
  Both values are reported with :math:`\sigma`, whichever is selected, and an
  estimate above its reliability limit is marked unreliable. For most
  protein-protein and protein-nucleic acid complexes :math:`\sigma` is far above
  both limits.
* ``quasiharmonic`` superimposes the analysed frames and computes the entropy
  of the ligand from its covariance matrix. With ``qha_analyze_complex: true``
  the receptor and the complex are analysed as well, which gives a binding
  entropy.
* ``nmode`` takes ``nmode_frames`` of the analysed frames, minimizes the
  complex, the receptor and the ligand in each, computes their vibrational
  frequencies from a finite-difference Hessian and gives the binding entropy of
  each frame. The mean over frames and its standard deviation are reported.

``none`` is the same as ``run_entropy_analysis: false``.

.. _opt-qha_analyze_complex:

``qha_analyze_complex``
-----------------------

:Default: ``false``
:Values: ``true``, ``false``

With ``entropy_method: quasiharmonic``, also analyse the receptor and the
complex, so that :math:`T\Delta S = TS_\mathrm{complex} - TS_\mathrm{receptor} -
TS_\mathrm{ligand}` is a binding entropy. Without it only the ligand's own
entropy is computed. The covariance matrix of the receptor and the complex has
three rows per atom, so this needs memory and time for large receptors.

Normal-mode options
-------------------

These options apply only with ``entropy_method: nmode``.

.. _opt-nmode_frames:

``nmode_frames``
~~~~~~~~~~~~~~~~

:Default: ``10``
:Values: a positive number of frames, or ``null`` for every analysed frame
:Amber: MMPBSA.py's ``&nmode`` ``nmstartframe``, ``nmendframe`` and ``nminterval``

The number of analysed frames used for the normal-mode estimate, spread evenly
over them: the first and the last analysed frames are always included, and
the others lie in between at equal spacing. Each frame takes minutes on a GPU,
so a few frames are usual.

.. _opt-nmode_drms:

``nmode_drms``
~~~~~~~~~~~~~~

:Default: ``0.001``
:Values: a positive number, in kcal/mol/Å
:Amber: ``drms``

Each minimization runs until the root-mean-square gradient is below this
value. A frame in which the complex, the receptor or the ligand does not reach
it is left out of the mean and reported in the log and in
``nmode_entropy.csv``.

.. _opt-nmode_maxcyc:

``nmode_maxcyc``
~~~~~~~~~~~~~~~~

:Default: ``10000``
:Values: a positive number of iterations
:Amber: ``maxcyc``

The maximum number of iterations of each minimization round.

.. _opt-nmode_precision:

``nmode_precision``
~~~~~~~~~~~~~~~~~~~

:Default: ``mixed``
:Values: ``mixed``, ``double``

The floating-point precision of the minimization and the Hessian on CUDA and
OpenCL. ``mixed`` minimizes first in mixed precision, finishes in double
precision and computes the Hessian in mixed precision. ``double`` uses double
precision throughout. CPU runs always use double precision. At the same minima
of a 4,712-atom complex, ``double`` reproduced Amber's entropies within 0.02
cal/(mol K) and ``mixed`` gave :math:`-T\Delta S` 0.16 kcal/mol higher
(:doc:`../theory/entropy_methods`). On consumer GPUs, whose double-precision
speed is a small fraction of their mixed-precision speed, ``double`` took three
to four times as long: 62 to 66 instead of 16 to 19 minutes per frame for that
complex on an NVIDIA GeForce RTX 2080 SUPER.

The temperature
---------------

``temperature`` (:ref:`opt-temperature`) is the temperature of every estimate:
:math:`k_\mathrm{B}T` in the interaction entropy and C2, and :math:`T` in
:math:`T\Delta S`. With ``amber_compatible: true`` the normal-mode estimate uses
MMPBSA.py's energy function and its fixed temperature of 298.15 K, so that the
two programs can be compared directly.
