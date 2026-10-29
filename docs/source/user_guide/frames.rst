Choosing frames
===============

MM/GBSA averages the binding enthalpy over frames of the trajectory. The
options on this page decide which frames are used and how they are prepared.
They belong to ``analysis_settings``.

.. code-block:: yaml

    analysis_settings:
      frame_start: 0                # first frame of the range, counted from 0
      frame_end: null               # end of the range (excluded); null = last frame
      frame_stride: 1               # every N-th frame of the range
      max_frames: 50                # at most this many frames; null = all
      frame_selection: sequential   # sequential, equidistant or random
      random_seed: 42               # seed of the random selection
      reimage_trajectory: true      # make molecules whole, one periodic image

How the frames are chosen
-------------------------

OpenGBSA chooses the frames in four steps:

#. The range: frames ``frame_start`` up to, but not including, ``frame_end``.
   A value outside the trajectory is moved to its first or last frame, with a
   warning.
#. The stride: every ``frame_stride``-th frame of the range, starting with the
   first.
#. The selection method picks ``max_frames`` of these frames.
#. The chosen frames are evaluated in the order of the trajectory.

For example, a trajectory of 1,000 frames with ``frame_start: 500``,
``frame_stride: 2`` and ``max_frames: 100`` leaves 250 frames after the
stride. ``sequential`` takes the frames 500, 502, …, 698. ``equidistant`` takes
every 2.5th of the 250 frames, rounded down: 500, 504, 510, 514, …, 994. The
column ``trajectory_frame`` of the per-frame output gives the index of each
evaluated frame.

.. _opt-frame_start:

``frame_start`` and ``frame_end``
---------------------------------

:Default: the first and the last frame
:Values: whole numbers, counted from 0; ``frame_end`` larger than ``frame_start``
:Amber: MMPBSA.py's ``startframe`` and ``endframe``, which count from 1 and include ``endframe``

The range of frames to consider. Leave out the equilibration part of the
simulation here. MMPBSA.py's ``startframe = 501, endframe = 1000`` corresponds
to ``frame_start: 500, frame_end: 1000``.

.. _opt-frame_stride:

``frame_stride``
----------------

:Default: ``1``
:Values: a positive whole number
:Amber: ``interval``

Use every N-th frame of the range. Frames close in time are correlated, so a
stride makes the frames more independent without changing the time span they
cover.

.. _opt-max_frames:

``max_frames``
--------------

:Default: ``50``
:Values: 1 to 10,000, or ``null`` for every frame of the range after the stride

The number of frames to evaluate. The standard error of the mean binding
enthalpy falls with the square root of the number of independent frames: the
quickstart's standard deviation of 4.60 kcal/mol gives a standard error of 0.65
kcal/mol with 50 frames and 1.46 kcal/mol with 10. On a GPU a frame takes
milliseconds, so the cost of more frames is small. When frames are correlated,
the standard error is underestimated; a stride or an equidistant selection
helps.

.. _opt-frame_selection:

``frame_selection``
-------------------

:Default: ``sequential``
:Values: ``sequential``, ``equidistant``, ``random``

* ``sequential``: the first ``max_frames`` frames after the stride.
* ``equidistant``: ``max_frames`` frames spread evenly over all frames after
  the stride, so that they cover the whole range. Of :math:`N` frames after the
  stride, the frames at positions :math:`\lfloor i N / \mathrm{max\_frames}
  \rfloor`, :math:`i = 0, 1, \ldots`, are taken; the first frame is always
  included.
* ``random``: ``max_frames`` different frames drawn at random, with
  ``random_seed``, from the frames after the stride.

``equidistant`` is the better choice when ``max_frames`` is smaller than the
number of frames in the range: ``sequential`` then covers only the beginning of
the range.

.. _opt-random_seed:

``random_seed``
---------------

:Default: ``42``
:Values: a whole number

The seed of ``frame_selection: random``. The same seed gives the same frames,
so a run can be repeated exactly.

.. _opt-reimage_trajectory:

``reimage_trajectory``
----------------------

:Default: ``true``
:Values: ``true``, ``false``

A simulation with periodic boundaries can write a molecule broken across the
box, or receptor and ligand in different periodic images. Their energies would
then be wrong by hundreds of kcal/mol. With ``true``, OpenGBSA repairs every
frame before computing energies:

#. Each molecule is made whole by following its bonds.
#. The largest molecule and every molecule at least half its size, for example
   the chains of a protein complex, are brought into one periodic image, largest
   first.
#. Every other molecule is placed in the image nearest to them.

Rectangular and triclinic boxes, such as GROMACS' rhombic dodecahedron, are
both handled, with the box vectors of each frame: for a triclinic box OpenGBSA
searches the 27 neighbouring images for the nearest one. A trajectory without
box information is used as it is, with a warning.

*When to change it.* Set ``false`` for a trajectory that was already imaged,
for example with ``gmx trjconv -pbc cluster`` or ``cpptraj autoimage``, to save
time.

Frames of the decomposition and the normal-mode entropy
-------------------------------------------------------

The per-residue decomposition and the normal-mode entropy use a subset of the
analysed frames, ``decomp_frames`` and ``nmode_frames``. The subset is always
taken from the frames chosen above (:doc:`decomposition`, :doc:`entropy`).
