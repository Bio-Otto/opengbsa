A first example
===============

This page goes through the configuration of the :doc:`../quickstart` line by
line and then follows a run from start to end. After it, the other pages of the
user guide describe every option in detail.

The configuration
-----------------

.. literalinclude:: ../../../examples/quickstart/quickstart.yaml
   :language: yaml
   :lines: 12-

A configuration is a YAML file with a few sections. Each section holds options
written as ``name: value``. Indentation shows which section an option belongs
to. A ``#`` starts a comment.

.. code-block:: yaml

    input_files:
      complex_pdb: test/data/6t1h_1151_comp/complex.prmtop

``input_files`` lists the files OpenGBSA reads. ``complex_pdb`` is the system:
here an Amber topology of the receptor and the ligand together. Because it is a
topology, OpenGBSA runs in Native Mode and uses the force field of the
simulation as it is stored in the file (:doc:`input_formats`). Paths are
relative to the folder from which ``opengbsa`` is started.

.. code-block:: yaml

      receptor_topology: test/data/6t1h_1151_comp/receptor.prmtop
      ligand_topology: test/data/6t1h_1151_comp/ligand.prmtop

The topologies of the receptor and of the ligand alone. In the complex topology
the receptor atoms come first and the ligand atoms follow, so the atom count of
``receptor.prmtop`` tells OpenGBSA where the ligand starts. Without these two
files OpenGBSA finds the ligand by its residue name
(:ref:`receptor-and-ligand`).

.. code-block:: yaml

      trajectory: examples/quickstart/trajectory.xtc

The MD trajectory, here 50 frames from which water and ions were already
removed. Any format that MDTraj reads works.

.. code-block:: yaml

    output_settings:
      output_directory: quickstart_results

Every run writes a new folder ``analysis_<date>_<time>`` inside this
directory, so a second run never overwrites the first (:doc:`outputs`).

.. code-block:: yaml

    analysis_settings:
      gb_model: OBC2
      salt_concentration: 0.15
      temperature: 300.0

The energy function: the OBC2 Generalized Born model for the polar solvation
energy, with the salt screening of 0.15 M at 300 K. Every other option of the
energy function keeps its default (:doc:`energy_function`).

.. code-block:: yaml

      max_frames: 50

Evaluate up to 50 frames. With the defaults these are the first 50 frames of
the trajectory, which here is all of them (:doc:`frames`).

.. code-block:: yaml

      run_per_residue_decomposition: true
      decomp_frames: 50

Split the binding enthalpy into the contributions of the receptor residues,
using 50 of the analysed frames (:doc:`decomposition`).

.. code-block:: yaml

      run_entropy_analysis: true
      entropy_method: interaction

Estimate the entropy term with the interaction entropy method
(:doc:`entropy`).

Everything not written in the file takes its default. Options that OpenGBSA
does not know are reported as warnings, with the closest known name, so a typing
error does not go unnoticed (:doc:`configuration`).

What happens during a run
-------------------------

``opengbsa examples/quickstart/quickstart.yaml`` goes through these steps. The
log announces each of them.

#. **Validation.** The configuration is checked before anything is computed:
   every file must exist, and every option must have an allowed value. An error
   stops the run with a message that names the option. Unknown options give a
   warning.
#. **Topology.** The complex topology is read. In Native Mode the force-field
   parameters come from it; in Coordinate Mode they come from the force field
   named in ``forcefield_settings``.
#. **Trajectory.** The trajectory is read with MDTraj. Water and ions are
   removed if it still contains them, and every frame is re-imaged so that each
   molecule is whole and receptor and ligand are in the same periodic image
   (:ref:`opt-reimage_trajectory`).
#. **Frames.** The frames to evaluate are chosen from the range, the stride
   and the selection method (:doc:`frames`).
#. **Receptor and ligand.** The atoms of the complex are divided into receptor
   and ligand.
#. **Systems.** OpenGBSA builds three OpenMM systems, for the complex, the
   receptor and the ligand, each with the Generalized Born force and the
   non-polar force described in :doc:`energy_function`.
#. **Energies.** For every frame the coordinates of the complex are set, and
   the receptor and the ligand take the coordinates of their atoms from the same
   frame. Each energy term is evaluated separately, and the binding enthalpy of
   the frame is complex minus receptor minus ligand.
#. **Statistics.** The mean binding enthalpy, its standard deviation over the
   frames and the standard error of the mean are computed.
#. **Entropy and decomposition**, when requested.
#. **Output.** Tables, plots, a report and the configuration as it was run are
   written to the results folder (:doc:`outputs`).

The run ends with a summary:

.. code-block:: text

    Mean binding energy (Enthalpy): -24.99 ± 4.60 kcal/mol (Standard Deviation)

.. TODO re-check this line against the 0.0.23 quickstart output
