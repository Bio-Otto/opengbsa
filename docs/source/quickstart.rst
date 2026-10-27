Quickstart
==========

This page runs a complete calculation on a protein-ligand complex that ships
with the repository: the OXA-family beta-lactamase 6T1H with the inhibitor
1151. It takes a few minutes on a GPU and needs no download besides the
repository itself.

1. Get the example
------------------

Install OpenGBSA (:doc:`installation`) and clone the repository, which holds the
example files:

.. code-block:: bash

    conda activate opengbsa
    git clone https://github.com/Bio-Otto/opengbsa.git
    cd opengbsa

The example consists of three Amber topologies (a CHARMM36 system written as
Amber chamber files) and a trajectory of 50 frames, one per nanosecond from the
last 50 ns of a 100 ns simulation, with water and ions removed:

.. code-block:: text

    test/data/6t1h_1151_comp/complex.prmtop     receptor + ligand, 3,917 atoms
    test/data/6t1h_1151_comp/receptor.prmtop    receptor, 3,840 atoms
    test/data/6t1h_1151_comp/ligand.prmtop      ligand 1151, 77 atoms
    examples/quickstart/trajectory.xtc          50 frames
    examples/quickstart/quickstart.yaml         the configuration

2. The configuration
--------------------

The whole calculation is defined in one YAML file. This is the
configuration-driven interface (CDI, see :doc:`user_guide/interfaces`):

.. literalinclude:: ../../examples/quickstart/quickstart.yaml
   :language: yaml
   :lines: 12-

``input_files`` names the topology and the trajectory. Because receptor and
ligand topologies are given, their atom counts define which atoms belong to the
ligand. ``analysis_settings`` selects the OBC2 GB model at 0.15 M salt and
300 K, evaluates all 50 frames, decomposes the binding enthalpy into residue
contributions and estimates the entropy. Every option is described in the
:doc:`user_guide/configuration`.

3. Run it
---------

From the root of the repository:

.. code-block:: bash

    opengbsa examples/quickstart/quickstart.yaml

The log names each step. The first lines show the platform that was chosen,
for example ``Running on OpenMM platform: CUDA``. On an NVIDIA RTX 2080 SUPER
the whole run takes about 4.5 minutes, most of it in the per-residue
decomposition. Without a GPU it runs on the CPU and takes longer.

4. The result
-------------

The run ends with a summary like this:

.. code-block:: text

    Mean binding energy (Enthalpy): -24.99 ± 4.60 kcal/mol (Standard Deviation)
    Standard Error of Mean: 0.65 kcal/mol
    sigma(dE_int) = 8.79 kcal/mol over 50 frames at 300.0 K
    Interaction entropy (-TdS) = 19.28 kcal/mol  (unreliable: sigma > 3.6 kcal/mol, Ekberg & Ryde 2021)
    C2 entropy (-TdS) = 64.86 kcal/mol  (unreliable: sigma > 6.0 kcal/mol, Ekberg & Ryde 2021)

The mean binding enthalpy is -24.99 kcal/mol. The standard deviation of 4.60
kcal/mol describes how much it varies from frame to frame; the standard error
of the mean is 0.65 kcal/mol. The entropy estimates are flagged as unreliable
because the interaction energy fluctuates too much. :doc:`theory/entropy_methods`
explains why. Small differences in the last digit between GPU and CPU runs are
normal.

All files are written to ``quickstart_results/analysis_<date>_<time>/``
(:doc:`user_guide/outputs` describes each one). The most useful ones are:

* ``fixed_enhanced_mmgbsa_results_obc2.csv``: every energy term for every frame.
* ``final_report.txt`` and ``results_summary.yaml``: the summary in text and YAML.
* ``energy_analysis.png``: the binding enthalpy over the frames, its
  distribution, its running average and its terms (below).
* ``per_residue_detailed.csv`` and ``per_residue_decomposition.png``: the
  contribution of every receptor residue.
* ``interactive_report.html``: all results in one page that opens in a browser.
* ``view_binding.pml``: a PyMOL script that colours the receptor by residue
  contribution.

.. figure:: _static/figures/quickstart_energy.png
   :width: 100%
   :alt: Binding enthalpy over the 50 frames, its distribution, running average and terms

   ``energy_analysis.png`` of the quickstart. Top left: the binding enthalpy of
   each frame. Top right: its distribution. Bottom left: the running average,
   which flattens out when the mean has converged. Bottom right: the gas-phase
   interaction energy (van der Waals plus electrostatic), the polar solvation
   (GB) and the non-polar solvation (SA) terms of each frame.

.. figure:: _static/figures/quickstart_residues.png
   :width: 100%
   :alt: Per-residue decomposition plots of the quickstart

   ``per_residue_decomposition.png``. The receptor residues that contribute
   most to binding (top left), their van der Waals, electrostatic and
   solvation terms (top right), the distribution over all residues (bottom
   left) and van der Waals against electrostatic contributions (bottom right).
   PHE111 and ALA126 contribute most.

.. figure:: _static/figures/quickstart_pymol.png
   :width: 80%
   :align: center
   :alt: PyMOL view of the binding pocket coloured by residue contribution

   ``view_binding.pml`` opened in PyMOL (``cd`` into the results folder, then
   ``pymol view_binding.pml``). The ligand is green. Receptor residues are
   coloured from blue (favourable) through white to red (unfavourable), and
   residues that contribute less than -0.5 kcal/mol are shown as sticks with
   their contribution.

5. Next steps
-------------

* Compare the result with Amber MMPBSA.py conventions: add
  ``amber_compatible: true`` under ``analysis_settings`` and run again
  (:doc:`theory/mmgbsa_theory` explains the difference).
* Run the unit tests and the third-party validation cases: :doc:`testing`.
* Set up your own system: :doc:`user_guide/input_formats` and
  :doc:`user_guide/configuration`.
