Protein-Ligand Batch: SARS-CoV-2 3CLpro
========================================

This example is the standard OpenGBSA case: a real Amber topology, a
dry (solvent-stripped) trajectory, and Native Mode -- no system
reconstruction needed. It also shows how to run a batch of related
ligands against the same receptor.

Dataset
-------

Seven independent complexes -- SARS-CoV-2 3CLpro bound to baicalein and
six 7-hydroxystilbene-coumarin hybrid inhibitors (compound7a-f) -- from
Zenodo `17926575 <https://doi.org/10.5281/zenodo.17926575>`_ (CC-BY-4.0).
Amber24 MD, 300 ns dry trajectories, 30000 frames each.

.. warning::
   The full dataset (all 7 ligands) ships as a single 9.1 GB archive with
   no selective per-file download available. Fetching even one ligand
   requires downloading the whole thing first. Budget real time for this
   one.

Run it
------

.. code-block:: bash

    cd test/configs/3clpro_zenodo_test
    ./fetch_data.sh          # downloads 9.1 GB, extracts baicalein only
    opengbsa baicalein_config.yaml

The other six ligands (``compound7a`` through ``compound7f``) follow the
identical pattern once extracted -- see the folder's own
``README.md`` for the exact extraction command.

What to notice in the config
-----------------------------

.. code-block:: yaml

    input_files:
      complex_pdb: raw/1_MD/baicalein/complex.prmtop
      trajectory: raw/1_MD/baicalein/dry_MD_baicalein.trj
      receptor_topology: raw/1_MD/baicalein/protein.prmtop
      ligand_topology: raw/1_MD/baicalein/ligand.prmtop

    analysis_settings:
      frame_start: 29700   # 0-indexed; last 300 of 30000 frames
      frame_end: 30000

- ``receptor_topology``/``ligand_topology`` are given explicitly, split
  from the same complex by the dataset's own authors -- OpenGBSA doesn't
  need to guess which atoms are the ligand.
- ``frame_start``/``frame_end`` select a window out of a much longer
  trajectory, the standard way to analyze only the equilibrated tail of a
  long production run.

Expected result
---------------

From the root of the repository, ``check_validation.py`` runs the
configuration and compares the mean binding enthalpy (kcal/mol) with the
recorded value; ``--amber-compatible`` runs it with ``amber_compatible: true``
and compares with Amber MMPBSA.py on the same frames:

.. code-block:: bash

    python test/configs/check_validation.py 3clpro_baicalein 3clpro_compound7a 3clpro_compound7b 3clpro_compound7c 3clpro_compound7d 3clpro_compound7e 3clpro_compound7f
    python test/configs/check_validation.py 3clpro_baicalein 3clpro_compound7a 3clpro_compound7b 3clpro_compound7c 3clpro_compound7d 3clpro_compound7e 3clpro_compound7f --amber-compatible

.. list-table::
   :header-rows: 1
   :widths: 22 26 26 26

   * - Case
     - OpenGBSA, default settings
     - OpenGBSA, ``amber_compatible``
     - Amber MMPBSA.py
   * - ``3clpro_baicalein``
     - -23.873
     - -23.888
     - -23.888
   * - ``3clpro_compound7a``
     - -22.034
     - -22.137
     - -22.138
   * - ``3clpro_compound7b``
     - -30.642
     - -30.699
     - -30.700
   * - ``3clpro_compound7c``
     - -40.916
     - -40.886
     - -40.886
   * - ``3clpro_compound7d``
     - -37.840
     - -37.879
     - -37.879
   * - ``3clpro_compound7e``
     - -27.277
     - -27.472
     - -27.472
   * - ``3clpro_compound7f``
     - -46.463
     - -46.399
     - -46.398

With ``amber_compatible`` OpenGBSA and MMPBSA.py agree within 0.01 kcal/mol.
At the default settings they differ by the Born-radius cutoff of MMPBSA.py
(:doc:`../theory/mmgbsa_theory`). Provenance and preparation
details: ``test/configs/3clpro_zenodo_test/README.md``.
