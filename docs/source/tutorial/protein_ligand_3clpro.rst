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

Full write-up (including the other six ligands' results and a
per-residue decomposition comparison): see
``test/configs/3clpro_zenodo_test/README.md`` in the repository.
