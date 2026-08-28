Protein-RNA: Amber-Native Topology Split
===========================================

This example shows a non-small-molecule ligand (an RNA chain) handled
through Native Mode's explicit ``receptor_topology``/``ligand_topology``
inputs, rather than through ``binding_mode: ppi`` chain selection.

Dataset
-------

A 288-residue protein bound to a 25-nucleotide RNA chain (system
**3dd2**), from Zenodo
`6973437 <https://doi.org/10.5281/zenodo.6973437>`_ (CC-BY-4.0),
"Investigating RNA-Protein Recognition Mechanisms through Supervised
Molecular Dynamics (SuMD) Simulations". Ships a real Amber ``.prmtop``,
already GB-radii-parameterized -- no topology reconstruction needed.

Run it
------

.. code-block:: bash

    cd test/configs/rna_protein_test
    ./fetch_data.sh          # ~44 MB
    python3 prepare_system.py
    opengbsa 3dd2_config.yaml

For a quick 3-frame smoke test instead of the full 300-frame run, use
``3dd2_config_test.yaml``.

What to notice in the config
-----------------------------

.. code-block:: yaml

    input_files:
      complex_pdb: raw/3dd2.prmtop
      trajectory: raw/3dd2_suMD1.dcd
      receptor_topology: prepared/protein.prmtop
      ligand_topology: prepared/ligand.prmtop

    forcefield_settings:
      protein_forcefield: amber

- ``prepare_system.py`` here does only one thing: split the single
  complex ``.prmtop`` into ``protein.prmtop`` (atoms 0-4627) and
  ``ligand.prmtop`` (atoms 4628-5433) via ParmEd's ``Structure.strip()``.
  Amber's residue-template force fields already recognize RNA residue
  names (``A``/``U``/``G``/``C`` and terminal variants) natively -- no
  extra configuration needed for the nucleic acid itself.
- This is the Native Mode equivalent of the DNA example's Coordinate Mode
  approach: same "ligand is a nucleic acid chain" problem, solved via an
  explicit topology split here instead of ``chainid`` selection, because
  a real topology was available to split.

Full write-up (including two OpenGBSA bugs found in nucleic-acid ligand
detection during this validation, both in how Native Mode identifies
which atoms are "the ligand" when a topology is already provided): see
``test/configs/rna_protein_test/README.md`` in the repository.
