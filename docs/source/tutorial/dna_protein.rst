Protein-DNA: Coordinate Mode on a Nucleic Acid
=================================================

This example combines two things seen separately in earlier tutorials:
Coordinate Mode (no topology provided) and a nucleic-acid "ligand" -- this
time DNA rather than RNA, and split via ``binding_mode: ppi`` chain
selection rather than an explicit topology split.

Dataset
-------

A 109-residue protein (IRF11) bound to a 42-nucleotide double-stranded DNA
fragment, from Zenodo
`14377950 <https://doi.org/10.5281/zenodo.14377950>`_ (CC-BY-4.0),
"Molecular dynamics simulation of the IRF11-DNA complex". Delivered as a
single 101-frame multi-model PDB trajectory -- no separate topology,
coordinate, or trajectory files.

Run it
------

.. code-block:: bash

    cd test/configs/dna_protein_test
    ./fetch_data.sh          # ~24 MB
    python3 prepare_system.py
    opengbsa irf11_dna_config.yaml

For a quick 3-frame smoke test instead of the full 101-frame run, use
``irf11_dna_config_test.yaml``.

What to notice in the config
-----------------------------

.. code-block:: yaml

    input_files:
      complex_pdb: prepared/complex_ref.pdb
      trajectory: prepared/complex_traj.dcd

    analysis_settings:
      binding_mode: ppi
      receptor_selection: 'chainid 0'
      ligand_selection: 'chainid 1'

    forcefield_settings:
      protein_forcefield: amber14-all.xml

- ``prepare_system.py`` here only splits the raw multi-model PDB into a
  single reference structure (frame 0) plus a ``.dcd`` trajectory, and
  relabels histidine residues -- no bond-inference workaround was needed
  for this particular file, unlike the GROMACS-sourced PPI datasets.
- Amber's residue-template force fields already recognize DNA residue
  names (``DA``/``DC``/``DG``/``DT`` and terminal variants) natively.
- This is the same ``binding_mode: ppi`` + Coordinate Mode combination
  used in the protein-protein tutorial, applied here to confirm it
  generalizes to a fourth distinct receptor/ligand kind (nucleic acid
  double helix) without any further code changes.

Full write-up (component-by-component analysis of the OpenGBSA-vs-Amber
gap, and why it's a parameterization-sensitivity effect rather than a
bug): see ``test/configs/dna_protein_test/README.md`` in the repository.
