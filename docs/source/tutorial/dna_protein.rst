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

Expected result
---------------

From the root of the repository, ``check_validation.py`` runs the
configuration and compares the mean binding enthalpy (kcal/mol) with the
recorded value; ``--amber-compatible`` runs it with ``amber_compatible: true``
and compares with Amber MMPBSA.py on the same frames:

.. code-block:: bash

    python test/configs/check_validation.py irf11
    python test/configs/check_validation.py irf11 --amber-compatible

.. list-table::
   :header-rows: 1
   :widths: 22 26 26 26

   * - Case
     - OpenGBSA, default settings
     - OpenGBSA, ``amber_compatible``
     - Amber MMPBSA.py
   * - ``irf11``
     - -103.845
     - -109.647
     - -109.623

For this complex OpenGBSA and MMPBSA.py differ by -0.024 kcal/mol even with
``amber_compatible``. OpenMM's Coulomb constant is 3.5e-5 larger than
Amber's. The electrostatic term of this complex (-7,225 kcal/mol) therefore
differs by -0.25 kcal/mol and the polar term (+7,220 kcal/mol) by +0.22
kcal/mol, and the two largely cancel. At the default settings the two programs
also differ by the Born-radius cutoff of MMPBSA.py (:doc:`../theory/mmgbsa_theory`).
Provenance and preparation
details: ``test/configs/dna_protein_test/README.md``.
