Protein-Peptide: NAMD/CHARMM Native Mode
==========================================

This example shows Native Mode for CHARMM-family input: a NAMD ``.psf``
topology and CHARMM36 parameter files, loaded directly with no OpenMM
``ForceField`` parameterization step -- and ``binding_mode: ppi`` for a
receptor/ligand split where the "ligand" is a peptide chain, not a small
molecule.

Dataset
-------

A CTLA4-derived peptide bound to its protein partner, from Zenodo
`7186684 <https://doi.org/10.5281/zenodo.7186684>`_ (CC-BY-4.0). Native
NAMD production run, ships as a 321 MB archive (only ~85 MB of it is
actually needed).

Run it
------

.. code-block:: bash

    cd test/configs/namd_charmm_test
    ./fetch_data.sh
    opengbsa ctla4_peptide_config.yaml

For a quick 3-frame smoke test instead of the full 300-frame run, use
``ctla4_peptide_config_test.yaml``.

What to notice in the config
-----------------------------

.. code-block:: yaml

    input_files:
      complex_pdb: raw/Zenodo_CTLA4_Peptide/ctla4_P16_wat.psf
      solvated_topology: raw/Zenodo_CTLA4_Peptide/ctla4_P16_wat.pdb
      trajectory: raw/Zenodo_CTLA4_Peptide/output/stride50_pro_ctla4_P16.D.dcd

    analysis_settings:
      binding_mode: ppi
      receptor_selection: 'chainid 1'
      ligand_selection: 'chainid 0'

    forcefield_settings:
      charmm_params:
        - raw/Zenodo_CTLA4_Peptide/par_all36m_prot.prm
        - raw/Zenodo_CTLA4_Peptide/toppar_water_ions_prot.str
      charmm_coordinates: raw/Zenodo_CTLA4_Peptide/ctla4_P16_wat.pdb

- A ``.psf`` as ``complex_pdb`` (Native Mode's entry point for CHARMM
  input, analogous to ``.prmtop`` for Amber or ``.tpr`` for GROMACS)
  triggers CHARMM parameterization instead of OpenMM's ``ForceField``.
  ``forcefield_settings.charmm_params``/``charmm_coordinates`` supply the
  parameter/stream files and reference coordinates CHARMM parameterization
  needs.
- ``receptor_selection``/``ligand_selection`` as ``chainid N`` strings
  (rather than a small-molecule ``ligand_resname``) is what makes
  ``binding_mode: ppi`` work for a peptide or protein "ligand" -- the same
  mechanism used in the two Coordinate Mode PPI examples.

Full write-up (including the independent Amber reference comparison and
six OpenGBSA bugs found and fixed while validating this path): see
``test/configs/namd_charmm_test/README.md`` in the repository.
