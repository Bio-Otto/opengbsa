Protein-Protein: GROMACS Coordinate Mode
===========================================

This example shows Coordinate Mode -- parameterizing a system from plain
coordinates via OpenMM's ``ForceField`` -- applied to a protein-protein
complex where the source dataset provides no topology at all, only
``.gro``/``.xtc``.

Dataset
-------

Two independent complexes from the same Zenodo record
`6638504 <https://doi.org/10.5281/zenodo.6638504>`_ (CC-BY-4.0), "A
dynamical view of protein-protein complexes":

- **1GCQ** (Vav/GRB2 SH3 domain complex), 1984 atoms
- **2OOB** (ubiquitin/ubiquitin-ligase complex), 1895 atoms

Neither ships a topology file -- ``prepare_system.py`` parameterizes both
from scratch (Amber ff14SB via OpenMM's ``amber14-all.xml``), the same
force field OpenGBSA's own Coordinate Mode uses internally.

Run it
------

.. code-block:: bash

    cd test/configs/ppi_1gcq_test
    ./fetch_data.sh                     # ~2.8 GB, both systems
    python3 prepare_system.py --system 1gcq
    python3 prepare_system.py --system 2oob
    opengbsa 1gcq_config.yaml
    opengbsa 2oob_config.yaml

For a quick 3-frame smoke test instead of the full 300-frame run, use
``1gcq_config_test.yaml``/``2oob_config_test.yaml``.

What to notice in the config
-----------------------------

.. code-block:: yaml

    input_files:
      complex_pdb: prepared/1gcq/complex_chains.pdb
      trajectory: prepared/1gcq/complex_chains_last300.dcd

    analysis_settings:
      binding_mode: ppi
      receptor_selection: 'chainid 1'
      ligand_selection: 'chainid 0'

    forcefield_settings:
      protein_forcefield: amber14-all.xml

- No ``.prmtop``/``.tpr``/``.psf`` anywhere in ``input_files`` -- a plain
  PDB plus ``forcefield_settings.protein_forcefield`` is what triggers
  Coordinate Mode.
- ``prepare_system.py`` handles two real-world complications any
  from-scratch topology build hits: bond-inference artifacts from naive
  PDB conversion (fixed by writing a bond-free PDB and letting the
  residue-template matcher infer standard connectivity), and ambiguous
  histidine protonation states (Amber-family force fields have no generic
  ``HIS`` template -- ``HID``/``HIE``/``HIP`` must be determined from
  which hydrogens are actually present).

Full write-up (including two OpenGBSA bugs specific to the
Coordinate-Mode + ``binding_mode: ppi`` combination, found and fixed
during this validation): see ``test/configs/ppi_1gcq_test/README.md`` in
the repository.
