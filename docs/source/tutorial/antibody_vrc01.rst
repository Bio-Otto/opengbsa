Antibody-antigen: CHARMM .psf with chain selections
====================================================

This example evaluates the binding of the broadly neutralising antibody VRC01
to the HIV-1 envelope glycoprotein gp120, for the wild type and four alanine
mutants. It reads NAMD/CHARMM files in Native Mode and defines antigen and
antibody by residue ranges.

Dataset
-------

Five complexes of gp120 with the VRC01 Fv fragment (wild type and the heavy
chain mutants D99A, G54A, N100aA and W100bA), simulated with CHARMM36m by
Ovchinnikov and coworkers, Zenodo
`6879091 <https://doi.org/10.5281/zenodo.6879091>`_. Each system has a
``.psf``, a ``.pdb`` and a ``.dcd`` of 4,601 frames (about 1.05 GB). The last
300 frames are evaluated.

Run it
------

.. code-block:: bash

    cd test/configs/antibody_vrc01_test
    ./fetch_data.sh          # 5 x 1.05 GB; also provides the CHARMM36m parameter files
    opengbsa vrc01_wt_config.yaml

``fetch_data.sh`` copies the two CHARMM36m parameter files into ``params/``
from the CTLA4-peptide dataset of ``../namd_charmm_test/``, which it downloads
first if needed. ``vrc01_wt_config_test.yaml`` evaluates three frames only.

What to notice in the config
----------------------------

.. code-block:: yaml

    input_files:
      complex_pdb: raw/wt.psf
      solvated_topology: raw/wt.pdb
      trajectory: raw/wt.dcd
    analysis_settings:
      binding_mode: ppi
      receptor_selection: 'resid 0 to 350'      # gp120
      ligand_selection: 'resid 351 to 570'      # VRC01 Fv, heavy and light chain
      frame_start: 4301
      frame_end: 4601
      max_frames: 300
    forcefield_settings:
      charmm_params:
        - params/par_all36m_prot.prm
        - params/toppar_water_ions_prot.str

- The antibody has two chains, so the ligand is selected by residue range
  rather than by one chain index. Residue indices count from 0 after water and
  ions are removed.
- CHARMM's pair-specific Lennard-Jones corrections (NBFIX) are applied from the
  parameter files.

Expected result
---------------

.. code-block:: bash

    python test/configs/check_validation.py vrc01_wt vrc01_d99a vrc01_g54a vrc01_n100aa vrc01_w100ba
    python test/configs/check_validation.py vrc01_wt vrc01_d99a vrc01_g54a vrc01_n100aa vrc01_w100ba --amber-compatible

.. list-table::
   :header-rows: 1
   :widths: 22 26 26 26

   * - Case
     - OpenGBSA, default settings
     - OpenGBSA, ``amber_compatible``
     - Amber MMPBSA.py
   * - ``vrc01_wt``
     - -36.820
     - -52.459
     - -52.461
   * - ``vrc01_d99a``
     - -43.599
     - -59.162
     - -59.167
   * - ``vrc01_g54a``
     - -44.912
     - -60.695
     - -60.699
   * - ``vrc01_n100aa``
     - -48.600
     - -64.280
     - -64.280
   * - ``vrc01_w100ba``
     - -25.427
     - -40.900
     - -40.903

With ``amber_compatible`` OpenGBSA and MMPBSA.py agree within 0.01 kcal/mol.
For this large interface the Born-radius cutoff of MMPBSA.py changes the
binding enthalpy by about 16 kcal/mol, the largest difference between the two
settings among the validation cases (:doc:`../theory/mmgbsa_theory`). Both
settings rank the five variants in the same order. Provenance and the
preparation of the Amber reference: ``test/configs/antibody_vrc01_test/README.md``.
