Tutorials
=========

Six worked examples, each built from an MD dataset published by another group
and each ready to run from ``test/configs/`` in the repository. Together they
cover protein-ligand, protein-peptide, protein-protein, protein-RNA,
protein-DNA and antibody-antigen complexes, the Amber, GROMACS, CHARMM/NAMD and
PDB input formats, and both Native and Coordinate Mode. Each example has a
``fetch_data.sh`` that downloads the data from Zenodo, a configuration file, a
three-frame variant for a quick test (``*_config_test.yaml``) and a
``README.md`` with the provenance.

Every example follows the same pattern:

.. code-block:: bash

    cd test/configs/<example>
    ./fetch_data.sh                 # plus a preparation step for some examples
    opengbsa <config>.yaml

``python test/configs/check_validation.py <case>`` runs an example and
compares the result with the expected value (:doc:`../testing`).

.. list-table::
   :header-rows: 1
   :widths: 30 26 26 18

   * - Example
     - Input
     - Mode
     - Download
   * - :doc:`protein_ligand_3clpro`
     - Amber ``.prmtop``
     - Native
     - 9.1 GB
   * - :doc:`namd_charmm`
     - NAMD/CHARMM ``.psf``
     - Native
     - 321 MB
   * - :doc:`ppi_gromacs`
     - GROMACS ``.gro``
     - Coordinate
     - 2.8 GB
   * - :doc:`rna_protein`
     - Amber ``.prmtop``
     - Native
     - 44 MB
   * - :doc:`dna_protein`
     - multi-model PDB
     - Coordinate
     - 25 MB
   * - :doc:`antibody_vrc01`
     - NAMD/CHARMM ``.psf``
     - Native
     - 5.3 GB

The protein-RNA example is the quickest to try: a 44 MB download, and
``check_validation.py 3dd2`` takes about 30 seconds on a GPU. The
configurations name the platform for a direct ``opengbsa <config>.yaml`` run;
several name the CPU, on which a run of 300 frames takes minutes to an hour.
``check_validation.py`` runs every case on CUDA (:doc:`../testing`).

.. toctree::
   :hidden:

   protein_ligand_3clpro
   namd_charmm
   ppi_gromacs
   rna_protein
   dna_protein
   antibody_vrc01

Before starting, make sure OpenGBSA is installed and working
(:doc:`../installation`).
