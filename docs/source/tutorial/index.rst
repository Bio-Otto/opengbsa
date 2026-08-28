Tutorial
========

Five worked examples, each built from an independent, published MD
dataset, each self-contained under ``test/configs/`` in the OpenGBSA
repository. Every example ships its own data-fetch script, its own YAML
config, and (where the source dataset provides no ready-to-use topology)
a system-preparation script -- so each one can be reproduced from a clean
checkout with three commands:

.. code-block:: bash

    cd test/configs/<example>_test
    ./fetch_data.sh
    opengbsa <config>.yaml

Together they cover every receptor/ligand *kind* combination OpenGBSA
supports and every native input format it accepts (Amber, GROMACS, and
CHARMM/NAMD), plus both of its two system-building paths, Native Mode
(load an existing topology) and Coordinate Mode (parameterize from plain
coordinates via OpenMM's ``ForceField``).

.. grid:: 2

    .. grid-item-card::  🧬  Protein-ligand batch (Amber)
        :link: protein_ligand_3clpro
        :link-type: doc

        Seven independent small-molecule inhibitors against one protease,
        from real Amber topologies. Native Mode, the standard case.

    .. grid-item-card::  🧵  Protein-peptide (NAMD/CHARMM)
        :link: namd_charmm
        :link-type: doc

        A peptide-protein complex loaded directly from native NAMD
        ``.psf``/CHARMM36 parameter files -- no topology reconstruction.

    .. grid-item-card::  🤝  Protein-protein (GROMACS)
        :link: ppi_gromacs
        :link-type: doc

        Two complexes with no topology provided at all -- built from
        scratch via OpenMM's ``ForceField`` (Coordinate Mode) and split by
        ``binding_mode: ppi`` chain selection instead of a ligand residue.

    .. grid-item-card::  🧶  Protein-RNA (Amber-native)
        :link: rna_protein
        :link-type: doc

        A protein bound to a 25-nucleotide RNA chain, using explicit
        ``receptor_topology``/``ligand_topology`` split from one complex
        prmtop.

    .. grid-item-card::  🧬  Protein-DNA (GROMACS-derived)
        :link: dna_protein
        :link-type: doc

        A protein bound to a 42-nucleotide DNA duplex, delivered as a
        single multi-model PDB trajectory -- Coordinate Mode again, this
        time on a nucleic acid.

.. toctree::
   :hidden:

   protein_ligand_3clpro
   namd_charmm
   ppi_gromacs
   rna_protein
   dna_protein

Before starting any of these, make sure OpenGBSA itself is installed and
working -- see :doc:`../installation`. Each tutorial page links back to
the corresponding ``test/configs/`` folder's own ``README.md``, which is
the authoritative, most detailed version of these instructions and is
kept in sync with the actual config files.
