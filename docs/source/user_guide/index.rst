User guide
==========

The user guide explains how to set up, run and read an OpenGBSA calculation.
It starts with one complete example and then describes every option: what it
does, how OpenGBSA applies it, how much it changes a result and when to
change it. :doc:`configuration` lists all options on one page.

.. grid:: 2

    .. grid-item-card::  A first example
        :link: first_example
        :link-type: doc

        The quickstart configuration line by line, and what happens during a
        run.

    .. grid-item-card::  CLI and CDI
        :link: interfaces
        :link-type: doc

        The ``opengbsa`` command and the configuration file, and how they work
        together.

    .. grid-item-card::  Input formats and modes
        :link: input_formats
        :link-type: doc

        Amber, GROMACS, CHARMM/NAMD, OpenMM and PDB input; Native and
        Coordinate Mode; choosing receptor and ligand.

    .. grid-item-card::  The energy function
        :link: energy_function
        :link-type: doc

        GB model, radii, salt, temperature, dielectric constants, non-polar
        term, Born-radius cutoff and the Amber-compatible mode.

    .. grid-item-card::  Choosing frames
        :link: frames
        :link-type: doc

        Frame range, stride, number and selection, and re-imaging.

    .. grid-item-card::  Per-residue decomposition
        :link: decomposition
        :link-type: doc

        The contribution of every receptor residue.

    .. grid-item-card::  Entropy
        :link: entropy
        :link-type: doc

        Interaction entropy, C2, quasi-harmonic and normal-mode estimates.

    .. grid-item-card::  Running an analysis
        :link: running_analysis
        :link-type: doc

        Platforms, GPUs, CPU workers, run time and logs.

    .. grid-item-card::  Output files
        :link: outputs
        :link-type: doc

        Every file a run writes and every column of the tables.

    .. grid-item-card::  Plots, report and 3D views
        :link: visualization
        :link-type: doc

        How to read the plots, the interactive report and the PyMOL view.

    .. grid-item-card::  Configuration reference
        :link: configuration
        :link-type: doc

        Every option on one page, with its default, its values and its
        equivalent in Amber MMPBSA.py.

.. toctree::
   :hidden:

   first_example
   interfaces
   input_formats
   energy_function
   frames
   decomposition
   entropy
   running_analysis
   outputs
   visualization
   configuration
