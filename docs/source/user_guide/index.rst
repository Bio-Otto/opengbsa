User guide
==========

The user guide explains how to set up, run and read an OpenGBSA calculation.

.. grid:: 2

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

    .. grid-item-card::  Configuration reference
        :link: configuration
        :link-type: doc

        Every option, with its default and its meaning.

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

.. toctree::
   :hidden:

   interfaces
   input_formats
   configuration
   running_analysis
   outputs
   visualization
