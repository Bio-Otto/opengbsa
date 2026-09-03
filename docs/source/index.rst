OpenGBSA Manual
===============

.. image:: _static/logo.png
   :align: center
   :alt: OpenGBSA Logo
   :width: 600px

**OpenGBSA** is a high-performance computational chemistry toolkit for estimating the relative binding free energies of protein-ligand and protein-protein complexes. It implements the **MM/GBSA** (Molecular Mechanics / Generalized Born Surface Area) method using **OpenMM** for rapid processing on CPUs and GPUs.

Designed for both ease of use and scientific rigor, OpenGBSA automates the complex workflows of topology preparation, trajectory processing, and energy decomposition.

.. grid:: 2

    .. grid-item-card::  🚀  Getting Started
        :link: installation
        :link-type: doc

        Installation guide and quickstart tutorial.

    .. grid-item-card::  📘  User Guide
        :link: user_guide/index
        :link-type: doc

        Detailed instructions on configuration, input structures, and output interpretation.

    .. grid-item-card::  📚  Tutorial
        :link: tutorial/index
        :link-type: doc

        Five worked examples, one per receptor/ligand combination and
        input format, each runnable end-to-end from the repository.

    .. grid-item-card::  🧠  Theory Guide
        :link: theory/index
        :link-type: doc

        Rigorous mathematical background, GB model definitions (OBC2, GBn), and entropy methods.

    .. grid-item-card::  ⚙️  API Reference
        :link: api_reference
        :link-type: doc

        Python API documentation for custom workflows and integration.

--------------------------------------------------------------------------------

Table of Contents
-----------------
.. toctree::
   :maxdepth: 2
   :caption: Introduction
   :includehidden:

   installation
   multi_engine_guide

.. toctree::
   :maxdepth: 2
   :caption: User Guide
   :includehidden:

   user_guide/index

.. toctree::
   :maxdepth: 2
   :caption: Tutorials
   :includehidden:

   tutorial/index

.. toctree::
   :maxdepth: 2
   :caption: Theory and Methods
   :includehidden:

   theory/index

.. toctree::
   :maxdepth: 2
   :caption: API and Development
   :includehidden:

   api_reference
   changelog
   contributing

Citing OpenGBSA
---------------
If you use OpenGBSA in your research, please cite:

*   **OpenGBSA**: see `CITATION.cff <https://github.com/bio-otto/opengbsa/blob/main/CITATION.cff>`_ in the repository root for citation metadata.
*   **OpenMM**: Eastman, P., et al. (2017). OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. *PLoS Computational Biology*, 13(7), e1005659.
*   **MDTraj**: McGibbon, R. T., et al. (2015). MDTraj: A modern open library for the analysis of molecular dynamics trajectories. *Biophysical Journal*, 109(8), 1528-1532.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
