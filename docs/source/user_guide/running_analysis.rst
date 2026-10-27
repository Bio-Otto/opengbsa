Running an analysis
===================

.. code-block:: bash

    opengbsa config.yaml
    opengbsa config.yaml --output-dir results     # results/mmgbsa_analysis_<date>_<time>/

A run goes through numbered steps that the log announces: configuration
validation, output setup, calculator setup, the MM/GBSA calculation, the
entropy estimate and the per-residue decomposition (when requested), the final
report and the saved results. It ends with the mean binding enthalpy and the
results folder.

Platforms
---------

All energies are computed by OpenMM on one of its platforms:

* **CUDA** on NVIDIA GPUs, in mixed precision.
* **OpenCL** on NVIDIA, AMD, Intel and Apple GPUs.
* **CPU**, with the frames divided among worker processes.
* **Reference**, a slow double-precision implementation for checks.

Without ``platform_settings.preferred_platform`` OpenGBSA takes the first that
works of CUDA, OpenCL and CPU. Before using a GPU platform it compiles and runs
a small test kernel, so a GPU that cannot run OpenMM's code (for example
because the CUDA runtime is newer than the driver) is detected at the start.
The log names the platform in the line ``Running on OpenMM platform: ...`` and
``results_summary.yaml`` records it as ``platform_used``. When a requested
platform cannot be used, the log shows a ``PLATFORM FALLBACK`` warning and the
run continues on the next platform; ``strict_platform: true`` stops it instead.
:doc:`../installation` explains how to match the CUDA runtime to the driver.

CPU runs
--------

On the CPU, frames are evaluated in parallel by ``analysis_settings.max_workers``
worker processes (default: one less than the number of cores, at most four).
Each worker loads the three OpenMM systems once. The per-residue decomposition
runs in parallel when ``platform_settings.decomposition_platform: CPU``, with
``advanced_settings.n_jobs`` processes.

Run time
--------

On a GPU a frame of a complex of 1,895 to 5,434 atoms takes 6 to 15 ms. The
per-residue decomposition is slower, because every receptor residue is
evaluated separately. The figure shows the calculation time for 100 frames of
two complexes from the OpenGBSA paper on an Intel Core i9-9900K (8 cores) and
an NVIDIA RTX 2080 SUPER, with Amber MMPBSA.py for comparison.

.. figure:: ../_static/figures/speed_benchmark.png
   :width: 100%
   :alt: Calculation time for 100 frames on CPU and GPU, OpenGBSA and MMPBSA.py

   Calculation time for 100 frames, median of three repeats. On the GPU,
   OpenGBSA was about 8 and 21 times faster than MMPBSA.py with eight MPI
   ranks on the two systems. On eight CPU cores MMPBSA.py was faster.

Logs
----

The console output of a run is the complete log. Keep it with

.. code-block:: bash

    opengbsa config.yaml 2>&1 | tee run.log

``--log-file FILE`` writes only the messages of Python's ``logging`` module.

Long runs and many systems
--------------------------

A configuration describes one complex. For a series, write one configuration
per complex and run them one after another, for example with a shell loop:

.. code-block:: bash

    for c in configs/*.yaml; do opengbsa "$c" 2>&1 | tee "${c%.yaml}.log"; done

On a cluster, submit one job per configuration. Each run writes its own
time-stamped folder, so runs never overwrite each other.
