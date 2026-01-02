Testing Guide
=============

OpenGBSA includes a comprehensive test suite to validate installation, configuration, and scientific accuracy.

Test Directory Structure
------------------------

All test files are located in the ``test/`` directory:

- ``test/unit/``: Unit tests for individual components (runners, parameterization).
- ``test/integration/``: End-to-end integration tests.
- ``test/diagnostics/``: Tools for checking system health (charges, PBC, etc.).
- ``test/utils/``: Helper scripts for file conversion and repair.
- ``test/configs/``: Sample YAML configurations for testing.
- ``test/data/``: Input data files (PDB, SDF, trajectories).

Running Tests
-------------

You can run the main analysis test runner:

.. code-block:: bash

   python test/run_analysis_test.py test/configs/test_complete_config.yaml

Running Diagnostics
-------------------

If you encounter issues, use the diagnostic scripts:

**Check Charges:**

.. code-block:: bash

   python test/diagnostics/check_charges.py --pdb complex.pdb

**Verify Surface Area (SASA):**

.. code-block:: bash

   python test/diagnostics/verify_sasa.py

**Check PBC Integrity:**

.. code-block:: bash

   python test/diagnostics/check_pbc_integrity.py

Unit Tests
----------

Run specific unit tests to verify components:

.. code-block:: bash

   python test/unit/test_runner_simple.py
