Contributing
============

Bug reports, questions and pull requests are welcome at
https://github.com/Bio-Otto/opengbsa.

Development setup
-----------------

Install from source in editable mode as described in :doc:`installation`
(section 6), then check that the tests pass:

.. code-block:: bash

    git clone https://github.com/Bio-Otto/opengbsa.git
    cd opengbsa
    # create and activate the development environment (installation, section 6)
    pip install --no-deps -e .
    python -m pytest test/unit -q

The tests that need MD trajectories are skipped until the test data are
downloaded with ``opengbsa --fetch-test-data`` (:doc:`testing`). Continuous
integration runs the unit tests on Linux and macOS, builds this documentation
with warnings treated as errors and builds the conda packages.

Building the documentation
--------------------------

.. code-block:: bash

    pip install -r requirements.txt
    python docs/scripts/make_figures.py      # only after changing a diagram or plot
    sphinx-build -b html -W docs/source docs/_build/html

Pull requests
-------------

* Add a regression test under ``test/unit/`` for every fix and feature. Tests
  must not download data; gate them on the presence of their fixture with
  ``pytest.mark.skipif``.
* Describe a change in ``CHANGELOG.md`` under "Unreleased".
* A change that alters computed energies needs evidence: a comparison with
  Amber MMPBSA.py or with the previous release on the validation cases
  (``test/configs/check_validation.py``).

Releases
--------

A release is a tag ``vX.Y.Z`` on ``main``. The tag builds and tests the conda
packages for Linux and macOS, uploads them to the ``bio-otto`` channel, creates
the GitHub release and archives it on Zenodo.
