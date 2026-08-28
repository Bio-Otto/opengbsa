"""
Shared pytest configuration for test/unit/.

Warns (once, at collection time) when `test/data/` fixtures required by
some tests are missing, pointing developers to the fetch command rather
than downloading automatically -- tests should never trigger network
access on their own. Individual test modules still gate themselves with
their own `pytest.mark.skipif(not os.path.exists(...))` (see
test_charmm_loading.py for the pattern); this hook only makes the reason
for those skips discoverable without reading test output line by line.
"""
import os
from pathlib import Path

TEST_DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def pytest_collection_modifyitems(config, items):
    if TEST_DATA_DIR.exists() and any(TEST_DATA_DIR.iterdir()):
        return
    print(
        "\n"
        "NOTE: test/data/ is missing or empty -- tests that need real MD "
        "fixtures will be skipped. Run `opengbsa --fetch-test-data` to "
        "download them (see `opengbsa --list-test-data` for what's "
        "available)."
    )
