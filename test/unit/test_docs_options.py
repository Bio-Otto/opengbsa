"""The documentation and the option table of mmgbsa/config.py agree.

Every option OpenGBSA reads is listed in the configuration reference, the
reference names no option that does not exist, and every YAML example in the
documentation uses only existing options in their right sections.
"""
import re
import textwrap
from pathlib import Path

import pytest
import yaml

from mmgbsa import config

DOCS = Path(__file__).resolve().parents[2] / "docs" / "source"
REFERENCE = DOCS / "user_guide" / "configuration.rst"
KNOWN = config.known_options()
ALL_OPTIONS = {name for options in KNOWN.values() for name in options}


def yaml_blocks():
    """(file, text) of every ``.. code-block:: yaml`` in the documentation."""
    for rst in sorted(DOCS.rglob("*.rst")):
        lines = rst.read_text().splitlines()
        for i, line in enumerate(lines):
            if line.strip() != ".. code-block:: yaml":
                continue
            block = []
            for body in lines[i + 1:]:
                if body.strip() and not body.startswith((" ", "\t")):
                    break
                block.append(body)
            yield rst.relative_to(DOCS), textwrap.dedent("\n".join(block)).strip("\n")


def test_every_option_is_in_the_reference():
    text = REFERENCE.read_text()
    missing = [f"{s}.{o}" for s, opts in KNOWN.items() for o in opts if not re.search(rf"\b{o}\b", text)]
    assert not missing, f"options missing from {REFERENCE.name}: {missing}"


def test_reference_tables_name_only_known_options():
    text = REFERENCE.read_text()
    tables = text.split("Equivalents in Amber MMPBSA.py")[0]
    named = set(re.findall(r":ref:`([a-z_0-9, ]+) <opt-", tables))
    named = {n.strip() for group in named for n in group.split(",")}
    named |= {m for m in re.findall(r"\* - ``([a-z_0-9]+)``", tables)}
    assert named - ALL_OPTIONS == set()


@pytest.mark.parametrize("where,block", list(yaml_blocks()), ids=lambda v: str(v)[:40])
def test_yaml_examples_use_known_options(where, block):
    data = yaml.safe_load(block)
    assert isinstance(data, dict), f"{where}: not a mapping"
    for key, value in data.items():
        if key in KNOWN:                                   # a section with its options
            assert isinstance(value, dict), f"{where}: section {key} is not a mapping"
            unknown = set(value) - set(KNOWN[key])
            assert not unknown, f"{where}: unknown options in {key}: {unknown}"
        else:                                              # a fragment: options alone
            assert key in ALL_OPTIONS, f"{where}: unknown option or section {key!r}"


def test_quickstart_configuration_is_clean():
    manager = config.ConfigManager(str(Path(__file__).resolve().parents[2] / "examples/quickstart/quickstart.yaml"))
    manager.validate_config()
    assert not [w for w in manager.get_validation_warnings() if "unknown" in w.lower() or "removed" in w.lower()]
