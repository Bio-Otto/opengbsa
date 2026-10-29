#!/usr/bin/env python3
"""
Configuration validation as a user meets it on the command line.

Up to 0.0.19 `temperature: 300` (an integer) failed validation, an omitted
gb_model was an error although every code path defaults to OBC2, a section
left empty in YAML crashed the run, and an unknown frame_selection,
entropy_method or decomposition_method surfaced only late in the run.
"""
import pytest

from mmgbsa.config import ConfigManager


@pytest.fixture
def make_config(tmp_path):
    for name in ("complex.prmtop", "traj.xtc"):
        (tmp_path / name).write_text("placeholder\n")

    def make(analysis_settings, extra=""):
        path = tmp_path / "config.yaml"
        path.write_text(
            "input_files:\n"
            f"  complex_pdb: {tmp_path / 'complex.prmtop'}\n"
            f"  trajectory: {tmp_path / 'traj.xtc'}\n"
            f"analysis_settings:\n{analysis_settings}{extra}")
        manager = ConfigManager(str(path))
        return manager, manager.validate_config()
    return make


def test_integer_temperature_is_accepted_as_float(make_config):
    manager, ok = make_config("  temperature: 300\n  salt_concentration: 0\n")
    assert ok, manager.get_validation_errors()
    settings = manager.get_config()["analysis_settings"]
    assert settings["temperature"] == 300.0 and isinstance(settings["temperature"], float)
    assert settings["salt_concentration"] == 0.0


def test_gb_model_defaults_to_obc2(make_config):
    manager, ok = make_config("  max_frames: 10\n")
    assert ok, manager.get_validation_errors()
    assert manager.get_config()["analysis_settings"]["gb_model"] == "OBC2"


def test_empty_section_is_an_empty_section(make_config):
    manager, ok = make_config("  max_frames: 10\n", "platform_settings:\n")
    assert ok, manager.get_validation_errors()
    # an empty section takes the default of each of its options
    assert manager.get_config()["platform_settings"] == {
        "preferred_platform": None, "decomposition_platform": None, "strict_platform": False}


@pytest.mark.parametrize("key,value", [("frame_selection", "every_other"),
                                       ("entropy_method", "gaussian"),
                                       ("decomposition_method", "pairwise")])
def test_unknown_option_value_is_an_error(make_config, key, value):
    manager, ok = make_config(f"  {key}: {value}\n")
    assert not ok
    assert any(key in error for error in manager.get_validation_errors())
