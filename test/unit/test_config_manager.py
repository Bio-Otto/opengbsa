#!/usr/bin/env python3
"""
Unit tests for mmgbsa.config.ConfigManager -- the real config-loading and
validation class used by the CLI (via mmgbsa.runner), not a reimplementation.

Complements test_yaml_config.py/test_complete_config.py, which exercise
their own local dict-building/validation logic without ever importing
ConfigManager itself.
"""
import os
import yaml
import pytest

from mmgbsa.config import ConfigManager


def write_yaml(path, data):
    with open(path, "w") as f:
        yaml.safe_dump(data, f)
    return str(path)


def make_valid_config(tmp_path, complex_pdb, trajectory):
    return {
        "input_files": {
            "complex_pdb": str(complex_pdb),
            "trajectory": str(trajectory),
        },
        "analysis_settings": {
            "temperature": 300.0,
            "gb_model": "OBC2",
            "salt_concentration": 0.15,
            "max_frames": 50,
        },
    }


@pytest.fixture
def existing_files(tmp_path):
    complex_pdb = tmp_path / "complex.pdb"
    trajectory = tmp_path / "traj.xtc"
    complex_pdb.write_text("dummy pdb")
    trajectory.write_bytes(b"dummy traj")
    return complex_pdb, trajectory


# --- load_config / _normalize_config ---------------------------------------

def test_load_config_missing_file_returns_false():
    mgr = ConfigManager()
    assert mgr.load_config("/no/such/file.yaml") is False


def test_load_config_invalid_yaml_returns_false(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("input_files: [unbalanced\n  brackets")
    mgr = ConfigManager()
    assert mgr.load_config(str(bad)) is False


def test_load_config_valid_file_returns_true(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    config_path = write_yaml(tmp_path / "config.yaml", make_valid_config(tmp_path, complex_pdb, trajectory))
    mgr = ConfigManager()
    assert mgr.load_config(config_path) is True
    assert mgr.config["input_files"]["complex_pdb"] == str(complex_pdb)


def test_normalize_legacy_input_section_maps_to_input_files(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    legacy = {
        "input": {
            "topology": str(complex_pdb),
            "trajectory": str(trajectory),
        },
        "analysis_settings": {"temperature": 300.0, "gb_model": "OBC2",
                               "salt_concentration": 0.15, "max_frames": 50},
    }
    config_path = write_yaml(tmp_path / "legacy.yaml", legacy)
    mgr = ConfigManager()
    mgr.load_config(config_path)
    # 'topology' under the legacy 'input' section maps to 'complex_pdb'
    assert mgr.config["input_files"]["complex_pdb"] == str(complex_pdb)
    assert mgr.config["input_files"]["trajectory"] == str(trajectory)


def test_normalize_does_not_override_existing_input_files(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    cfg = make_valid_config(tmp_path, complex_pdb, trajectory)
    cfg["input"] = {"topology": "/should/be/ignored.pdb"}
    config_path = write_yaml(tmp_path / "both.yaml", cfg)
    mgr = ConfigManager()
    mgr.load_config(config_path)
    assert mgr.config["input_files"]["complex_pdb"] == str(complex_pdb)


# --- validate_config: required sections -------------------------------------

def test_validate_missing_required_sections():
    mgr = ConfigManager()
    mgr.config = {}
    assert mgr.validate_config() is False
    errors = mgr.get_validation_errors()
    assert any("input_files" in e for e in errors)
    assert any("analysis_settings" in e for e in errors)


def test_validate_valid_config_passes(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    mgr.config = make_valid_config(tmp_path, complex_pdb, trajectory)
    assert mgr.validate_config() is True
    assert mgr.get_validation_errors() == []


# --- _validate_input_files ---------------------------------------------------

def test_validate_input_files_missing_required_key(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    cfg = make_valid_config(tmp_path, complex_pdb, trajectory)
    del cfg["input_files"]["trajectory"]
    mgr.config = cfg
    assert mgr.validate_config() is False
    assert any("trajectory" in e for e in mgr.get_validation_errors())


def test_validate_input_files_nonexistent_path_is_invalid(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    cfg = make_valid_config(tmp_path, complex_pdb, trajectory)
    cfg["input_files"]["complex_pdb"] = str(tmp_path / "does_not_exist.pdb")
    mgr.config = cfg
    assert mgr.validate_config() is False
    assert any("Invalid file path" in e for e in mgr.get_validation_errors())


def test_validate_large_file_produces_warning_not_error(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    # Sparse-allocate a >1GB file without actually writing 1GB of data.
    big = tmp_path / "big.pdb"
    with open(big, "wb") as f:
        f.seek(1024 * 1024 * 1024 + 1)
        f.write(b"\0")
    mgr = ConfigManager()
    cfg = make_valid_config(tmp_path, big, trajectory)
    mgr.config = cfg
    assert mgr.validate_config() is True
    assert any("Large file" in w for w in mgr.get_validation_warnings())


# --- _validate_analysis_settings ---------------------------------------------

def test_validate_analysis_settings_injects_defaults_for_missing_optional_params(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    cfg = make_valid_config(tmp_path, complex_pdb, trajectory)
    del cfg["analysis_settings"]["temperature"]
    del cfg["analysis_settings"]["salt_concentration"]
    del cfg["analysis_settings"]["max_frames"]
    mgr.config = cfg
    assert mgr.validate_config() is True
    settings = mgr.get_analysis_settings()
    assert settings["temperature"] == 300.0
    assert settings["salt_concentration"] == 0.15
    assert settings["max_frames"] == 50


def test_validate_analysis_settings_missing_gb_model_is_an_error(tmp_path, existing_files):
    # gb_model has no injected default (unlike temperature/salt/max_frames).
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    cfg = make_valid_config(tmp_path, complex_pdb, trajectory)
    del cfg["analysis_settings"]["gb_model"]
    mgr.config = cfg
    assert mgr.validate_config() is False
    assert any("gb_model" in e for e in mgr.get_validation_errors())


@pytest.mark.parametrize("gb_model", ["OBC1", "OBC2", "HCT", "GBn", "GBn2"])
def test_validate_all_documented_gb_models_are_accepted(tmp_path, existing_files, gb_model):
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    cfg = make_valid_config(tmp_path, complex_pdb, trajectory)
    cfg["analysis_settings"]["gb_model"] = gb_model
    mgr.config = cfg
    assert mgr.validate_config() is True


def test_validate_unknown_gb_model_is_rejected(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    cfg = make_valid_config(tmp_path, complex_pdb, trajectory)
    cfg["analysis_settings"]["gb_model"] = "NOT_A_REAL_MODEL"
    mgr.config = cfg
    assert mgr.validate_config() is False


def test_validate_temperature_out_of_range_is_rejected(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    cfg = make_valid_config(tmp_path, complex_pdb, trajectory)
    cfg["analysis_settings"]["temperature"] = -10.0
    mgr.config = cfg
    assert mgr.validate_config() is False


def test_validate_wrong_type_is_rejected(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    cfg = make_valid_config(tmp_path, complex_pdb, trajectory)
    cfg["analysis_settings"]["temperature"] = "hot"
    mgr.config = cfg
    assert mgr.validate_config() is False


@pytest.mark.parametrize("binding_mode", ["standard", "dimer_ligand", "ppi"])
def test_validate_documented_binding_modes_are_accepted(tmp_path, existing_files, binding_mode):
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    cfg = make_valid_config(tmp_path, complex_pdb, trajectory)
    cfg["analysis_settings"]["binding_mode"] = binding_mode
    mgr.config = cfg
    assert mgr.validate_config() is True
    assert mgr.get_analysis_settings()["binding_mode"] == binding_mode


def test_validate_unknown_binding_mode_is_rejected(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    cfg = make_valid_config(tmp_path, complex_pdb, trajectory)
    cfg["analysis_settings"]["binding_mode"] = "not_a_real_mode"
    mgr.config = cfg
    assert mgr.validate_config() is False


def test_validate_default_binding_mode_is_standard_when_unset(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    mgr.config = make_valid_config(tmp_path, complex_pdb, trajectory)
    assert mgr.validate_config() is True
    assert mgr.get_analysis_settings()["binding_mode"] == "standard"


# --- _validate_cross_fields ---------------------------------------------------

def test_validate_frame_end_before_frame_start_is_rejected(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    cfg = make_valid_config(tmp_path, complex_pdb, trajectory)
    cfg["analysis_settings"]["frame_start"] = 100
    cfg["analysis_settings"]["frame_end"] = 50
    mgr.config = cfg
    assert mgr.validate_config() is False
    assert any("frame_end" in e for e in mgr.get_validation_errors())


def test_validate_decomp_frames_greater_than_max_frames_is_rejected(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    cfg = make_valid_config(tmp_path, complex_pdb, trajectory)
    cfg["analysis_settings"]["max_frames"] = 10
    cfg["analysis_settings"]["decomp_frames"] = 20
    mgr.config = cfg
    assert mgr.validate_config() is False
    assert any("decomp_frames" in e for e in mgr.get_validation_errors())


def test_validate_random_frame_selection_requires_seed(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    cfg = make_valid_config(tmp_path, complex_pdb, trajectory)
    cfg["analysis_settings"]["frame_selection"] = "random"
    mgr.config = cfg
    assert mgr.validate_config() is False
    assert any("random_seed" in e for e in mgr.get_validation_errors())


def test_validate_random_frame_selection_with_seed_passes(tmp_path, existing_files):
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    cfg = make_valid_config(tmp_path, complex_pdb, trajectory)
    cfg["analysis_settings"]["frame_selection"] = "random"
    cfg["analysis_settings"]["random_seed"] = 42
    mgr.config = cfg
    assert mgr.validate_config() is True


# --- get_* accessors return shallow copies ------------------------------------

def test_get_config_returns_a_shallow_copy(tmp_path, existing_files):
    """get_config() is `self.config.copy()` -- a shallow copy, so replacing
    a top-level key on the returned dict doesn't affect the original, but
    mutating a nested dict-in-place (e.g. ['input_files']['x'] = ...) does,
    since both dicts still share the same nested `input_files` object."""
    complex_pdb, trajectory = existing_files
    mgr = ConfigManager()
    mgr.config = make_valid_config(tmp_path, complex_pdb, trajectory)
    snapshot = mgr.get_config()

    # Replacing a top-level key on the copy leaves the original untouched.
    snapshot["input_files"] = {"replaced": True}
    assert "replaced" not in mgr.config["input_files"]

    # But nested dicts are shared references (shallow copy), so mutating
    # them in-place DOES affect the original -- documenting this behavior
    # rather than asserting a deep-copy guarantee the code doesn't provide.
    snapshot2 = mgr.get_config()
    snapshot2["input_files"]["complex_pdb"] = "mutated"
    assert mgr.config["input_files"]["complex_pdb"] == "mutated"


# --- create_default_config / create_complete_config ---------------------------

def test_create_default_config_writes_valid_loadable_yaml(tmp_path):
    out_path = tmp_path / "default_config.yaml"
    mgr = ConfigManager()
    assert mgr.create_default_config(str(out_path)) is True
    assert out_path.exists()
    with open(out_path) as f:
        loaded = yaml.safe_load(f)
    assert "input_files" in loaded
    assert "analysis_settings" in loaded
    assert loaded["analysis_settings"]["gb_model"] == "OBC2"


def test_create_complete_config_writes_valid_loadable_yaml(tmp_path):
    out_path = tmp_path / "complete_config.yaml"
    mgr = ConfigManager()
    assert mgr.create_complete_config(str(out_path)) is True
    assert out_path.exists()
    with open(out_path) as f:
        loaded = yaml.safe_load(f)
    assert "input_files" in loaded
    assert "analysis_settings" in loaded
    assert "forcefield_settings" in loaded
    assert "platform_settings" in loaded


def test_create_default_config_fails_gracefully_on_unwritable_path():
    mgr = ConfigManager()
    assert mgr.create_default_config("/nonexistent_dir_xyz/config.yaml") is False
