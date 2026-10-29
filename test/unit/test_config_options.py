#!/usr/bin/env python3
"""
The option table (mmgbsa.config.OPTIONS) and what validation does with it:
unknown and removed options, deprecated spellings, value checks, the
defaults a run gets, the generated template and the CLI path argument.
"""
import sys
from pathlib import Path
from unittest import mock

import pytest
import yaml

from mmgbsa import config as cfgmod
from mmgbsa.config import (ConfigManager, OPTIONS, complete_config_text, default_config,
                           effective_config, known_options)


@pytest.fixture
def files(tmp_path):
    complex_pdb = tmp_path / "complex.prmtop"
    trajectory = tmp_path / "traj.xtc"
    complex_pdb.write_text("placeholder\n")
    trajectory.write_text("placeholder\n")
    return complex_pdb, trajectory


@pytest.fixture
def validate(files):
    """validate(**sections) -> (ok, manager) for a config with these sections added."""
    complex_pdb, trajectory = files

    def run(analysis=None, **sections):
        config = {"input_files": {"complex_pdb": str(complex_pdb), "trajectory": str(trajectory)},
                  "analysis_settings": dict(analysis or {})}
        for name, value in sections.items():
            if name == "input_files":
                config["input_files"].update(value)
            else:
                config[name] = value
        manager = ConfigManager()
        manager.config = config
        return manager.validate_config(), manager
    return run


# --- unknown, removed and deprecated options ----------------------------------

def test_misspelled_option_warns_with_the_closest_name(validate):
    ok, mgr = validate({"gb_modle": "GBn2"})
    assert ok
    assert ("Unknown option 'gb_modle' in analysis_settings is ignored; did you mean 'gb_model'?"
            in mgr.get_validation_warnings())
    assert mgr.get_analysis_settings()["gb_model"] == "OBC2"
    assert "gb_modle" not in mgr.get_analysis_settings()


def test_unknown_section_warns_with_the_closest_name(validate):
    ok, mgr = validate(platform_setting={"preferred_platform": "CPU"})
    assert ok
    assert any("Unknown section 'platform_setting'" in w and "'platform_settings'" in w
               for w in mgr.get_validation_warnings())


def test_option_in_the_wrong_section_names_its_section(validate):
    ok, mgr = validate({"protein_forcefield": "charmm"})
    assert ok
    assert any("'protein_forcefield'" in w and "belongs in forcefield_settings" in w
               for w in mgr.get_validation_warnings())


@pytest.mark.parametrize("section,key", [
    ("analysis_settings", "use_cache"), ("analysis_settings", "parallel_processing"),
    ("analysis_settings", "dimer_mode"), ("analysis_settings", "receptor_subunits"),
    ("analysis_settings", "compound_name"), ("reporting_settings", "entropy_approximation"),
    ("reporting_settings", "include_plots"), ("output_settings", "analysis_name"),
    ("input_files", "receptor_resname"), ("input_files", "ligand_charge"),
])
def test_removed_options_warn_and_are_dropped(validate, section, key):
    ok, mgr = validate(**{section: {key: "x"}}) if section != "analysis_settings" \
        else validate({key: "x"})
    assert ok, mgr.get_validation_errors()
    assert any(f"Option '{key}' in {section} was removed" in w for w in mgr.get_validation_warnings())
    assert key not in mgr.get_config()[section]


@pytest.mark.parametrize("section", ["reproducibility_settings", "visualization_settings"])
def test_removed_sections_warn(validate, section):
    ok, mgr = validate(**{section: {"anything": True}})
    assert ok
    assert any(f"Section '{section}' was removed" in w for w in mgr.get_validation_warnings())
    assert section not in mgr.get_config()


def test_params_output_directory_is_a_deprecated_alias(validate):
    ok, mgr = validate(params={"output_directory": "from_params"},
                       output_settings={"output_directory": "from_output_settings"})
    assert ok
    assert mgr.get_config()["output_settings"]["output_directory"] == "from_params"
    assert "params" not in mgr.get_config()
    assert any("params.output_directory is deprecated; use output_settings.output_directory" in w
               for w in mgr.get_validation_warnings())


@pytest.mark.parametrize("section,key", [("analysis_settings", "per_residue_decomposition"),
                                         ("advanced_settings", "run_per_residue_decomposition")])
def test_per_residue_decomposition_aliases(validate, section, key):
    ok, mgr = (validate({key: True}) if section == "analysis_settings"
               else validate(advanced_settings={key: True}))
    assert ok, mgr.get_validation_errors()
    assert mgr.get_analysis_settings()["run_per_residue_decomposition"] is True
    assert any(f"{section}.{key} is deprecated; use analysis_settings.run_per_residue_decomposition"
               in w for w in mgr.get_validation_warnings())


def test_legacy_input_section(tmp_path, files):
    complex_pdb, trajectory = files
    path = tmp_path / "legacy.yaml"
    path.write_text(yaml.safe_dump({"input": {"topology": str(complex_pdb), "trajectory": str(trajectory)},
                                    "analysis_settings": {}}))
    mgr = ConfigManager(str(path))
    assert mgr.validate_config(), mgr.get_validation_errors()
    assert mgr.get_input_files()["complex_pdb"] == str(complex_pdb)
    assert any("top-level 'input' section is deprecated" in w for w in mgr.get_validation_warnings())


# --- value checks ---------------------------------------------------------------

@pytest.mark.parametrize("key,value", [
    ("max_workers", 0), ("max_workers", 2.5), ("print_interval", 0), ("print_interval", None),
    ("solute_dielectric", 0), ("solute_dielectric", None), ("solute_dielectric", "1"),
    ("solvent_dielectric", 0), ("solvent_dielectric", -80),
    ("temperature", 0), ("temperature", 1001), ("temperature", True),
    ("nmode_drms", None), ("nmode_maxcyc", None), ("nmode_maxcyc", 0), ("entropy_method", None),
    ("run_entropy_analysis", "yes"), ("run_per_residue_decomposition", "true"),
    ("reimage_trajectory", "false"), ("amber_compatible", 1), ("report_raw_energies", 0),
    ("energy_decomposition", "no"), ("qha_analyze_complex", None),
    ("frame_start", -1), ("frame_start", 1.5), ("frame_end", 0), ("frame_stride", 0),
    ("max_frames", 50.0), ("decomp_frames", 0), ("plot_top_residues", 0),
    ("charge_method", "resp"), ("binding_mode", "dimer"), ("verbose", -1), ("verbose", "yes"),
    ("random_seed", "42"),
])
def test_invalid_analysis_value_is_an_error_that_says_what_is_allowed(validate, key, value):
    ok, mgr = validate({key: value})
    assert not ok
    errors = mgr.get_validation_errors()
    assert any(f"analysis_settings.{key}" in e and ("must be" in e) for e in errors), errors


@pytest.mark.parametrize("section,key,value", [
    ("advanced_settings", "n_jobs", 0), ("advanced_settings", "n_jobs", -2),
    ("platform_settings", "strict_platform", "true"), ("platform_settings", "preferred_platform", "metal"),
    ("platform_settings", "decomposition_platform", 1), ("reporting_settings", "export_format", "pdf"),
    ("reporting_settings", "include_all_plots", "yes"), ("forcefield_settings", "protein_forcefield", "gaff"),
])
def test_invalid_value_in_other_sections(validate, section, key, value):
    ok, mgr = validate(**{section: {key: value}})
    assert not ok
    assert any(f"{section}.{key}" in e and "must be" in e for e in mgr.get_validation_errors())


@pytest.mark.parametrize("section,key,value,expected", [
    ("analysis_settings", "gb_radii_set", None, "auto"),
    ("analysis_settings", "verbose", 0, False),
    ("analysis_settings", "verbose", 1, True),
    ("analysis_settings", "verbose", 2, True),
    ("forcefield_settings", "protein_forcefield", "amber14-all.xml", "amber14"),
    ("forcefield_settings", "protein_forcefield", "CHARMM", "charmm"),
    ("analysis_settings", "temperature", 310, 310.0),
    ("analysis_settings", "solvent_dielectric", None, None),
    ("analysis_settings", "charge_method", "MMFF94", "mmff94"),
    ("analysis_settings", "nmode_precision", None, "double"),
    ("analysis_settings", "random_seed", None, 42),
    ("platform_settings", "preferred_platform", "cuda", "CUDA"),
    ("platform_settings", "preferred_platform", "opencl", "OpenCL"),
    ("platform_settings", "decomposition_platform", "cpu", "CPU"),
    ("platform_settings", "preferred_platform", "AUTO", None),
    ("advanced_settings", "n_jobs", -1, -1),
])
def test_values_are_normalised(validate, section, key, value, expected):
    ok, mgr = (validate({key: value}) if section == "analysis_settings"
               else validate(**{section: {key: value}}))
    assert ok, mgr.get_validation_errors()
    assert mgr.get_config()[section][key] == expected


def test_random_selection_without_seed_uses_42(validate):
    ok, mgr = validate({"frame_selection": "random"})
    assert ok
    assert mgr.get_analysis_settings()["random_seed"] == 42


def test_entropy_method_null_says_how_to_switch_entropy_off(validate):
    ok, mgr = validate({"entropy_method": None})
    assert not ok
    assert any("use 'none' or run_entropy_analysis: false" in e for e in mgr.get_validation_errors())


@pytest.mark.parametrize("value", [8.0, 12, 0])
def test_a_nonbonded_cutoff_is_an_error_pointing_to_rgbmax(validate, value):
    ok, mgr = validate({"nonbonded_cutoff": value})
    assert not ok
    assert any("nonbonded_cutoff is no longer supported" in e and "rgbmax" in e
               for e in mgr.get_validation_errors())


@pytest.mark.parametrize("value", [None, 999, 1000.0])
def test_nonbonded_cutoff_meaning_no_cutoff_only_warns(validate, value):
    ok, mgr = validate({"nonbonded_cutoff": value})
    assert ok
    assert any("nonbonded_cutoff is no longer supported" in w for w in mgr.get_validation_warnings())
    assert "nonbonded_cutoff" not in mgr.get_analysis_settings()


def test_decomp_frames_default_does_not_conflict_with_small_max_frames(validate):
    ok, mgr = validate({"max_frames": 5, "run_per_residue_decomposition": True})
    assert ok, mgr.get_validation_errors()


def test_missing_optional_input_file_is_an_error(validate, tmp_path):
    ok, mgr = validate(input_files={"ligand_mol": str(tmp_path / "missing.sdf")})
    assert not ok
    assert any("input_files.ligand_mol" in e for e in mgr.get_validation_errors())


# --- defaults, table and template ---------------------------------------------

def test_validated_config_has_every_option_with_its_default(validate):
    ok, mgr = validate({"gb_model": "GBn2"})
    assert ok
    config = mgr.get_config()
    defaults = default_config()
    assert set(config) == set(OPTIONS)
    for section, options in OPTIONS.items():
        assert set(config[section]) == set(options)
    assert config["analysis_settings"]["gb_model"] == "GBn2"
    assert config["analysis_settings"]["max_frames"] == defaults["analysis_settings"]["max_frames"] == 50
    assert config["reporting_settings"] == defaults["reporting_settings"]


def test_effective_config_of_the_python_interface_fills_defaults_and_maps_aliases():
    raw = {"input_files": {"complex_pdb": "c.prmtop", "trajectory": "t.xtc"},
           "analysis_settings": {"per_residue_decomposition": True, "use_cache": False},
           "params": {"output_directory": "out"}, "platform_settings": {"preferred_platform": "cuda"}}
    eff = effective_config(raw)
    assert eff["analysis_settings"]["run_per_residue_decomposition"] is True
    assert "use_cache" not in eff["analysis_settings"]
    assert eff["output_settings"]["output_directory"] == "out"
    assert eff["platform_settings"]["preferred_platform"] == "CUDA"
    assert raw["analysis_settings"] == {"per_residue_decomposition": True, "use_cache": False}


def test_option_table_is_importable_and_consistent():
    table = known_options()
    assert "gb_model" in table["analysis_settings"]
    assert "n_jobs" in table["advanced_settings"]
    assert "nonbonded_cutoff" not in table["analysis_settings"]
    for removed_section, options in cfgmod.REMOVED_OPTIONS.items():
        assert not set(options) & set(table[removed_section])


def test_complete_template_lists_exactly_the_table_with_its_defaults(validate, files):
    template = yaml.safe_load(complete_config_text())
    defaults = default_config()
    assert {s: set(o) for s, o in template.items()} == {s: set(o) for s, o in defaults.items()}
    for section, options in OPTIONS.items():
        for name, opt in options.items():
            if opt.example is None:
                assert template[section][name] == defaults[section][name], (section, name)
    # the template itself is a valid configuration once the two required paths exist
    complex_pdb, trajectory = files
    template["input_files"].update(complex_pdb=str(complex_pdb), trajectory=str(trajectory))
    mgr = ConfigManager()
    mgr.config = template
    assert mgr.validate_config(), mgr.get_validation_errors()
    assert mgr.get_validation_warnings() == []


def test_create_complete_config_writes_to_the_given_path(tmp_path, monkeypatch):
    from mmgbsa import cli
    monkeypatch.chdir(tmp_path)
    target = tmp_path / "sub" / "my_template.yaml"
    monkeypatch.setattr(sys, "argv", ["opengbsa", "--create-complete-config", str(target)])
    assert cli.main() == 0
    assert target.read_text() == complete_config_text()
    assert not (tmp_path / "complete_mmgbsa_config.yaml").exists()

    monkeypatch.setattr(sys, "argv", ["opengbsa", "--create-complete-config"])
    assert cli.main() == 0
    assert (tmp_path / "complete_mmgbsa_config.yaml").exists()

    short = tmp_path / "short.yaml"
    monkeypatch.setattr(sys, "argv", ["opengbsa", "--create-config", str(short)])
    assert cli.main() == 0
    assert yaml.safe_load(short.read_text())["analysis_settings"]["gb_model"] == "OBC2"


# --- runner side -------------------------------------------------------------

def test_run_folders_started_in_the_same_second_are_distinct(tmp_path):
    from mmgbsa.runner import make_run_directory
    from datetime import datetime as real_datetime

    class Frozen(real_datetime):
        @classmethod
        def now(cls, tz=None):
            return real_datetime(2026, 10, 29, 12, 0, 0)

    with mock.patch("mmgbsa.runner.datetime", Frozen):
        first = make_run_directory(tmp_path, "analysis")
        second = make_run_directory(tmp_path, "analysis")
        third = make_run_directory(tmp_path, "analysis")
    assert first.name == "analysis_20261029_120000"
    assert second.name == "analysis_20261029_120000_2"
    assert third.name == "analysis_20261029_120000_3"
    assert first.is_dir() and second.is_dir() and third.is_dir()


def test_runner_checks_only_known_input_files(tmp_path, files):
    from mmgbsa.runner import MMGBSARunner
    complex_pdb, trajectory = files
    runner = MMGBSARunner({"input_files": {"complex_pdb": str(complex_pdb), "trajectory": str(trajectory),
                                           "receptor_resname": "PRO"},
                           "analysis_settings": {}}, output_dir=str(tmp_path / "run"))
    assert runner._validate_input_files()
    runner.config["input_files"]["ligand_topology"] = str(tmp_path / "missing.prmtop")
    assert not runner._validate_input_files()


def test_runner_rejects_an_invalid_configuration_file(tmp_path, files):
    from mmgbsa.runner import MMGBSARunner
    complex_pdb, trajectory = files
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump({"input_files": {"complex_pdb": str(complex_pdb),
                                                    "trajectory": str(trajectory)},
                                    "analysis_settings": {"print_interval": 0}}))
    with pytest.raises(ValueError, match="print_interval"):
        MMGBSARunner(str(path))


# --- charge_method and platform names in the calculator --------------------------

def test_calculator_rejects_an_unknown_charge_method():
    from mmgbsa.mmgbsa_core import GBSACalculator
    with pytest.raises(ValueError, match="charge_method"):
        GBSACalculator(charge_method="resp")
    assert GBSACalculator(charge_method="MMFF94").charge_method == "mmff94"


class _Molecule:
    def __init__(self, charges=None, fail=False):
        self.partial_charges = charges
        self.fail = fail
        self.calls = []

    def assign_partial_charges(self, partial_charge_method):
        self.calls.append(partial_charge_method)
        if self.fail:
            raise RuntimeError("sqm not found")
        self.partial_charges = [0.0]


@pytest.mark.parametrize("method", ["am1bcc", "mmff94", "gasteiger"])
def test_ligand_without_charges_gets_charge_method(method):
    from mmgbsa.mmgbsa_core import GBSACalculator
    calc = GBSACalculator(charge_method=method)
    mol = _Molecule()
    assert calc.assign_ligand_charges(mol, "lig.sdf") == method
    assert mol.calls == [method]
    assert any(method in a for a in calc.physics_assumptions)


def test_ligand_file_charges_are_kept():
    from mmgbsa.mmgbsa_core import GBSACalculator
    calc = GBSACalculator()
    mol = _Molecule(charges=[0.1, -0.1])
    assert calc.assign_ligand_charges(mol, "lig.sdf") == "file"
    assert mol.calls == []


def test_failing_am1bcc_stops_without_fallback():
    from mmgbsa.mmgbsa_core import GBSACalculator
    calc = GBSACalculator()
    mol = _Molecule(fail=True)
    with pytest.raises(RuntimeError, match="mmff94"):
        calc.assign_ligand_charges(mol, "lig.sdf")
    assert mol.calls == ["am1bcc"]


def test_calculator_platform_names_in_any_case():
    from mmgbsa.mmgbsa_core import GBSACalculator
    calc = GBSACalculator()
    calc.set_platform_settings({"preferred_platform": "cpu", "decomposition_platform": "Cuda"})
    assert calc.preferred_platform == "CPU" and calc.platform_preference == "CPU"
    assert calc.decomposition_platform == "CUDA"
