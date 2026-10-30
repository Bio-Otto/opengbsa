#!/usr/bin/env python3
"""
A failed per-residue decomposition is reported as a failure (before, the
MM/GBSA results were stored under 'decomposition' and the run said
'Per-residue decomposition completed'); with strict_platform the run stops
with an error. Also: analysis_config.yaml holds the effective configuration.
"""
import pandas as pd
import pytest
import yaml

import mmgbsa.runner as runner_module
from mmgbsa.runner import MMGBSARunner


def _fake_mmgbsa(**kwargs):
    df = pd.DataFrame({"frame": [0, 1], "binding_energy": [-10.0, -12.0],
                       "delta_vdw": [-20.0, -21.0], "delta_elec": [-5.0, -6.0]})
    return {"mean_binding_energy": -11.0, "std_error": 1.0, "std_dev": 1.4, "n_frames": 2,
            "gb_model": "OBC2", "dataframe": df, "entropy": None, "binding_energies": [-10.0, -12.0],
            "interaction_energies": [-25.0, -27.0], "physics_assumptions": []}


class _FailingDecomposition:
    mode = "raise"

    def __init__(self, *args, **kwargs):
        pass

    def run_per_residue_analysis(self, **kwargs):
        if self.mode == "raise":
            raise RuntimeError("decomposition platform unusable")
        return kwargs["baseline_results"]       # what decomposition.py returns when it fails


@pytest.fixture
def make_runner(tmp_path, monkeypatch):
    for name in ("complex.prmtop", "traj.xtc"):
        (tmp_path / name).write_text("placeholder\n")
    monkeypatch.setattr(runner_module, "PerResidueDecomposition", _FailingDecomposition)
    monkeypatch.setattr(runner_module.GBSACalculator, "run_comprehensive",
                        lambda self, **kwargs: _fake_mmgbsa(**kwargs))

    def make(strict=False, mode="raise"):
        _FailingDecomposition.mode = mode
        config = {"input_files": {"complex_pdb": str(tmp_path / "complex.prmtop"),
                                  "trajectory": str(tmp_path / "traj.xtc")},
                  "analysis_settings": {"run_per_residue_decomposition": True, "decomp_frames": 2,
                                        "max_frames": 2, "gb_model": "GBn2"},
                  "platform_settings": {"strict_platform": strict}}
        runner = MMGBSARunner(config, output_dir=str(tmp_path / f"run_{strict}_{mode}"))
        # the decomposition frames need the analysed frames of a real run
        monkeypatch.setattr(runner_module.GBSACalculator, "decomposition_frames",
                            lambda self, *a, **k: [0, 1])
        return runner
    return make


@pytest.mark.parametrize("mode", ["raise", "fallback"])
def test_failed_decomposition_is_marked_failed(make_runner, mode):
    runner = make_runner(mode=mode)
    results = runner.run_analysis()
    assert results is not None and results["decomposition"]["status"] == "failed"
    summary = yaml.safe_load((runner.output_dir / "results_summary.yaml").read_text())
    assert summary["decomposition"]["status"] == "failed"
    assert summary["decomposition"]["error"]
    report = (runner.output_dir / "final_report.txt").read_text()
    assert "FAILED" in report and "completed successfully" not in report


def test_failed_decomposition_with_strict_platform_stops(make_runner):
    runner = make_runner(strict=True)
    with pytest.raises(RuntimeError, match="strict_platform"):
        runner.run_analysis()
    # the results are written before the run stops
    summary = yaml.safe_load((runner.output_dir / "results_summary.yaml").read_text())
    assert summary["decomposition"]["status"] == "failed"


def test_cli_exits_non_zero_on_strict_failure(make_runner, tmp_path, monkeypatch):
    import sys
    from mmgbsa import cli
    runner = make_runner(strict=True)
    path = tmp_path / "strict.yaml"
    path.write_text(yaml.safe_dump(runner.config))
    monkeypatch.setattr(sys, "argv", ["opengbsa", str(path), "--output-dir", str(tmp_path / "cli")])
    assert cli.main() == 1


def test_analysis_config_is_the_effective_configuration(make_runner):
    runner = make_runner(mode="fallback")
    runner.run_analysis()
    text = (runner.output_dir / "analysis_config.yaml").read_text()
    assert text.startswith("# Effective configuration")
    written = yaml.safe_load(text)
    settings = written["analysis_settings"]
    assert settings["gb_model"] == "GBn2"                       # as given
    assert settings["temperature"] == 300.0                     # default filled in
    assert settings["max_workers"] == runner.calculator.max_workers >= 1   # resolved at run time
    assert settings["sa_model"] == "LCPO" and settings["lcpo_typing"] == "element"
    assert settings["solvent_dielectric"] == 78.5
    assert written["reporting_settings"]["export_width"] == 4200
    # a run can be repeated from it: it validates without warnings
    from mmgbsa.config import ConfigManager
    manager = ConfigManager(str(runner.output_dir / "analysis_config.yaml"))
    assert manager.validate_config(), manager.get_validation_errors()
    assert manager.get_validation_warnings() == []
