#!/usr/bin/env python3
"""
Run the third-party validation cases and compare the binding enthalpies with
the expected values in expected_results.yaml.

    python test/configs/check_validation.py --list
    python test/configs/check_validation.py 3dd2
    python test/configs/check_validation.py 3dd2 --amber-compatible
    python test/configs/check_validation.py all

Each case runs `opengbsa <config>` in its own folder (the configurations use
paths relative to it), after its data were downloaded with that folder's
./fetch_data.sh. At default settings the mean binding enthalpy is compared
with the value recorded for this configuration. With --amber-compatible the
run uses amber_compatible: true and is compared with Amber MMPBSA.py on the
same frames. A case passes when the difference is within the tolerance.
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent


def load_cases():
    with open(HERE / "expected_results.yaml") as fh:
        return yaml.safe_load(fh)


def inputs_present(folder, config):
    files = config.get("input_files", {})
    missing = [v for k, v in files.items()
               if v and k not in ("ligand_resname",) and not (folder / v).exists()]
    for v in config.get("forcefield_settings", {}).get("charmm_params") or []:
        if not (folder / v).exists():
            missing.append(v)
    return missing


def run_case(name, case, amber_compatible, opengbsa, keep):
    folder = HERE / case["folder"]
    with open(folder / case["config"]) as fh:
        config = yaml.safe_load(fh)
    missing = inputs_present(folder, config)
    if missing:
        print(f"  {name}: input data missing ({missing[0]} ...). Prepare them first:\n"
              f"      cd {HERE / case['folder']} && {case['prepare']}")
        return None

    config["analysis_settings"]["amber_compatible"] = bool(amber_compatible)
    out = Path(tempfile.mkdtemp(prefix=f"opengbsa_check_{name}_"))
    run_config = folder / f".check_{name}.yaml"   # in the case folder: its paths are relative
    with open(run_config, "w") as fh:
        yaml.safe_dump(config, fh, sort_keys=False)
    try:
        log = out / "run.log"
        with open(log, "w") as fh:
            done = subprocess.run([opengbsa, run_config.name, "--output-dir", str(out)],
                                  cwd=folder, stdout=fh, stderr=subprocess.STDOUT)
        summaries = sorted(out.glob("*/results_summary.yaml"))
        if done.returncode != 0 or not summaries:
            print(f"  {name}: the run failed; see {log}")
            return False
        with open(summaries[-1]) as fh:
            mmgbsa = yaml.safe_load(fh)["mmgbsa"]
        value = float(mmgbsa["mean_binding_energy"])
    finally:
        run_config.unlink(missing_ok=True)
        if not keep:
            shutil.rmtree(out, ignore_errors=True)

    if amber_compatible:
        reference, label = case.get("mmpbsa"), "Amber MMPBSA.py"
    else:
        reference, label = case["opengbsa_defaults"], "expected"
    if reference is None:
        print(f"  {name}: dH = {value:.3f} kcal/mol ({case.get('note', 'no reference value')})")
        return True
    diff = value - reference
    ok = abs(diff) <= case.get("tolerance", 0.02)
    print(f"  {name}: dH = {value:9.3f}   {label} {reference:9.3f}   "
          f"difference {diff:+.3f} kcal/mol   {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    cases = load_cases()
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("cases", nargs="*", help="case names, or 'all'")
    parser.add_argument("--list", action="store_true", help="list the cases and exit")
    parser.add_argument("--amber-compatible", action="store_true",
                        help="run with amber_compatible: true and compare with Amber MMPBSA.py")
    parser.add_argument("--opengbsa", default=shutil.which("opengbsa") or "opengbsa",
                        help="the opengbsa executable (default: the one on PATH)")
    parser.add_argument("--keep", action="store_true", help="keep the output folders")
    args = parser.parse_args()

    if args.list or not args.cases:
        print(f"{'case':16s} {'folder':22s} {'frames':>6s}  {'download':>9s}  description")
        for name, case in cases.items():
            print(f"{name:16s} {case['folder']:22s} {case['frames']:6d}  {case['download']:>9s}  {case['description']}")
        return 0

    names = list(cases) if args.cases == ["all"] else args.cases
    unknown = [n for n in names if n not in cases]
    if unknown:
        print(f"Unknown case(s): {unknown}. See --list.")
        return 2
    setting = "amber_compatible" if args.amber_compatible else "default settings"
    print(f"OpenGBSA validation cases, {setting} ({args.opengbsa}):")
    results = [run_case(n, cases[n], args.amber_compatible, args.opengbsa, args.keep) for n in names]
    ran = [r for r in results if r is not None]
    print(f"{sum(ran)} of {len(ran)} cases passed" + (f"; {len(results) - len(ran)} skipped (no data)"
                                                   if len(ran) < len(results) else ""))
    return 0 if ran and all(ran) else 1


if __name__ == "__main__":
    sys.exit(main())
