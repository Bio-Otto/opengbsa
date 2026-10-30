"""How much each option of the energy function and the frame choice changes the binding enthalpy.

    python docs/scripts/option_effects.py [--platform CUDA] [--work DIR]

Runs two systems of the repository once with their base settings and once
for every variant below, each changing one option:

* 6T1H with ligand 1151 (examples/quickstart, Native Mode, 50 frames), the
  quickstart settings;
* IRF11 with DNA (test/configs/dna_protein_test, Coordinate Mode, 101 frames),
  a highly charged protein-nucleic acid complex; run its fetch_data.sh and
  prepare_system.py first.

Writes docs/data/option_effects.csv: the mean binding enthalpy, its standard
deviation over the frames and the change from the base settings, in kcal/mol.
From it, the two tables the user guide shows: option_effects_energy.csv (the
energy function) and option_effects_frames.csv (the frame choice).
``--tables-only`` rebuilds these two from option_effects.csv without running.
Run it from the root of the repository.
"""
import argparse
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "data" / "option_effects.csv"
DNA = "test/configs/dna_protein_test/prepared"

SYSTEMS = {
    "6T1H-1151": {
        "input_files": {"complex_pdb": "test/data/6t1h_1151_comp/complex.prmtop",
                        "receptor_topology": "test/data/6t1h_1151_comp/receptor.prmtop",
                        "ligand_topology": "test/data/6t1h_1151_comp/ligand.prmtop",
                        "trajectory": "examples/quickstart/trajectory.xtc"},
        "analysis_settings": {"gb_model": "OBC2", "salt_concentration": 0.15, "temperature": 300.0,
                              "max_frames": 50},
    },
    "IRF11-DNA": {
        "input_files": {"complex_pdb": f"{DNA}/complex_ref.pdb", "trajectory": f"{DNA}/complex_traj.dcd"},
        "analysis_settings": {"gb_model": "OBC2", "salt_concentration": 0.15, "temperature": 300.0,
                              "max_frames": 101, "receptor_selection": "chainid 0",
                              "ligand_selection": "chainid 1", "reimage_trajectory": False},
        "forcefield_settings": {"protein_forcefield": "amber14"},
    },
}

# (option, value[, other options set with it]); each changes one analysis_settings option
VARIANTS = [
    ("gb_model", "HCT"), ("gb_model", "OBC1"), ("gb_model", "GBn"), ("gb_model", "GBn2"),
    ("gb_radii_set", "mbondi"), ("gb_radii_set", "mbondi2"), ("gb_radii_set", "mbondi3"),
    ("gb_radii_set", "bondi"),
    ("salt_concentration", 0.0), ("salt_concentration", 0.5),
    ("temperature", 310.0),
    ("solute_dielectric", 2.0), ("solute_dielectric", 4.0),
    ("solvent_dielectric", 80.0),
    ("sa_model", "ACE"),
    ("lcpo_typing", "amber"),
    ("rgbmax", 25.0), ("rgbmax", 15.0),
    ("amber_compatible", True),
    ("max_frames", 10),
    ("frame_selection", "equidistant", {"max_frames": 10}),
    ("frame_selection", "random", {"max_frames": 10}),
    ("frame_stride", 5),
]


def run(system, settings, platform, work):
    base = SYSTEMS[system]
    cfg = {**base, "analysis_settings": {**base["analysis_settings"], **settings},
           "platform_settings": {"preferred_platform": platform, "strict_platform": True}}
    out = Path(tempfile.mkdtemp(dir=work))
    (out / "config.yaml").write_text(yaml.safe_dump(cfg))
    # the opengbsa of the running interpreter's environment
    exe = Path(sys.executable).with_name("opengbsa")
    subprocess.run([str(exe), str(out / "config.yaml"), "--output-dir", str(out)], cwd=ROOT,
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    df = pd.read_csv(next(out.glob("*analysis_*/fixed_enhanced_mmgbsa_results_*.csv")))
    return df["binding_energy"].mean(), df["binding_energy"].std()


FRAME_OPTIONS = ("max_frames", "frame_selection", "frame_stride")


def write_tables():
    """The two tables of the user guide, from option_effects.csv."""
    t = pd.read_csv(OUT).fillna("")
    cols = {"Option": t["option"].map(lambda o: o if o.startswith("(") else f"``{o}``"),
            "Value": t["value"]}
    for system in SYSTEMS:
        cols[f"{system} ΔH"] = [f"{h:.2f} ± {s:.2f}" for h, s in zip(t[f"{system} dH"], t[f"{system} SD"])]
        cols[f"{system} change"] = [f"{c:+.2f}" if i else "" for i, c in enumerate(t[f"{system} change"])]
    shown = pd.DataFrame(cols)
    frames = t["option"].isin(FRAME_OPTIONS)
    shown[~frames].to_csv(OUT.with_name("option_effects_energy.csv"), index=False)
    shown[frames | (t.index == 0)].to_csv(OUT.with_name("option_effects_frames.csv"), index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--platform", default="CUDA")
    ap.add_argument("--work", default=None, help="directory for the runs (default: a temporary one)")
    ap.add_argument("--tables-only", action="store_true", help="rebuild the two tables from option_effects.csv")
    a = ap.parse_args()
    if a.tables_only:
        write_tables()
        return
    work = Path(a.work or tempfile.mkdtemp())
    work.mkdir(parents=True, exist_ok=True)

    rows = [{"option": "(base settings)", "value": ""}]
    for v in VARIANTS:
        option, value, extra = (*v, {})[:3]
        label = str(value).lower() if isinstance(value, bool) else str(value)
        if extra:
            label += ", " + ", ".join(f"{k} {x}" for k, x in extra.items())
        rows.append({"option": option, "value": label, "settings": {option: value, **extra}})
    for system in SYSTEMS:
        t0 = time.time()
        base, base_sd = run(system, {}, a.platform, work)
        rows[0].update({f"{system} dH": base, f"{system} SD": base_sd, f"{system} change": 0.0})
        for row in rows[1:]:
            dh, sd = run(system, row["settings"], a.platform, work)
            row.update({f"{system} dH": dh, f"{system} SD": sd, f"{system} change": dh - base})
            print(f"{system}: {row['option']} = {row['value']}: dH {dh:.2f} kcal/mol", file=sys.stderr)
        print(f"{system}: {time.time() - t0:.0f} s", file=sys.stderr)
    t = pd.DataFrame(rows).drop(columns="settings")
    t.round(2).to_csv(OUT, index=False)
    write_tables()
    print(f"wrote {OUT} and the two tables of the user guide", file=sys.stderr)


if __name__ == "__main__":
    main()
