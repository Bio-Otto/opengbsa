"""Figures of the OpenGBSA documentation.

    python docs/scripts/make_figures.py

Writes PNG files to docs/source/_static/figures/. The diagrams are drawn here;
the plots read the small tables in docs/data/ (validation_summary.csv: mean
binding enthalpies of the validation runs of the OpenGBSA paper, OpenGBSA at
default settings and with amber_compatible against Amber MMPBSA.py;
entropy_summary.csv: interaction-energy entropy estimates of the same runs).
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "source" / "_static" / "figures"
DATA = ROOT / "data"

INK, MUTED, LINE = "#1f2933", "#52606d", "#9aa5b1"
BLUE, ORANGE, GREEN, VERMILION, PURPLE = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#6F4E9C"
FILL = {"input": "#e8f1fb", "process": "#eef2f5", "native": "#e3f4ee", "coord": "#fff4e0",
        "output": "#fbeee6", "config": "#f1ecf8"}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": LINE,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False})


def box(ax, x, y, w, h, text, fill, size=10, weight="normal", ha="center"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                fc=fill, ec=LINE, lw=1.0))
    tx = x + w / 2 if ha == "center" else x + 0.12
    ax.text(tx, y + h / 2, text, ha=ha, va="center", fontsize=size, color=INK, weight=weight,
            linespacing=1.35)


def arrow(ax, x0, y0, x1, y1, text=None, color=MUTED, rad=0.0):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=12,
                                 color=color, lw=1.2, connectionstyle=f"arc3,rad={rad}"))
    if text:
        ax.text((x0 + x1) / 2, (y0 + y1) / 2 + 0.12, text, ha="center", va="bottom",
                fontsize=9, color=MUTED)


def canvas(w, h):
    fig, ax = plt.subplots(figsize=(w, h))
    ax.set_xlim(0, w)
    ax.set_ylim(0, h)
    ax.axis("off")
    return fig, ax


def save(fig, name):
    fig.savefig(OUT / name, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", OUT / name)


def thermodynamic_cycle():
    fig, ax = canvas(9.0, 3.8)
    box(ax, 0.3, 2.6, 2.6, 1.0, "complex\nin water (GB)", FILL["input"], weight="bold")
    box(ax, 5.9, 2.6, 2.8, 1.0, "receptor + ligand\nin water (GB)", FILL["input"], weight="bold")
    box(ax, 0.3, 0.4, 2.6, 1.0, "complex\nin the gas phase", FILL["process"])
    box(ax, 5.9, 0.4, 2.8, 1.0, "receptor + ligand\nin the gas phase", FILL["process"])
    arrow(ax, 5.85, 3.1, 2.95, 3.1)
    ax.text(4.4, 3.25, r"$\Delta G_\mathrm{bind}$", ha="center", fontsize=12, color=INK)
    arrow(ax, 5.85, 0.9, 2.95, 0.9)
    ax.text(4.4, 1.05, r"$\Delta E_\mathrm{MM} = \Delta E_\mathrm{vdW} + \Delta E_\mathrm{elec}$",
            ha="center", fontsize=10, color=INK)
    arrow(ax, 1.6, 1.45, 1.6, 2.55)
    ax.text(1.75, 2.0, r"$G_\mathrm{solv}$(complex)", fontsize=9, color=MUTED, va="center")
    arrow(ax, 7.3, 1.45, 7.3, 2.55)
    ax.text(7.15, 2.0, r"$G_\mathrm{solv}$(receptor) + $G_\mathrm{solv}$(ligand)", fontsize=9,
            color=MUTED, va="center", ha="right")
    ax.text(4.5, -0.25,
            r"$\Delta H_\mathrm{bind} = \Delta E_\mathrm{vdW} + \Delta E_\mathrm{elec} + "
            r"\Delta G_\mathrm{GB} + \Delta G_\mathrm{SA}$   (single trajectory: bonded terms cancel)",
            ha="center", fontsize=10, color=INK)
    save(fig, "thermodynamic_cycle.png")


def interfaces():
    fig, ax = canvas(10.6, 4.6)
    box(ax, 0.2, 3.3, 3.2, 1.0, "CLI\ncommand-line interface", FILL["input"], weight="bold")
    box(ax, 0.2, 0.3, 3.2, 2.7,
        "opengbsa --create-config\nopengbsa --create-complete-config\nopengbsa --check-deps\n"
        "opengbsa --list-test-data\nopengbsa --fetch-test-data\nopengbsa --version",
        "white", size=9, ha="left")
    box(ax, 4.0, 3.3, 3.5, 1.0, "CDI\nconfiguration-driven interface", FILL["config"], weight="bold")
    box(ax, 4.0, 0.3, 3.5, 2.7,
        "config.yaml\n  input_files\n  output_settings\n  analysis_settings\n  forcefield_settings\n"
        "  platform_settings\n  advanced_settings",
        "white", size=9, ha="left")
    box(ax, 8.2, 2.1, 2.3, 1.0, "opengbsa config.yaml\n[--output-dir DIR]", FILL["process"], size=9)
    box(ax, 8.2, 0.3, 2.3, 1.3, "results folder\nCSV, YAML, text,\nplots, HTML, PyMOL", FILL["output"], size=9)
    arrow(ax, 3.45, 2.5, 3.95, 2.5, "writes")
    arrow(ax, 7.55, 2.6, 8.15, 2.6, "read by")
    arrow(ax, 9.35, 2.05, 9.35, 1.65)
    save(fig, "interfaces.png")


def input_modes():
    fig, ax = canvas(10.0, 5.0)
    box(ax, 3.4, 4.0, 3.2, 0.8, "Is the MD topology available?", FILL["process"], weight="bold")
    box(ax, 0.2, 1.6, 4.6, 1.8,
        "Native Mode\nthe force field of the simulation is read as is\n\n"
        ".prmtop (Amber)   .top / .tpr (GROMACS)\n.psf + parameter files (CHARMM/NAMD)\n"
        ".xml + .pdb (OpenMM System)", FILL["native"], size=9)
    box(ax, 5.2, 1.6, 4.6, 1.8,
        "Coordinate Mode\nOpenGBSA builds the topology\n\n"
        ".pdb or .gro of the complex + trajectory\nproteins and nucleic acids: amber14-all\n"
        "small-molecule ligand: .sdf, OpenFF", FILL["coord"], size=9)
    arrow(ax, 4.2, 3.95, 2.5, 3.45, "yes")
    arrow(ax, 5.8, 3.95, 7.5, 3.45, "no")
    box(ax, 0.2, 0.1, 9.6, 1.0,
        "Receptor and ligand: ligand_selection / receptor_selection (MDTraj), or binding_mode: ppi "
        "(chain 0 / chain 1),\nor receptor_topology + ligand_topology (Amber), "
        "or the ligand residue name (LIG, UNL, UNK, DRUG, MOL)", "white", size=9)
    save(fig, "input_modes.png")


def validation_agreement():
    d = pd.read_csv(DATA / "validation_summary.csv")
    d = d[~d.case.str.startswith("1gcq_") | (d.case == "1gcq_obc2")]
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.axhline(0, color=LINE, lw=1)
    ax.scatter(d.amber, d.preset - d.amber, s=26, color=BLUE, zorder=3,
               label="amber_compatible: true")
    ax.scatter(d.amber, d.defaults - d.amber, s=26, color=ORANGE, marker="s", zorder=3,
               label="default settings")
    ax.set_xlabel(r"Amber MMPBSA.py $\Delta H_\mathrm{bind}$ (kcal mol$^{-1}$)")
    ax.set_ylabel(r"OpenGBSA $-$ MMPBSA.py (kcal mol$^{-1}$)")
    ax.text(ax.get_xlim()[0] + 1, 0.4, "identical result", fontsize=8, color=MUTED, va="bottom")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.14), ncol=2)
    save(fig, "validation_agreement.png")


def entropy_reliability():
    e = pd.read_csv(DATA / "entropy_summary.csv")
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.scatter(e.sigma, e.ie, s=22, color=BLUE, label="interaction entropy (IE)", zorder=3)
    ax.scatter(e.sigma, e.c2, s=22, color=VERMILION, marker="s", label="C2", zorder=3)
    for limit, text in ((15 / 4.184, "IE limit 3.6"), (25 / 4.184, "C2 limit 6.0")):
        ax.axvline(limit, color=LINE, lw=1, ls="--")
        ax.text(limit * 1.05, 2e4, text, fontsize=8, color=MUTED, rotation=90, va="top")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\sigma$ of the interaction energy (kcal mol$^{-1}$)")
    ax.set_ylabel(r"$-T\Delta S$ (kcal mol$^{-1}$)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.14), ncol=2)
    save(fig, "entropy_reliability.png")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    thermodynamic_cycle()
    interfaces()
    input_modes()
    validation_agreement()
    entropy_reliability()
