#!/usr/bin/env python3
"""
The HTML report draws its time series over the trajectory frames that were
decomposed, matched to the MM/GBSA rows by trajectory index (before, row
positions were used, which mismatched with frame_start > 0 or equidistant/
random selection), and plot_top_residues sets the number of residues in its
bar chart, heat map and hot-spot table.
"""
import json
import re

import numpy as np
import pandas as pd

from mmgbsa.reporting import HTMLReportGenerator


def _frame_data(frame_index):
    rows = []
    for k, f in enumerate(frame_index):
        row = {"frame_index": f, "frame_number": k + 1}
        for r in range(6):
            row[f"RES{r}_total"] = -float(r) - 0.01 * k
        rows.append(row)
    return rows


def _global(trajectory_frames):
    n = len(trajectory_frames)
    return {"dataframe": pd.DataFrame({
        "frame": np.arange(n), "trajectory_frame": trajectory_frames,
        "binding_energy": [-100.0 - f for f in trajectory_frames],
        "delta_vdw": [-1.0 * f for f in trajectory_frames], "delta_elec": [0.0] * n,
        "delta_gb": [0.0] * n, "delta_sa": [0.0] * n}), "entropy": None}


def _plot_x_values(html):
    """x arrays of the Plotly traces in the HTML."""
    return [json.loads(m) for m in re.findall(r'"x":(\[[^\]]*\])', html)]


def test_trajectory_frames_from_the_decomposition_rows():
    df = pd.DataFrame(_frame_data([110, 120, 130]))
    assert HTMLReportGenerator._trajectory_frames(df, [110, 120, 130]) == [110, 120, 130]
    # parallel decomposition: frame_index is the position among the decomposed frames
    df = pd.DataFrame(_frame_data([0, 1, 2]))
    assert HTMLReportGenerator._trajectory_frames(df, [110, 120, 130]) == [110, 120, 130]
    assert HTMLReportGenerator._trajectory_frames(df, None) == [0, 1, 2]


def test_time_series_match_the_mmgbsa_rows_by_trajectory_frame(tmp_path):
    report = HTMLReportGenerator(str(tmp_path), config={"analysis_settings": {"plot_top_residues": 3}})
    analysed = list(range(100, 150, 5))                     # frame_start 100, stride 5
    decomposed = [110, 125, 145]
    html = report._generate_plotly_div(_frame_data(decomposed), global_results=_global(analysed),
                                       decomposition_frames=decomposed)
    assert "Plot Generation Failed" not in html
    xs = _plot_x_values(html)
    assert decomposed in xs                                  # time series and heat map x axis
    # the total binding series is that of the decomposed frames, not of rows 0-2
    ys = [json.loads(m) for m in re.findall(r'"y":(\[[^\]]*\])', html)]
    assert [-210.0, -225.0, -245.0] in ys
    assert "Top 3 Contributing Residues" in html


def test_unmatched_frames_leave_out_the_mmgbsa_series(tmp_path):
    report = HTMLReportGenerator(str(tmp_path), config={})
    glob = _global(list(range(10)))
    glob["dataframe"] = glob["dataframe"].drop(columns="trajectory_frame")
    html = report._generate_plotly_div(_frame_data([3, 5]), global_results=glob)
    ys = [json.loads(m) for m in re.findall(r'"y":(\[[^\]]*\])', html)]
    assert [-100.0, -101.0] not in ys                        # not the first two rows


def test_hotspot_table_uses_plot_top_residues(tmp_path):
    report = HTMLReportGenerator(str(tmp_path), config={"analysis_settings": {"plot_top_residues": 4}})
    residues = pd.DataFrame({"residue_id": [f"R{i}" for i in range(12)], "total": -np.arange(12.0),
                             "vdw": -np.arange(12.0), "electrostatic": np.zeros(12),
                             "solvation": np.zeros(12)})
    html = report._generate_stats_html({"dataframe": residues}, global_results=None)
    assert "Top 4 Binding Hotspots" in html
    assert html.count("<td>R") == 4 or len(re.findall(r">R\d+<", html)) == 4
