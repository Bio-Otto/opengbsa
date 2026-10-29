#!/usr/bin/env python3
"""
Frame selection (mmgbsa.frames): range, then stride, then the method; clamping
warns; the decomposition frames are a subset of the analysed frames.
"""
import random

import pytest

from mmgbsa import frames
from mmgbsa.frames import evenly_spaced, select_frames, subset_frames


def old_select(trajectory_length, max_frames=None, frame_start=None, frame_end=None,
               frame_stride=None, frame_selection='sequential', random_seed=42):
    """The selection of OpenGBSA <= 0.0.22 (GBSACalculator._select_frames)."""
    frame_start = 0 if frame_start is None else frame_start
    frame_end = trajectory_length if frame_end is None else frame_end
    frame_start = max(0, min(frame_start, trajectory_length - 1))
    frame_end = max(frame_start + 1, min(frame_end, trajectory_length))
    if frame_selection == "sequential":
        if frame_stride is None or frame_stride <= 0:
            idx = list(range(frame_start, frame_end))
        else:
            idx = list(range(frame_start, frame_end, frame_stride))
    elif frame_selection == "equidistant":
        if max_frames is None or max_frames <= 0:
            max_frames = frame_end - frame_start
        step = (frame_end - frame_start) / max_frames
        idx = [min(frame_start + int(i * step), frame_end - 1) for i in range(max_frames)]
    else:
        if max_frames is None or max_frames <= 0:
            max_frames = frame_end - frame_start
        random.seed(random_seed)
        available = list(range(frame_start, frame_end))
        if frame_stride is not None and frame_stride > 1:
            available = available[::frame_stride]
        idx = sorted(random.sample(available, min(max_frames, len(available))))
    if max_frames is not None and max_frames > 0:
        idx = idx[:max_frames]
    return sorted(set(idx))


@pytest.mark.parametrize("kwargs", [
    dict(trajectory_length=50, max_frames=50),
    dict(trajectory_length=1000, max_frames=50, frame_start=100, frame_end=500, frame_stride=5),
    dict(trajectory_length=1000, max_frames=30, frame_start=50, frame_selection='equidistant'),
    dict(trajectory_length=1000, max_frames=None, frame_selection='equidistant'),
    dict(trajectory_length=37, max_frames=60, frame_selection='equidistant'),
    dict(trajectory_length=1000, max_frames=25, frame_start=200, frame_end=800, frame_stride=2,
         frame_selection='random', random_seed=123),
    dict(trajectory_length=100, max_frames=None, frame_selection='random', random_seed=7),
    dict(trajectory_length=2000, max_frames=100, frame_start=1000, frame_stride=5),
])
def test_unchanged_where_the_old_rule_was_right(kwargs):
    """Without frame_stride for equidistant, the frames are those of 0.0.22."""
    n = kwargs.pop('trajectory_length')
    assert select_frames(n, **kwargs) == old_select(n, **kwargs)


def test_equidistant_applies_the_stride():
    chosen = select_frames(1000, max_frames=10, frame_stride=10, frame_selection='equidistant')
    assert chosen == list(range(0, 1000, 100))
    assert all((f - 5) % 10 == 0 for f in select_frames(1000, 7, 5, 905, 10, 'equidistant'))
    # the old rule ignored the stride
    assert old_select(1000, 10, frame_stride=10, frame_selection='equidistant') == list(range(0, 1000, 100))
    assert old_select(1000, 30, 0, 1000, 7, 'equidistant') != select_frames(1000, 30, 0, 1000, 7, 'equidistant')


def test_range_then_stride_then_method():
    assert select_frames(100, 3, 10, 50, 4) == [10, 14, 18]
    assert select_frames(100, None, 10, 50, 4) == list(range(10, 50, 4))
    rnd = select_frames(100, 5, 10, 50, 4, 'random', 3)
    assert len(rnd) == 5 and set(rnd) <= set(range(10, 50, 4))


def test_clamping_warns():
    messages = []
    assert select_frames(10, None, 50, None, warn=messages.append) == [9]
    assert "frame_start 50 is beyond the last frame" in messages[0]
    messages.clear()
    assert select_frames(10, None, 2, 99, warn=messages.append) == list(range(2, 10))
    assert "frame_end 99 is beyond the trajectory (10 frames)" in messages[0]
    messages.clear()
    select_frames(10, None, None, None, warn=messages.append)
    assert messages == []


def test_unknown_method_is_an_error():
    with pytest.raises(ValueError, match="Unknown frame selection"):
        select_frames(10, frame_selection='every_other')


def test_random_default_seed_is_42():
    assert select_frames(200, 10, frame_selection='random', random_seed=None) == \
        select_frames(200, 10, frame_selection='random', random_seed=42)


ANALYSED = select_frames(1000, 30, frame_selection='equidistant')    # 0, 33, 66, ...


@pytest.mark.parametrize("method", ['sequential', 'equidistant', 'random'])
@pytest.mark.parametrize("n_select", [1, 7, 29, 30, 31, None])
def test_decomposition_frames_are_a_subset_of_the_analysed_frames(method, n_select):
    chosen = subset_frames(ANALYSED, n_select, method, 11)
    assert set(chosen) <= set(ANALYSED)
    assert chosen == sorted(chosen)
    expected = len(ANALYSED) if n_select is None or n_select >= len(ANALYSED) else n_select
    assert len(chosen) == expected or method == 'equidistant'


def test_decomposition_frame_rules():
    analysed = list(range(100, 200, 5))                      # 20 frames
    assert subset_frames(analysed, 4, 'sequential') == [100, 105, 110, 115]
    # same rounding as the normal-mode frames: first and last included
    assert subset_frames(analysed, 4, 'equidistant') == [analysed[i] for i in (0, 6, 13, 19)]
    assert evenly_spaced(20, 4) == [0, 6, 13, 19]
    r1, r2 = subset_frames(analysed, 6, 'random', 5), subset_frames(analysed, 6, 'random', 5)
    assert r1 == r2 and len(r1) == 6
    # 0.0.22 re-selected equidistant frames over the whole range: not a subset
    old = old_select(1000, 7, frame_selection='equidistant')
    assert not set(old) <= set(ANALYSED)


def test_nmode_uses_the_same_even_spacing():
    from mmgbsa import nmode_entropy
    for n, k in [(50, 10), (7, 3), (100, 1), (5, 5), (5, None)]:
        assert nmode_entropy.select_frames(n, k) == evenly_spaced(n, k)


def test_calculator_decomposition_frames():
    from mmgbsa.mmgbsa_core import GBSACalculator
    calc = GBSACalculator()
    with pytest.raises(RuntimeError):
        calc.decomposition_frames(5)
    calc.analysed_frames = [3, 6, 9, 12]
    assert calc.decomposition_frames(2, 'sequential') == [3, 6]
    assert calc.decomposition_frames(10, 'random') == [3, 6, 9, 12]
    assert calc._select_frames(20, 3, 4, None, 2, 'sequential') == [4, 6, 8]
