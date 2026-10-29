"""
Which trajectory frames a run analyses, and which of them it decomposes.

`select_frames` applies the options in this order: the range
[frame_start, frame_end), then frame_stride, then frame_selection over the
strided frames (sequential: the first max_frames; equidistant: max_frames
spread evenly; random: max_frames chosen with random_seed). Out-of-range
frame_start/frame_end are clamped with a warning.

`subset_frames` chooses the decomp_frames frames of the per-residue
decomposition among the analysed frames, so they are always a subset of
them: sequential -> the first ones, equidistant -> evenly spaced with the
rounding rule of the normal-mode frames (nmode_entropy.select_frames),
random -> a seeded random subset.
"""
import random

import numpy as np

from .logger import ToolLogger

log = ToolLogger()

METHODS = ('sequential', 'equidistant', 'random')


def evenly_spaced(n_frames, n_select):
    """Positions of `n_select` items spread evenly over `n_frames` (first and last
    included; rounding may give fewer). None or n_select >= n_frames: all."""
    if n_select is None or n_select >= n_frames:
        return list(range(n_frames))
    return sorted({int(round(x)) for x in np.linspace(0, n_frames - 1, n_select)})


def frame_range(n_total, frame_start=None, frame_end=None, warn=None):
    """[start, end) of the trajectory after clamping, with a warning for each clamp."""
    warn = warn or log.warning
    if n_total < 1:
        raise ValueError("The trajectory has no frames")
    start = 0 if frame_start is None else int(frame_start)
    end = n_total if frame_end is None else int(frame_end)
    if start < 0:
        warn(f"frame_start {start} is negative; using 0.")
        start = 0
    if start > n_total - 1:
        warn(f"frame_start {start} is beyond the last frame of the trajectory ({n_total - 1}); "
             f"using frame {n_total - 1} only.")
        start = n_total - 1
    if end > n_total:
        warn(f"frame_end {end} is beyond the trajectory ({n_total} frames); using {n_total}.")
        end = n_total
    if end <= start:
        warn(f"frame_end {end} is not after frame_start {start}; using frame {start} only.")
        end = start + 1
    return start, end


def select_frames(n_total, max_frames=None, frame_start=None, frame_end=None, frame_stride=None,
                  frame_selection='sequential', random_seed=42, warn=None):
    """Trajectory frame indices to analyse (sorted, unique).

    n_total : number of frames of the trajectory
    max_frames : how many frames to choose (None or <= 0: every frame of the
        range after striding)
    frame_stride : every Nth frame of the range (None or <= 0: 1)
    """
    if frame_selection not in METHODS:
        raise ValueError(f"Unknown frame selection method: {frame_selection} "
                         f"(one of {list(METHODS)})")
    start, end = frame_range(n_total, frame_start, frame_end, warn)
    stride = int(frame_stride) if frame_stride and frame_stride > 0 else 1
    candidates = list(range(start, end, stride))
    n = len(candidates) if (max_frames is None or max_frames <= 0) else min(int(max_frames), len(candidates))

    if frame_selection == 'sequential':
        chosen = candidates[:n]
    elif frame_selection == 'equidistant':
        # n frames spread over the strided frames; the first is always included
        # (with stride 1 this is the rule of earlier versions)
        if n >= len(candidates):
            chosen = candidates
        else:
            step = len(candidates) / n
            chosen = [candidates[int(i * step)] for i in range(n)]
    else:
        seed = 42 if random_seed is None else random_seed
        chosen = random.Random(seed).sample(candidates, n)
    return sorted(set(chosen))


def subset_frames(analysed_frames, n_select, frame_selection='sequential', random_seed=42):
    """The frames of the per-residue decomposition: `n_select` of `analysed_frames`.

    analysed_frames : trajectory indices analysed by the run (select_frames)
    n_select : decomp_frames (None, or more than the analysed frames: all of them)
    """
    analysed = sorted(int(f) for f in analysed_frames)
    if n_select is None or n_select >= len(analysed):
        return analysed
    n_select = int(n_select)
    if frame_selection == 'sequential':
        return analysed[:n_select]
    if frame_selection == 'equidistant':
        return [analysed[i] for i in evenly_spaced(len(analysed), n_select)]
    if frame_selection == 'random':
        seed = 42 if random_seed is None else random_seed
        return sorted(random.Random(seed).sample(analysed, n_select))
    raise ValueError(f"Unknown frame selection method: {frame_selection} (one of {list(METHODS)})")
