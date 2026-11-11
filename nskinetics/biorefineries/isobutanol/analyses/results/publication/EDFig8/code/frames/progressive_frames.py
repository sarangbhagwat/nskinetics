# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Progressive frames of a finished matplotlib figure, for presentations.

A figure folder's make_figure_frames.py draws its figure ONCE with the
folder's unmodified plot script, sorts the figure's artists into named groups
and saves one image per frame, each frame showing the groups revealed so far.
Nothing is redrawn or re-laid-out between frames: an artist is only switched
visible, and every frame is saved on the same canvas, so the frames overlay
exactly on a slide and the last frame is the complete figure.

The same file is copied into every figure folder that makes frames (each
folder regenerates its figure from its own code/ only). Sim-safe: imports
only the standard library and matplotlib.

    capture_figures(func, ...)   run a plot function that saves and closes its
                                 own figure; get the figure back instead
    recording(ax, ...)           collect the artists a block adds to a container
    legend_entry(legend, label)  the handle + text artists of one legend row
    tight_bbox(fig, dpi)         the crop bbox_inches='tight' gives the
                                 complete figure, to reuse for every frame
    save_frames(fig, groups, frames, ...)   write the frames
"""
import contextlib
import glob
import os
import re

import matplotlib.pyplot as plt
from matplotlib.figure import Figure

__all__ = ('capture_figures', 'recording', 'legend_entry', 'tight_bbox',
           'save_frames')


def capture_figures(func, *args, **kwargs):
    """Call func(*args, **kwargs), a plot function that saves and closes the
    figure(s) it draws, with saving and closing switched off. Returns
    (func's return value, [the figures it created, still open]).

    Figure.savefig and pyplot.close are replaced only for the duration of the
    call, so the plot function runs unmodified and writes no file.
    """
    before = set(plt.get_fignums())
    savefig, close = Figure.savefig, plt.close
    Figure.savefig = lambda self, *a, **k: None
    plt.close = lambda *a, **k: None
    try:
        result = func(*args, **kwargs)
    finally:
        Figure.savefig, plt.close = savefig, close
    figures = [plt.figure(n) for n in plt.get_fignums() if n not in before]
    return result, figures


@contextlib.contextmanager
def recording(*containers):
    """Collect the artists added to the given containers (Axes or Figure:
    anything with get_children()) inside the with-block.

        with recording(ax) as new:
            ax.plot(x, y)
        groups['curve'] = new          # filled when the block exits
    """
    before = [{id(a) for a in c.get_children()} for c in containers]
    new = []
    yield new
    for c, seen in zip(containers, before):
        new.extend(a for a in c.get_children() if id(a) not in seen)


def legend_entry(legend, label):
    """[handle artist, text artist] of the legend row(s) labelled `label`, so
    a legend row can be revealed with the data it describes (the row keeps its
    place in the legend while hidden)."""
    found = []
    for handle, text in zip(legend.legend_handles, legend.get_texts()):
        if text.get_text() == label:
            found += [handle, text]
    if not found:
        raise KeyError(f'no legend entry labelled {label!r}; have '
                       f'{[t.get_text() for t in legend.get_texts()]}')
    return found


def tight_bbox(fig, dpi, pad_inches=None):
    """The bounding box [in] that savefig(dpi=dpi, bbox_inches='tight') crops
    the figure to AS IT IS NOW (call it while every artist is visible). Pass
    it as save_frames(bbox_inches=...) so all frames share the complete
    figure's crop; a per-frame 'tight' would crop each frame to what it shows.
    """
    if pad_inches is None:
        pad_inches = plt.rcParams['savefig.pad_inches']
    dpi0 = fig.get_dpi()
    fig.set_dpi(dpi)      # text extents are measured at the save resolution
    try:
        fig.canvas.draw()
        bbox = fig.get_tightbbox(fig.canvas.get_renderer())
    finally:
        fig.set_dpi(dpi0)
    return bbox.padded(pad_inches)


def _slug(name):
    return re.sub(r'[^0-9A-Za-z]+', '-', name).strip('-')


def save_frames(fig, groups, frames, out_dir, stem, dpi, formats=('png',),
                bbox_inches=None, **savefig_kwargs):
    """Save one image per frame of `fig` into out_dir; returns the paths.

    groups: {group name: [artists]}. Every artist of the figure that is in no
            group is part of the fixed background and is drawn in all frames.
    frames: [(frame name, [names of the groups this frame adds])], in order.
            Frames are cumulative: frame i shows the groups of frames 1..i.
            Every group must be added by exactly one frame, so the last frame
            is the complete figure.

    Files are <stem>_frame<NN>_<frame name>.<ext> (NN from 01); frames of
    `stem` left in out_dir by an earlier run are removed first. An artist that
    is already invisible in the finished figure stays invisible. The figure
    is left as it was (everything visible) on return.
    """
    added = [g for _, names in frames for g in names]
    unknown = sorted(set(added) - set(groups))
    unused = sorted(set(groups) - set(added))
    repeated = sorted({g for g in added if added.count(g) > 1})
    if unknown or unused or repeated:
        raise ValueError(f'frames do not partition the groups: unknown '
                         f'{unknown}, never shown {unused}, shown twice '
                         f'{repeated}')
    owner = {}
    for name, artists in groups.items():
        for a in artists:
            if owner.setdefault(id(a), name) != name:
                raise ValueError(f'{a!r} is in two groups: '
                                 f'{owner[id(a)]!r} and {name!r}')
    # only artists that are visible in the finished figure are toggled
    shown = {name: [a for a in artists if a.get_visible()]
             for name, artists in groups.items()}

    os.makedirs(out_dir, exist_ok=True)
    for ext in formats:
        for old in glob.glob(os.path.join(
                glob.escape(out_dir), f'{glob.escape(stem)}_frame[0-9][0-9]_*.{ext}')):
            os.remove(old)

    paths = []
    try:
        for artists in shown.values():
            for a in artists:
                a.set_visible(False)
        for i, (name, names) in enumerate(frames, start=1):
            for g in names:
                for a in shown[g]:
                    a.set_visible(True)
            for ext in formats:
                path = os.path.join(out_dir,
                                    f'{stem}_frame{i:02d}_{_slug(name)}.{ext}')
                fig.savefig(path, dpi=dpi, bbox_inches=bbox_inches,
                            **savefig_kwargs)
                paths.append(path)
                print('wrote', path)
    finally:
        for artists in shown.values():
            for a in artists:
                a.set_visible(True)
    return paths
