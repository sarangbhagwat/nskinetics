# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Progressive presentation frames of Extended Data Figure 4 (PI, ethanol
yield and IRR across the 1-D k_3 slice through the TRY-informed optimum).

The figure is drawn once by the unmodified
../plots/plot_axis_try_informed_sweep.py (plot_param('k_3'), exactly as
make_figure.py runs it), then edited as in the
manuscript: the bottom legend row is hidden (its blank canvas is kept) and a
rotated "growth derated at higher k_3 values" label is added beside the
dotted burden-derating line. The frames reveal it one y axis at a time; each
frame is saved on the script's own fixed canvas (no tight crop), so the
frames overlay exactly and the last one is the manuscript figure.

Artists are grouped by what the plot code made them: each curve is found by
the colour plot_param gave it (PI, ethanol yield, IRR, derate line) and its
Axes then names the y axis the curve belongs to. A y-axis group is that
Axes' y axis (ticks, tick labels, title) plus its side spines; the x-axis
group is the x axis plus the top / bottom spines of all three Axes. Any other
visible artist, or a legend entry other than the derating one, raises, so a
change in the plot code or the data cannot leak into frame 01 (or be hidden
with the legend) unnoticed.

Usage (through ../../make_figure_frames.py, which points the results root
NSKINETICS_IBO_RESULTS_DIR at the folder's ./data):
    python plot_k_3_try_informed_sweep_frames.py
           [--dpi N] [--formats png [pdf svg ...]]
Frames go to ../../output/frames/.
"""
import argparse
import contextlib
import importlib.util
import io
import os
import sys

sys.dont_write_bytecode = True   # no __pycache__ under code/ when run directly

from matplotlib.colors import same_color
from matplotlib.text import Text
from matplotlib.transforms import blended_transform_factory, offset_copy

HERE = os.path.dirname(os.path.abspath(__file__))
FOLDER = os.path.dirname(os.path.dirname(HERE))
OUT_DIR = os.path.join(FOLDER, 'output', 'frames')
STEM = 'k_3_try_informed_sweep_etoh_yield'

# Frame sequence: (frame name, [groups it reveals]); frames are cumulative.
# Groups: x-axis, yield-axis, yield-curve, pi-axis, pi-curve, irr-axis,
#         irr-curve, derate-line, derate-label
FRAMES = [
    ('axes', ['x-axis']),
    ('irr', ['irr-axis', 'irr-curve']),
    ('profitability-index', ['pi-axis', 'pi-curve']),
    ('ethanol-yield', ['yield-axis', 'yield-curve']),
    ('growth-derating', ['derate-line', 'derate-label']),
]

#: the label added in the manuscript beside the derate line, set to match it:
#: offset of the label's lower-left corner from the line's foot (points; the
#: text reads bottom to top, its first line next to the line) and font size.
#: The label follows the line's x, but its place right of the line and below
#: the curves was chosen for the manuscript's line at k_3 = 4.85: check the
#: last frame if the sweep data change.
DERATE_LABEL = 'growth derated at\nhigher $k_3$ values'
DERATE_LABEL_OFFSET_PT = (3.6, 9.3)
DERATE_LABEL_SIZE = 9.5


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pf = _load('_progressive_frames', os.path.join(HERE, 'progressive_frames.py'))


def draw_figure(dpi):
    """The complete figure from plot_param('k_3'), unsaved. plot_param's
    console summary (which ends in a 'wrote' line for files it does not
    write here) is silenced, and the empty output folders it creates under
    the results root are removed again."""
    sweep = _load('_axis_sweep_plot', os.path.join(
        FOLDER, 'code', 'plots', 'plot_axis_try_informed_sweep.py'))
    created = [d for d in (sweep.OUT_DIR, os.path.dirname(sweep.OUT_DIR))
               if not os.path.isdir(d)]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            _, figures = pf.capture_figures(sweep.plot_param, 'k_3', dpi)
    finally:
        for d in created:
            if os.path.isdir(d) and not os.listdir(d):
                os.rmdir(d)
    (fig,) = figures
    return sweep, fig


def curve(fig, color):
    """The one line of the figure drawn in `color` (a plot_param style)."""
    found = [line for ax in fig.axes for line in ax.get_lines()
             if same_color(line.get_color(), color)]
    if len(found) != 1:
        raise RuntimeError(f'expected one line of colour {color}, '
                           f'found {len(found)}')
    return found[0]


def y_axis(ax):
    """Everything that makes an Axes' y axis visible: the axis (ticks, tick
    labels, title) and its side spines."""
    return [ax.yaxis, ax.spines['left'], ax.spines['right']]


def add_derate_label(line):
    """The manuscript's rotated label, its lower-left corner tied to the
    derate line's x (data) and the axes bottom, offset in points."""
    ax = line.axes
    x = line.get_xdata()[0]
    base = blended_transform_factory(ax.transData, ax.transAxes)
    return ax.text(x, 0.0, DERATE_LABEL, rotation=90, rotation_mode='anchor',
                   ha='left', va='top', multialignment='left',
                   linespacing=1.25, fontsize=DERATE_LABEL_SIZE, color='black',
                   transform=offset_copy(base, fig=ax.figure,
                                         x=DERATE_LABEL_OFFSET_PT[0],
                                         y=DERATE_LABEL_OFFSET_PT[1],
                                         units='points'))


def groups_of(sweep, fig):
    """{group name: [artists]} of the finished figure, after the manuscript
    edits: the figure legend is hidden and the derate label added."""
    k13p = sweep.k13p
    pi, etoh, irr = (curve(fig, c) for c in
                     (sweep.PI_COLOR, sweep.YIELD_COLOR, sweep.IRR_COLOR))
    derate = curve(fig, k13p.DERATE_KW['color'])
    ax, ax_r, ax_irr = pi.axes, etoh.axes, irr.axes
    if len({ax, ax_r, ax_irr}) != 3 or len(fig.axes) != 3:
        raise RuntimeError('expected the PI, yield and IRR curves on three '
                           'separate Axes')
    for legend in fig.legends:   # cropped away in the manuscript
        entries = [text.get_text() for text in legend.get_texts()]
        if len(entries) != 1 or not entries[0].startswith('Growth derated'):
            raise RuntimeError('expected a legend with only the growth-'
                               f'derating entry (the label replaces it), '
                               f'got {entries}')
        legend.set_visible(False)
    groups = {
        'x-axis': [a.xaxis for a in fig.axes]
                  + [a.spines[s] for a in fig.axes for s in ('top', 'bottom')],
        'yield-axis': y_axis(ax_r), 'yield-curve': [etoh],
        'pi-axis': y_axis(ax), 'pi-curve': [pi],
        'irr-axis': y_axis(ax_irr), 'irr-curve': [irr],
        'derate-line': [derate],
        'derate-label': [add_derate_label(derate)],
    }
    grouped = {id(a) for artists in groups.values() for a in artists}
    for a in fig.axes:
        for child in a.get_children():
            if (child.get_visible() and id(child) not in grouped
                    and child is not a.patch      # the white axes background
                    and not (isinstance(child, Text) and not child.get_text())):
                raise RuntimeError(f'ungrouped visible artist {child!r}')
    for child in fig.get_children():   # nothing else drawn on the figure
        if (child.get_visible() and child is not fig.patch
                and child not in fig.axes):
            raise RuntimeError(f'ungrouped visible artist {child!r}')
    return groups


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--dpi', type=int, default=300,
                    help='resolution of the frames (default: %(default)s, '
                         "as make_figure.py's figure)")
    ap.add_argument('--formats', nargs='+', default=['png'],
                    help='file formats to write (default: png)')
    args = ap.parse_args()
    sweep, fig = draw_figure(args.dpi)
    pf.save_frames(fig, groups_of(sweep, fig), FRAMES, OUT_DIR, STEM,
                   args.dpi, formats=args.formats)


if __name__ == '__main__':
    main()
