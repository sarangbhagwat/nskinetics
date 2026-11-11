# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Progressive presentation frames of Extended Data Figure 8 (k_13 x
isobutanol-inhibition IRR contour, opt_PI_TRY_informed).

Loads plots/plot_k13_inhib_isobutanol_opt_PI_TRY_informed.py by file path and
draws the complete figure once with its engine's make_figure, exactly as
make_figure.py does. While the panel is drawn, the axes' contourf / contour /
plot / annotate calls and the colourbar's axhline are recorded, so each
artist is filed under what the plot code drew it for:

    IRR map                     the filled IRR contours
    break-even contour          the white IRR = 0 contour on the map
    colour bar                  the colourbar axes (bar, ticks, title)
    colour bar break-even line  the white 0 line on the colourbar
    <optimum label>             one optimum's marker + its callout (text and
                                arrow), one group per row of the plot
                                script's OPTIMA, named by its label with line
                                breaks as spaces: 'profitability', 'TCI',
                                'cell density', 'ethanol titer', 'ethanol
                                yield', 'ethanol productivity', 'isobutanol
                                titer', 'isobutanol yield', 'isobutanol
                                productivity'

Everything else (the plot axes' spines, ticks, tick labels and axis titles)
is the background shown in every frame. Frames are cumulative and saved with
one fixed crop (that of the complete figure), so they overlay exactly and the
last frame is the figure make_figure.py saves. Re-cut the sequence by editing
FRAMES only.

Usage: python plot_k13_inhib_isobutanol_opt_PI_TRY_informed_frames.py
           [--dpi N] [--formats png [pdf svg ...]]
(run through ../../make_figure_frames.py, which sets the results directory).
Frames go to ../../output/frames/.

Sim-safe: reads the sweep CSVs only and never imports the package.
"""

import argparse
import contextlib
import importlib.util
import os

from matplotlib import pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FOLDER = os.path.dirname(os.path.dirname(HERE))
OUT_DIR = os.path.join(FOLDER, 'output', 'frames')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pf = _load('_progressive_frames', os.path.join(HERE, 'progressive_frames.py'))
plot = _load('_plot_k13_inhib_isobutanol',
             os.path.join(FOLDER, 'code', 'plots',
                          'plot_k13_inhib_isobutanol_opt_PI_TRY_informed.py'))
ksf = plot.ksf

# Frame sequence: (frame name, [groups the frame adds]). Groups (see the
# module docstring): 'IRR map', 'break-even contour', 'colour bar',
# 'colour bar break-even line', 'profitability', 'ethanol titer',
# 'ethanol yield', 'ethanol productivity', 'isobutanol titer',
# 'isobutanol yield', 'isobutanol productivity', 'cell density', 'TCI'.
FRAMES = [
    ('axes', []),
    ('irr-map', ['IRR map', 'colour bar']),
    ('break-even', ['break-even contour', 'colour bar break-even line']),
    ('profitability', ['profitability']),
    ('ethanol-optima', ['ethanol titer', 'ethanol yield', 'ethanol productivity']),
    ('isobutanol-optima', ['isobutanol titer', 'isobutanol yield',
                           'isobutanol productivity']),
    ('cell-density-tci', ['cell density', 'TCI']),
]

# = ksf.main's savefig options
DPI = 600
FACECOLOR = 'white'


def _group_name(label):
    return label.replace('\n', ' ')


@contextlib.contextmanager
def _recorded_draw_panel(calls):
    """Run ksf.draw_panel unchanged, appending (kind, returned artist) to
    `calls` for each ax.contourf / ax.contour / ax.plot / ax.annotate and
    cax.axhline call it makes."""
    draw_panel = ksf.draw_panel

    def instrumented(figure, fig, ax, cax):
        def record(obj, method, kind):
            original = getattr(obj, method)

            def wrapper(*args, **kwargs):
                artist = original(*args, **kwargs)
                calls.append((kind, artist))
                return artist
            setattr(obj, method, wrapper)

        for method in ('contourf', 'contour', 'plot', 'annotate'):
            record(ax, method, method)
        record(cax, 'axhline', 'cbar_line')
        try:
            return draw_panel(figure, fig, ax, cax)
        finally:
            for method in ('contourf', 'contour', 'plot', 'annotate'):
                delattr(ax, method)
            delattr(cax, 'axhline')

    ksf.draw_panel = instrumented
    try:
        yield
    finally:
        ksf.draw_panel = draw_panel


def _only(calls, kind):
    found = [artist for k, artist in calls if k == kind]
    if len(found) != 1:
        raise RuntimeError(f'expected one {kind} call while drawing the panel, '
                           f'got {len(found)}')
    return found[0]


def build_figure():
    """Draw the complete figure; returns (fig, {group name: [artists]})."""
    calls = []
    with _recorded_draw_panel(calls):
        fig, _ = ksf.make_figure(plot.FIGURE)
    _, cax = fig.axes
    groups = {
        'IRR map': [_only(calls, 'contourf')],
        'break-even contour': [_only(calls, 'contour')],
        'colour bar': [cax],
        'colour bar break-even line': [_only(calls, 'cbar_line')],
    }
    # each optimum: ax.plot draws its marker, then ax.annotate its callout
    markers = [a for k, a in calls if k == 'plot']
    callouts = [a for k, a in calls if k == 'annotate']
    expected = [_group_name(row[2]) for row in plot.FIGURE.optima]
    if len(markers) != len(callouts) or len(callouts) != len(expected):
        raise RuntimeError(f'{len(markers)} markers and {len(callouts)} callouts '
                           f'drawn for {len(expected)} optima')
    for (marker,), callout in zip(markers, callouts):
        groups[_group_name(callout.get_text())] = [marker, callout]
    missing = sorted(set(expected) - set(groups))
    if missing:
        raise RuntimeError(f'no marker + callout drawn for {missing}')
    return fig, groups


def main(dpi=DPI, formats=('png',)):
    fig, groups = build_figure()
    bbox = pf.tight_bbox(fig, dpi)
    paths = pf.save_frames(fig, groups, FRAMES, OUT_DIR, plot.FIGURE.output_stem,
                           dpi, formats, bbox_inches=bbox, facecolor=FACECOLOR)
    plt.close(fig)
    return paths


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--dpi', type=int, default=DPI)
    parser.add_argument('--formats', nargs='+', default=['png'])
    args = parser.parse_args()
    main(args.dpi, args.formats)
