# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Build GIFs of Extended Data Figure 8 (IRR over ALS capacity k_13 x
isobutanol inhibition, opt_PI_TRY_informed), one GIF per step towards the
full figure, for presentations.

The figure is drawn by the unmodified snapshot
plots/plot_k13_inhib_isobutanol_opt_PI_TRY_informed.py (through
plots/_kinetic_sweep_figure.py); this script only sorts its artists into
groups and reveals them step by step with gif_engine.Storyboard:

    start   the empty axes: frame, ticks and axis titles
    1       the colour bar, then the IRR map swept in along k_13
    2       the break-even (IRR = 0) contour and its line on the colour bar
    3       the ethanol titer, yield and productivity optima
    4       the isobutanol titer, yield and productivity optima
    5       the cell-density and TCI optima
    6       the profitability (IRR) optimum

Each landmark (marker, label and arrow) fades in as one, the landmarks of a
step one after another. The last frame is the complete figure (checked
pixel for pixel against this script's own render of it). Run through
../make_gifs.py, which points
the results root at ../data; see gif_engine.py for the files written and the
GIF options.

Sim-safe: reads the saved sweep CSVs only.
"""
import argparse
import importlib.util
import os
import sys

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import gif_engine as ge  # noqa: E402
from matplotlib.colors import same_color  # noqa: E402
from matplotlib.contour import ContourSet  # noqa: E402
from matplotlib.text import Annotation  # noqa: E402

STEM = 'k13_inhib_isobutanol_opt_PI_TRY_informed_IRR'
#: how the plot script saves the figure (besides the file name and dpi)
SAVEFIG_KWARGS = dict(bbox_inches='tight', facecolor='white')
#: seconds: the colour bar fading in, the map sweep, the contour fading in,
#: one landmark fading in
CBAR_S, MAP_S, CONTOUR_S, LANDMARK_S = 0.8, 3.0, 1.2, 0.7

# landmark steps: the label of each optimum, as the plot script writes it,
# in the order they appear
ETHANOL = ('ethanol titer', 'ethanol yield', 'ethanol\nproductivity')
ISOBUTANOL = ('isobutanol\ntiter', 'isobutanol\nyield',
              'isobutanol\nproductivity')
CELLS_AND_TCI = ('cell density', 'TCI')
PROFITABILITY = ('profitability',)


def load_plot_module():
    path = os.path.join(HERE, 'plots',
                        'plot_k13_inhib_isobutanol_opt_PI_TRY_informed.py')
    spec = importlib.util.spec_from_file_location('_k13_inhib_plot', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def landmarks(ax, optima):
    """{label: (marker, annotation)} for every optimum of the figure. The
    annotation draws the label and its arrow; the marker is the Line2D with
    the row's marker shape and face colour (unique per row)."""
    annotations = {a.get_text(): a for a in ax.texts
                   if isinstance(a, Annotation)}
    found = {}
    for _, _, label, marker, face, *_ in optima:
        (line,) = [line for line in ax.lines
                   if line.get_marker() == marker
                   and same_color(line.get_markerfacecolor(), face)]
        found[label] = (line, annotations[label])
    return found


def build_storyboard(fig, figure):
    ax, cax = fig.axes                      # the panel, then the colour bar
    filled, zero_contour = [c for c in ax.collections
                            if isinstance(c, ContourSet)]  # drawing order
    (cbar_zero_line,) = cax.lines           # the axhline at IRR = 0

    board = ge.Storyboard(fig, savefig_kwargs=SAVEFIG_KWARGS)
    # the colour-bar axes as a whole (colours, ticks, title); its IRR = 0
    # line is grouped with the contour and stays hidden until then
    board.group('colour bar', cax)
    board.group('IRR map', filled)
    board.group('break-even', zero_contour, cbar_zero_line)
    marks = landmarks(ax, figure.optima)
    staged = ETHANOL + ISOBUTANOL + CELLS_AND_TCI + PROFITABILITY
    if set(marks) != set(staged):
        raise ValueError(
            'the landmarks of the plot script and of this storyboard differ: '
            f'{sorted(set(marks) ^ set(staged))}')
    for label in staged:
        board.group(label, marks[label])
    # anything else the plot script draws on the panel (a new marker, a
    # callout) must be given a step rather than show from the first frame
    board.require_grouped(ax.lines, ax.texts, ax.collections, ax.patches,
                          ax.images, fig.texts, fig.lines, fig.patches,
                          fig.images, fig.legends)

    def one_by_one(labels):
        return [ge.fade(label, seconds=LANDMARK_S) for label in labels]

    board.step('irr-map',
               'IRR is highest at moderate ALS capacity and low isobutanol '
               'inhibition, and drops off at high ALS capacity',
               ge.fade('colour bar', seconds=CBAR_S),
               ge.wipe('IRR map', region=ge.Region(ax, pad=0),
                       seconds=MAP_S))
    board.step('break-even',
               'The white contour marks break-even (IRR = 0): the process is '
               'profitable mainly from just above zero up to about 9-13 g/g/h '
               'of ALS capacity, depending on inhibition, plus small marginal '
               'pockets (IRR below 4 %) at 14-16 g/g/h and low inhibition',
               ge.fade('break-even', seconds=CONTOUR_S))
    board.step('ethanol-optima',
               'Ethanol titer, yield and productivity peak at almost no ALS '
               'capacity, away from the most profitable region',
               *one_by_one(ETHANOL))
    board.step('isobutanol-optima',
               'Isobutanol titer, yield and productivity peak at about '
               '5-6 g/g/h of ALS capacity, at or near the lowest inhibition',
               *one_by_one(ISOBUTANOL))
    board.step('cell-density-and-tci',
               'Cell density peaks with no ALS capacity, and capital cost '
               '(TCI) is lowest at the highest ALS capacity, where the '
               'process loses money',
               *one_by_one(CELLS_AND_TCI))
    board.step('profitability-optimum',
               'Profitability peaks at about 4 g/g/h of ALS capacity and the '
               'lowest inhibition: near, but not at, the isobutanol optima',
               *one_by_one(PROFITABILITY))
    return board


def main(argv=None):
    parser = argparse.ArgumentParser(prog='make_gifs.py',
                                     description=__doc__.split('\n\n')[0])
    ge.add_cli_arguments(parser)
    args = parser.parse_args(argv)
    out_dir = args.out_dir or os.path.join(os.path.dirname(HERE), 'output',
                                           'gifs')
    font = ge.use_arial_or_metric_twin()
    if font != 'Arial':
        print(f'Arial is not installed; using {font or "the fallback font"}')
    plot = load_plot_module()
    fig = ge.capture_figure(plot.ksf.main, plot.FIGURE)
    build_storyboard(fig, plot.FIGURE).render(out_dir, STEM,
                                              **ge.render_kwargs(args))


if __name__ == '__main__':
    main()
