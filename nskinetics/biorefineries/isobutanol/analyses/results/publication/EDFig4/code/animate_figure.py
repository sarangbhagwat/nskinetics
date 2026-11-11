# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Build GIFs of the Extended Data Figure 4 plot panel (ethanol yield, PI and
IRR along the 1-D k_3 slice through the TRY-informed optimum), one GIF per
step towards the full figure, for presentations.

The figure is drawn by the unmodified snapshot plots/plot_axis_try_informed_sweep.py
(plot_param); this script only sorts its artists into groups and reveals them
step by step with gif_engine.Storyboard:

    start   the frame and the Pdc-capacity axis
    1       the ethanol-yield axis, then the yield curve traced left to right
    2       the profitability-index axis, then the PI curve
    3       the IRR axis, then the IRR curve
    4       the burden-derating marker and its legend entry

The last frame is the complete figure (checked pixel for pixel against this
script's own render of it). Run through ../make_gifs.py, which points the
results root at ../data;
see gif_engine.py for the files written and the GIF options.

Sim-safe: reads the saved sweep CSV + anchor JSON only.
"""
import argparse
import importlib.util
import os
import sys

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import gif_engine as ge  # noqa: E402

STEM = 'k_3_try_informed_sweep_etoh_yield'
PARAM = 'k_3'
#: seconds: an axis fading in, a curve being traced
AXIS_S, CURVE_S = 0.6, 2.4


def load_plot_module():
    path = os.path.join(HERE, 'plots', 'plot_axis_try_informed_sweep.py')
    spec = importlib.util.spec_from_file_location('_axis_sweep_plot', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def y_axis_artists(ax, side):
    """The y axis of one of the three overlaid axes: its coloured spine, the
    ticks, the tick labels and the axis title."""
    return [ax.spines[side], ax.yaxis]


def build_storyboard(fig):
    ax, ax_yield, ax_irr = fig.axes      # creation order in plot_param
    (derate_line, pi_line), (yield_line,), (irr_line,) = (
        ax.lines, ax_yield.lines, ax_irr.lines)
    (legend,) = fig.legends

    board = ge.Storyboard(fig)
    board.group('yield axis', y_axis_artists(ax_yield, 'right'))
    board.group('yield curve', yield_line)
    board.group('PI axis', y_axis_artists(ax, 'left'))
    board.group('PI curve', pi_line)
    board.group('IRR axis', y_axis_artists(ax_irr, 'right'))
    board.group('IRR curve', irr_line)
    board.group('derating', derate_line, legend)
    for axes in fig.axes:   # nothing the plot script draws is left ungrouped
        board.require_grouped(axes.lines, axes.collections, axes.texts)

    panel = ge.Region(ax)
    board.step('ethanol-yield',
               'Ethanol yield rises smoothly with Pdc capacity',
               ge.fade('yield axis', seconds=AXIS_S),
               ge.wipe('yield curve', region=panel, seconds=CURVE_S))
    board.step('profitability-index',
               'The profitability index over the same slice is rugged, with '
               'cliffs the yield does not show',
               ge.fade('PI axis', seconds=AXIS_S),
               ge.wipe('PI curve', region=panel, seconds=CURVE_S))
    board.step('irr',
               'IRR follows the profitability index',
               ge.fade('IRR axis', seconds=AXIS_S),
               ge.wipe('IRR curve', region=panel, seconds=CURVE_S))
    board.step('burden-derating',
               'Above this Pdc capacity the enzyme burden derates growth',
               ge.wipe('derating', region=ge.Region(ax, pad=0),
                       direction='up', seconds=0.9))
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
    fig = ge.capture_figure(plot.plot_param, PARAM, 300)
    build_storyboard(fig).render(out_dir, STEM, **ge.render_kwargs(args))


if __name__ == '__main__':
    main()
