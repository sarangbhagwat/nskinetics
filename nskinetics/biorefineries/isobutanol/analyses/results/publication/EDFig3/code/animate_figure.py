# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Build GIFs of Extended Data Figure 3 (panel A: the six process-level TRY
trajectories; panel B: final proteome allocation), one GIF per step towards
the full figure, for presentations.

The figure is drawn by the unmodified snapshot
nskinetics/biorefineries/isobutanol/plots/plot_kin_opt_parameter_sets.py
(--variants process), recorded and grouped with kin_opt_animation.py and
animated with gif_engine.Storyboard. Only the empty band below panel B is
trimmed; everything else is drawn as make_figure.py draws it.

Storyboard:

    start   the panel-A title and the six empty cells (isobutanol yield,
            titer, productivity; ethanol yield, titer, productivity)
    1  baseline                 the baseline reference line in every cell
    2  isobutanol-campaigns     each isobutanol cell's own campaign: its
                                trials and solid incumbent, drawn along the
                                trial axis
    3  ethanol-campaigns        the same for the ethanol cells
    4  isobutanol-crosscheck    the isobutanol campaigns' incumbents in the
                                five cells each did not optimize, dashed
    5  ethanol-crosscheck       the same for the ethanol campaigns
    6  profitability-uninformed the uninformed profitability campaign's
                                incumbent in every cell, dashed
    7  profitability-try-informed  the TRY-informed campaign's incumbent from
                                trial 1000, dashed
    8-12                        panel B, as in kin_opt_animation.py

Legend entries appear with what they describe. Run through ../make_gifs.py;
see gif_engine.py for the files written and the GIF options.

Sim-safe: reads the campaign CSVs, two workbooks and the antimony text only.
"""
import os
import sys

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import gif_engine as ge  # noqa: E402
import kin_opt_animation as ka  # noqa: E402

STEM = 'kinBO_parameter_sets_process'


def build_storyboard(fig, recorder):
    cells = recorder.cells
    isobutanol = [c for c in cells if c.col.startswith('IBO')]
    ethanol = [c for c in cells if c.col.startswith('EtOH')]
    if len(isobutanol) != 3 or len(ethanol) != 3:
        raise RuntimeError('expected three isobutanol and three ethanol cells')

    baseline = recorder.labels('baseline')[0]
    uninformed = recorder.labels('uninformed')[0]
    try_informed = recorder.labels('try-informed')[0]
    entry = lambda *labels: ka.legend_entries(fig, *labels)  # noqa: E731

    def in_cells(cells, kind, labels=None):
        return [c.find(kind, labels) for c in cells]

    board = ge.Storyboard(fig, ka.bottom_crop(fig, recorder))
    board.group('baseline', in_cells(cells, 'baseline'),
                ka.campaign_key_frame(fig, recorder), entry(baseline))
    board.group('isobutanol own', in_cells(isobutanol, 'cloud'),
                in_cells(isobutanol, 'incumbent'),
                entry(*recorder.labels('isobutanol'), 'Individual trials',
                      'Incumbent (campaign optimizing the shown metric)'))
    board.group('ethanol own', in_cells(ethanol, 'cloud'),
                in_cells(ethanol, 'incumbent'),
                entry(*recorder.labels('ethanol')))
    board.group('isobutanol other',
                in_cells(cells, 'incumbent-other',
                         recorder.labels('isobutanol')),
                entry('Incumbent (other campaign)'))
    board.group('ethanol other',
                in_cells(cells, 'incumbent-other', recorder.labels('ethanol')))
    board.group('uninformed', in_cells(cells, 'incumbent-other', {uninformed}),
                entry(uninformed))
    board.group('try-informed',
                in_cells(cells, 'incumbent-other', {try_informed}),
                entry(try_informed))
    ka.group_panel_b(board, fig, recorder)
    ka.require_data_grouped(board, recorder)

    n_preloaded = recorder.set_of(try_informed)['n_preloaded']
    whole = lambda cells: [ge.Region(c.ax) for c in cells]  # noqa: E731
    after_preload = [ge.Region(c.ax, x=(n_preloaded, c.ax.get_xlim()[1]))
                     for c in cells]
    T = ka.TRACE_S

    board.step('baseline',
               'The dotted grey line is the baseline strain before any '
               'optimization; it makes no isobutanol, so the line lies on '
               'the axis in the isobutanol cells',
               ge.fade('baseline', seconds=1.0))
    board.step('isobutanol-campaigns',
               'Each isobutanol cell shows the campaign that optimized that '
               'metric: its trials and its incumbent',
               ge.wipe('isobutanol own', region=whole(isobutanol), seconds=T))
    board.step('ethanol-campaigns',
               'Likewise for the three ethanol metrics',
               ge.wipe('ethanol own', region=whole(ethanol), seconds=T))
    board.step('isobutanol-crosscheck',
               'Dashed: the isobutanol campaigns on the metrics they did not '
               'optimize; they end up making almost no ethanol',
               ge.wipe('isobutanol other', region=whole(cells), seconds=T))
    board.step('ethanol-crosscheck',
               'Dashed: the ethanol campaigns on the metrics they did not '
               'optimize; they make almost no isobutanol',
               ge.wipe('ethanol other', region=whole(cells), seconds=T))
    board.step('profitability-uninformed',
               'Like the ethanol campaigns, the uninformed profitability '
               'campaign makes ethanol but almost no isobutanol',
               ge.wipe('uninformed', region=whole(cells), seconds=T))
    board.step('profitability-try-informed',
               'The TRY-informed profitability campaign, from trial 1000, '
               'co-produces isobutanol and ethanol',
               ge.wipe('try-informed', region=after_preload, seconds=2.0))
    ka.panel_b_steps(board, recorder)
    return board


def main(argv=None):
    ka.run(__doc__.split('\n\n')[0], 'process', STEM, build_storyboard, argv)


if __name__ == '__main__':
    main()
