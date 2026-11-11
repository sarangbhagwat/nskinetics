# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Build GIFs of Figure 4 (Monte Carlo MESP against TCI, ethanol production,
DDGS revenue, ethanol yield, titer and productivity; panels A-F), one GIF per
step towards the full figure, for presentations.

The figure is drawn by the unmodified snapshot
plots/plot_uncertainty_MPSP_vs_TCI.py; this script only sorts its artists into
groups and reveals them step by step with gif_engine.Storyboard:

    start   the six empty panels: axes, axis titles and panel letters
    1       panel A: the baseline diamond and its callout
    2       panel A: the density of the Monte Carlo samples and Spearman's rho
    3       panel A: the marginal box plots
    4       panel A: the Pareto frontier, traced along TCI, and its callout
    5       panel A: the conventional corn ethanol box and its callout
    6       panel A: the gasoline price band and its callout
    7-11    panels B to F, one per step: the samples (density, baseline, rho
            and box plots), then the binned-median line traced along x (named
            by its callout in panel B)
    12      panels B to F: the two grey ranges, panel after panel

Panel A is built piece by piece because its callouts name every element of the
figure; panels B to F then repeat those elements, in the same order. The last
frame is the complete figure (checked pixel for pixel against this script's
own render of it). Run through ../make_gifs.py, which points the results root
at ../data; see gif_engine.py for the files written and the GIF options.

Sim-safe: reads the two saved uncertainty workbooks only.
"""
import argparse
import importlib.util
import os
import sys

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import gif_engine as ge  # noqa: E402
from matplotlib.contour import ContourSet  # noqa: E402
from matplotlib.text import Annotation  # noqa: E402

#: the 6000-simulation scenario-A run make_figure.py plots
TAG = "_IBO_2026.9.29-23.45_['A']_6000sims_A"
STEM = 'uncertainty_bivariate_A_2026.9.29-23.45_6000sims'
#: how the plot script saves the figure (besides the file name and dpi)
SAVEFIG_KWARGS = dict(facecolor='white')
#: seconds: a small element fading in, the density, a line being traced
FADE_S, DENSITY_S, TRACE_S = 0.7, 1.2, 1.4

# panels B to F: slug and caption of each panel's step. The numbers are those
# the plot script reports for this run: Spearman's rho, and the change of the
# binned-median MESP from the first to the last tenth of the samples along x.
LATER_PANELS = {
    'B': ('ethanol-production',
          'Panel B: MESP falls with ethanol production (rho = -0.43); the '
          'line is the binned median, which drops 0.21 $/GGE from the lowest '
          'to the highest tenth of the samples'),
    'C': ('ddgs-revenue',
          'Panel C: MESP falls most clearly with DDGS revenue (rho = -0.51; '
          'the binned median drops 0.23 $/GGE)'),
    'D': ('ethanol-yield',
          'Panel D: MESP depends only weakly on ethanol yield (rho = -0.14; '
          'the binned median drops 0.07 $/GGE)'),
    'E': ('ethanol-titer',
          'Panel E: MESP depends only weakly on ethanol titer (rho = -0.19; '
          'the binned median drops 0.09 $/GGE)'),
    'F': ('ethanol-productivity',
          'Panel F: MESP depends only weakly on ethanol productivity '
          '(rho = -0.16; the binned median drops 0.07 $/GGE)'),
}
# the last step; the shares are the samples inside each light box
RANGES_CAPTION = (
    'Panels B-F: the light boxes are the ranges of conventional corn ethanol '
    'biorefineries (B, C) and of high-gravity fermentation (D-F), over the '
    'same gasoline price band; they hold 93 % and 83 % of the samples in B '
    'and C, all but one in D, over 99 % in E, and 56 % in F, where the rest '
    'and the baseline lie below 2.0 g/L/h')


def load_plot_module():
    path = os.path.join(HERE, 'plots', 'plot_uncertainty_MPSP_vs_TCI.py')
    spec = importlib.util.spec_from_file_location('_uncertainty_plot', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def one(items, what):
    """The only item of a list; raises with `what` otherwise."""
    items = list(items)
    if len(items) != 1:
        raise ValueError(f'expected one {what}, found {len(items)}')
    return items[0]


class Panel:
    """The artists of one joint panel and its two marginal axes, by what they
    are: `box` is the light box of conventional values, `band` the gasoline
    price band, `labels` any plain text besides the panel letter and rho
    (panel A's gasoline label)."""

    def __init__(self, letter, ax, ax_top, ax_right):
        self.letter, self.ax = letter, ax
        self.marginals = [ax_top, ax_right]
        self.density = one((c for c in ax.collections
                            if isinstance(c, ContourSet)), 'density')
        self.box = one((c for c in ax.collections
                        if not isinstance(c, ContourSet)), 'light box')
        self.band = one(ax.patches, 'gasoline band')
        self.baseline = one((line for line in ax.lines
                             if line.get_marker() == 'D'), 'baseline diamond')
        # the Pareto frontier (panel A) or the binned-median line (B to F)
        self.curve = one((line for line in ax.lines
                          if line is not self.baseline), 'curve')
        annotations = [a for a in ax.texts if isinstance(a, Annotation)]
        self.callouts = {a.get_text(): a for a in annotations}
        if len(self.callouts) != len(annotations):
            raise ValueError(f'panel {letter}: two callouts share a text')
        plain = [t for t in ax.texts if not isinstance(t, Annotation)]
        self.letter_text = one((t for t in plain if t.get_text() == letter),
                               'panel letter')
        self.rho = one((t for t in plain if r'\rho' in t.get_text()), 'rho')
        self.labels = [t for t in plain
                       if t is not self.letter_text and t is not self.rho]

    def callout(self, text):
        if text not in self.callouts:
            raise KeyError(f'panel {self.letter} has no callout {text!r}')
        return self.callouts.pop(text)


def build_storyboard(fig, plot):
    axes = fig.axes                     # joint, top box, right box per panel
    if len(axes) != 3 * len(plot.PANELS):
        raise ValueError(f'expected three axes per panel, found {len(axes)}')
    panels = [Panel(chr(ord('A') + i), *axes[3 * i:3 * i + 3])
              for i in range(len(plot.PANELS))]
    first, later = panels[0], panels[1:]
    if [p.letter for p in later] != list(LATER_PANELS):
        raise ValueError('the panels of the plot script and of this '
                         'storyboard differ')

    board = ge.Storyboard(fig, savefig_kwargs=SAVEFIG_KWARGS)
    a = first
    board.group('A baseline', a.baseline,
                a.callout(plot.BASELINE_CALLOUT['text']))
    board.group('A density', a.density, a.rho)
    board.group('A marginals', a.marginals)
    board.group('A frontier', a.curve)
    board.group('A frontier callout', a.callout(plot.PARETO_CALLOUT['text']))
    board.group('A conventional', a.box,
                a.callout(plot.CONVENTIONAL_CALLOUT['text']))
    # the double arrow across the band is an annotation without text; its
    # label is a plain text on top of it
    board.group('A gasoline', a.band, a.callout(''),
                one(a.labels, 'gasoline label in panel A'))
    for p in later:
        if p.labels:
            raise ValueError(f'panel {p.letter}: texts in no step: '
                             f'{[t.get_text() for t in p.labels]}')
        board.group(f'{p.letter} samples', p.density, p.baseline, p.rho,
                    p.marginals)
        board.group(f'{p.letter} trend', p.curve)
        board.group(f'{p.letter} ranges', p.box, p.band)
    trend_callout = later[0].callout(plot.TREND_CALLOUT['text'])
    board.group('B trend callout', trend_callout)

    # anything else the plot script draws (a new callout, another line) must
    # be given a step rather than show from the first frame; the panel
    # letters belong to the empty panels the build starts from
    for p in panels:
        if p.callouts:
            raise ValueError(f'panel {p.letter}: callouts in no step: '
                             f'{sorted(p.callouts)}')
        board.require_grouped(p.ax.lines, p.ax.collections, p.ax.patches,
                              p.ax.images,
                              [t for t in p.ax.texts
                               if t is not p.letter_text])
    board.require_grouped(fig.lines, fig.patches, fig.images, fig.legends)

    board.step('baseline',
               'Panel A: the baseline design has a minimum ethanol selling '
               'price (MESP) of 3.83 $/GGE at a total capital investment of '
               '139 MM$',
               ge.fade('A baseline', seconds=FADE_S))
    board.step('uncertainty',
               'Across 6000 Monte Carlo samples (darker = denser) MESP '
               'ranges from 3.40 to 4.36 $/GGE and is weakly anticorrelated '
               'with capital investment (Spearman rho = -0.26)',
               ge.fade('A density', seconds=DENSITY_S))
    board.step('box-plots',
               'The box plots summarize each axis (box 25th-75th percentile, '
               'whiskers 5th-95th, dots minimum and maximum): the middle '
               '90 % of MESP is 3.64-4.08 $/GGE, and of capital investment '
               '125-151 MM$',
               ge.fade('A marginals', seconds=FADE_S))
    board.step('pareto-frontier',
               'The Pareto frontier joins the samples that no other sample '
               'beats on both capital investment and MESP',
               ge.wipe('A frontier', seconds=TRACE_S,
                       region=ge.Region(artists=[a.curve])),
               ge.fade('A frontier callout', seconds=FADE_S))
    board.step('conventional-ranges',
               'The light box spans the ethanol market price range '
               '(2.31-5.15 $/GGE) and the capital investment of conventional '
               'corn ethanol biorefineries (110-196 MM$): every sample lies '
               'inside it',
               ge.fade('A conventional', seconds=FADE_S))
    board.step('gasoline-range',
               'The dark band is the gasoline market price range '
               '(2.16-4.84 $/GGE): every sample lies inside it too',
               ge.fade('A gasoline', seconds=FADE_S))
    for p in later:
        slug, caption = LATER_PANELS[p.letter]
        phases = [ge.fade(f'{p.letter} samples', seconds=DENSITY_S),
                  ge.wipe(f'{p.letter} trend', seconds=TRACE_S,
                          region=ge.Region(artists=[p.curve]))]
        if p is later[0]:
            phases.append(ge.fade('B trend callout', seconds=FADE_S))
        board.step(slug, caption, *phases)
    board.step('conventional-ranges-b-to-f', RANGES_CAPTION,
               *[ge.fade(f'{p.letter} ranges', seconds=FADE_S)
                 for p in later])
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
    results = os.path.join(plot.RESULTS_DIR, TAG + '_1_full_evaluation.xlsx')
    baseline = os.path.join(plot.RESULTS_DIR, TAG + '_0_baseline.xlsx')
    os.makedirs(out_dir, exist_ok=True)   # the plot function expects it
    fig = ge.capture_figure(plot.plot_uncertainty_MPSP_vs_TCI, results,
                            baseline, out_dir=out_dir)
    build_storyboard(fig, plot).render(out_dir, STEM,
                                       **ge.render_kwargs(args))


if __name__ == '__main__':
    main()
