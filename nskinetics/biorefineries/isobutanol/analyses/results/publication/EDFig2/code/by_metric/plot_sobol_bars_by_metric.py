# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Sobol' / Shapley bars of Extended Data Figure 2, one metric per figure,
with the decision variables grouped on the x axis.

Extended Data Figure 2 stacks five metrics over one x axis sorted by the
Shapley effect on the profitability index. This script draws ONE metric per
figure and arranges the twelve decision variables by what they are (GROUPS):
the six enzyme expression levels in pathway order (glycolysis, ethanol
production, isobutanol production), the three inhibition levels and the three
substrate feeding variables, each group under its own header. The order is the
same for every metric and every figure is saved on the same canvas, so the
figures of different metrics overlay exactly on consecutive slides.

The bars themselves (first-order Sobol' index, Shapley effect, total-order
Sobol' index, each with its replicate sd) are drawn by the folder's unmodified
analyses/sensitivity/analyze_sobol_split12d.py, loaded by file path, from the
same saved <study>_sobol_indices.csv: same values, colours and tick style as
the manuscript figure. Any metric of that table can be drawn (--list); a
metric whose surrogate is unreliable (CV Q2 < 0.8) is marked as such.

Usage (through ../../make_figure_by_metric.py, which points the results root
NSKINETICS_IBO_RESULTS_DIR at the folder's ./data):
    python plot_sobol_bars_by_metric.py --study-name STUDY
           [--metric NAME [NAME ...]] [--list] [--dpi N]
           [--formats png [pdf svg ...]]
Without --metric: the five metrics of the manuscript figure. Figures go to
../../output/by_metric/EDFig2_bars_<metric>.<ext>.

Sim-safe: reads the saved index table only; no surrogate fit, no simulation.
"""
import argparse
import importlib.util
import os
import re
import sys

sys.dont_write_bytecode = True   # no __pycache__ under code/ when run directly

HERE = os.path.dirname(os.path.abspath(__file__))
FOLDER = os.path.dirname(os.path.dirname(HERE))
OUT_DIR = os.path.join(FOLDER, 'output', 'by_metric')

# The x axis: (group header, [(sub-header, [decision variables])]), left to
# right. The enzyme expression levels follow the pathway: glycolysis, the
# ethanol branch (Pdc, Adh1), the isobutanol branch (ALS, Ilv5 / Ilv3 /
# Aro10, Adh6). A sub-header of None draws none.
GROUPS = [
    ('Enzyme expression levels', [
        ('glycolysis', ['glycolysis']),
        ('ethanol production', ['k_3', 'k_6']),
        ('isobutanol production', ['k_13', 'ehrlich_downstream', 'k_17'])]),
    ('Inhibition levels', [
        (None, ['inhib_ethanol', 'inhib_isobutanol', 'inhib_acetate'])]),
    ('Substrate feeding', [
        (None, ['threshold_conc', 'target_delta', 'max_n_spikes'])]),
]

# Metric names for the y-axis title: the plot script's STACK_TITLES (the five
# metrics of the manuscript figure) plus the table's other metrics
MORE_TITLES = {'IBO productivity': 'Isobutanol productivity',
               'tau': 'Fermentation time',
               'n_glu_spikes': 'Number of glucose spikes',
               'TCI': 'Total capital investment'}

# Canvas [in]: the same for every metric (no tight crop), so figures overlay
FIGSIZE = (7.8, 5.2)
MARGINS = dict(left=0.165, right=0.985, bottom=0.32, top=0.835)
LEGEND_Y = 0.995                 # top of the legend row [figure fraction]
# value axis: the manuscript panels' (-0.04, 1.0) with a little more room, so
# that no error bar is cut off for any metric of the table (a share of ~1 with
# its sd reaches past 1; a first-order index near 0 can dip below -0.04)
YLIM = (-0.06, 1.06)
DPI = 600                        # as the manuscript figure
# headers above the axes [axes fraction]: bracket line, group header,
# sub-header; the bracket stops this far [x units] short of a slot's edge
BRACKET_Y, HEADER_Y, SUBHEADER_Y, BRACKET_INSET = 1.075, 1.095, 1.02, 0.08
DIVIDER = dict(color='0.55', linewidth=0.8)                    # between groups
SUBDIVIDER = dict(color='0.75', linewidth=0.6, linestyle=(0, (3, 3)))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


sobol = _load('_analyze_sobol_split12d', os.path.join(
    FOLDER, 'code', 'analyses', 'sensitivity', 'analyze_sobol_split12d.py'))
plt = sobol.plt
TITLES = {**sobol.STACK_TITLES, **MORE_TITLES}


def load_table(study_name):
    """The study's saved index table (as the plot script's --replot reads it)."""
    return sobol.pd.read_csv(sobol.paths.require_input(
        os.path.join(sobol.RESULTS, study_name + '_sobol_indices.csv'),
        produced_by='a full run of analyses/sensitivity/'
                    'analyze_sobol_split12d.py --study-name <study>'))


def grouped_order(table):
    """The decision variables in GROUPS order; raises unless GROUPS holds
    exactly the table's parameters, each once."""
    order = [p for _, subgroups in GROUPS for _, ps in subgroups for p in ps]
    if sorted(order) != sorted(set(table.parameter)):
        raise ValueError(f'GROUPS {sorted(order)} does not match the index '
                         f"table's parameters {sorted(set(table.parameter))}")
    return order


def resolve_metric(name, table):
    """The table's metric called `name`: its name in the table or its title,
    whatever the case ('EtOH yield' or 'ethanol yield')."""
    known = list(dict.fromkeys(table.metric))
    for m in known:
        if name.lower() in (m.lower(), TITLES.get(m, m).lower()):
            return m
    raise SystemExit(f'no metric {name!r} in the index table (names or '
                     'titles; --list shows both); it has: '
                     + ', '.join(repr(m) for m in known))


def draw_headers(ax):
    """Group headers with brackets and sub-headers above the axes, a
    divider between neighbouring groups and a fainter one between the
    sub-groups of a group (x in data units, one slot per variable; y in axes
    fractions)."""
    at = ax.get_xaxis_transform()
    x = 0
    for i, (header, subgroups) in enumerate(GROUPS):
        x0 = x
        for j, (subheader, parameters) in enumerate(subgroups):
            if subheader:
                ax.text(x + (len(parameters) - 1) / 2, SUBHEADER_Y, subheader,
                        transform=at, ha='center', va='bottom', fontsize=8,
                        color='0.25', clip_on=False)
            x += len(parameters)
            if j < len(subgroups) - 1:
                ax.axvline(x - 0.5, zorder=0.4, **SUBDIVIDER)
        left, right = x0 - 0.5 + BRACKET_INSET, x - 0.5 - BRACKET_INSET
        ax.plot([left, right], [BRACKET_Y] * 2, transform=at, color='0.25',
                linewidth=0.9, solid_capstyle='butt', clip_on=False)
        ax.text((left + right) / 2, HEADER_Y, header, transform=at,
                ha='center', va='bottom', fontsize=10, fontweight='bold',
                clip_on=False)
        if i < len(GROUPS) - 1:
            ax.axvline(x - 0.5, zorder=0.4, **DIVIDER)


def plot_metric(table, metric):
    """The figure of one metric: its three bars per decision variable, in
    GROUPS order under the group headers."""
    sobol._style()
    order = grouped_order(table)
    t = table[table.metric == metric]
    fig, ax = plt.subplots(figsize=FIGSIZE)
    fig.subplots_adjust(**MARGINS)
    sobol._draw_index_bars(ax, t, order)
    ax.set_xlim(-0.5, len(order) - 0.5)
    ax.set_ylim(*YLIM)
    ax.set_ylabel('Share of variance for\n' + TITLES.get(metric, metric),
                  fontsize=12, fontweight='bold')
    sobol._ticks(ax)
    draw_headers(ax)
    # the key in one row above the group headers, centred on the axes, so the
    # panel is free of it for every metric (inside, it would straddle a group
    # divider). The surrogate fit goes top right, over the feeding variables:
    # the manuscript panels' top-left corner now holds the glycolysis bars.
    fig.legend(*ax.get_legend_handles_labels(), loc='upper center', ncol=3,
               bbox_to_anchor=((MARGINS['left'] + MARGINS['right']) / 2,
                               LEGEND_Y),
               frameon=False, columnspacing=1.8)
    reliable = bool(t.reliable.iloc[0])
    ax.text(0.985, 0.95, f'Q$^2$ = {t.Q2.iloc[0]:.2f}'
            + ('' if reliable else ' (unreliable)'),
            transform=ax.transAxes, ha='right', va='top', fontsize=10)
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--study-name', required=True,
                    help='the Sobol study whose <study>_sobol_indices.csv to '
                         'plot (set by make_figure_by_metric.py)')
    ap.add_argument('--metric', nargs='+', default=None, metavar='NAME',
                    help='metric(s) to plot, one figure each, by table name '
                         "or title, e.g. 'EtOH yield' or 'ethanol yield' "
                         '(default: the five of the manuscript figure)')
    ap.add_argument('--list', action='store_true',
                    help='list the metrics of the index table and exit')
    ap.add_argument('--dpi', type=int, default=DPI,
                    help='resolution of the figures (default: %(default)s, '
                         "as make_figure.py's figure)")
    ap.add_argument('--formats', nargs='+', default=['png'],
                    help='file formats to write (default: png)')
    args = ap.parse_args()

    table = load_table(args.study_name)
    if args.list:
        for m in dict.fromkeys(table.metric):
            t = table[table.metric == m]
            print(f'{m!r:22} {TITLES.get(m, m):26} Q2 = {t.Q2.iloc[0]:.2f}'
                  + ('' if bool(t.reliable.iloc[0]) else '  (unreliable)'))
        return
    metrics = [resolve_metric(name, table)
               for name in (args.metric or sobol.STACK_METRICS)]
    os.makedirs(OUT_DIR, exist_ok=True)
    for metric in dict.fromkeys(metrics):
        fig = plot_metric(table, metric)
        slug = re.sub(r'[^0-9A-Za-z]+', '-', metric).strip('-')
        for ext in args.formats:
            path = os.path.join(OUT_DIR, f'EDFig2_bars_{slug}.{ext}')
            fig.savefig(path, dpi=args.dpi)
            print('wrote', path)
        plt.close(fig)


if __name__ == '__main__':
    main()
