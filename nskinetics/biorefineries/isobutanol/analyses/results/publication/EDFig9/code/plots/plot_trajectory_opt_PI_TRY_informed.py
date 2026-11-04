# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Concentration vs time of the TRY-informed profitability optimum.

The V406 fed-batch fermentation of scenario ``opt_PI_TRY_informed`` (trial
#1912 of the PI (log-tail) relay campaign; scenarios.SCENARIOS feeding
strategy, enzyme burden on), in two panels:

A  cell density, active cell density, glucose and the products (ethanol,
   isobutanol, acetate), with the glucose spikes and the ends of aeration
   and of fermentation marked.
B  the pathway intermediates (pyruvate, acetaldehyde, and the Ehrlich-branch
   acetolactate, DHI, KIV and isobutyraldehyde) on their own scale.

Two stages, as in plot_feeding_strategy_baseline_A.py:
  1. (simulation) one isobutanol.load(separation_processes=('IBO_EtOH',))
     (smoke case 10's configuration) + scenarios.load_scenario(SCENARIO);
     the MPSPs are checked against the scenario pins and the V406 time
     course is cached to
     <results root>/trajectory_opt_PI_TRY_informed.{csv,json}.
     Runs only when that cache is missing or with --resimulate.
  2. (sim-safe) reads the cached trajectory and draws.

Output: <results root>/publication/Kinetic-trajectories/trajectory_opt_PI_TRY_informed.{png,pdf}
"""

import argparse
import importlib.util
import json
import math
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import AutoMinorLocator, MaxNLocator

#%% Paths and settings

HERE = os.path.dirname(os.path.abspath(__file__))


def _load_by_path(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


paths = _load_by_path('_ibo_paths', os.path.join(os.path.dirname(HERE), 'paths.py'))
# house style (fonts, tick styling, axis titles) of the scenario-A trajectory
# panel; sim-safe to load by path (it imports nothing from the package)
_style = _load_by_path('_ibo_feed_strat_A',
                       os.path.join(HERE, 'plot_feeding_strategy_baseline_A.py'))

SCENARIO = 'opt_PI_TRY_informed'
SEPARATION_PROCESSES = ('IBO_EtOH',)
RESULTS_DIR = paths.results_dir()
OUTPUT_DIR = paths.publication_dir('Kinetic-trajectories')
OUTPUT_STEM = f'trajectory_{SCENARIO}'
TRAJECTORY_STEM = os.path.join(RESULTS_DIR, f'trajectory_{SCENARIO}')

# (nskinetics results key, legend label); default colour cycle
MAIN_SERIES = [
    ('[x]', 'cell density'),
    ('curr_a', 'active cell density'),
    ('[s_glu]', 'glucose'),
    ('[s_EtOH]', 'ethanol'),
    ('[s_IBO]', 'isobutanol'),
    ('[s_acetate]', 'acetate'),
]
INTERMEDIATE_SERIES = [
    ('[s_pyr]', 'pyruvate'),
    ('[s_acetald]', 'acetaldehyde'),
    ('[s_AL]', 'acetolactate'),
    ('[s_DHI]', 'DHI'),
    ('[s_KIV]', 'KIV'),
    ('[s_isobutald]', 'isobutyraldehyde'),
]
SPIKE_KEY = 'curr_n_glu_spikes'

FONTS = _style.FONTS
G_PER_L = _style.G_PER_L
bold_title = _style.bold_title

#%% Stage 1: the fermentation time course (simulation)

def simulate_trajectory():
    """One load + the scenario baseline; checks the pins, caches V406's
    time course (every series of nsk_results_dict on the time grid)."""
    from nskinetics.biorefineries import isobutanol
    isobutanol.load(separation_processes=SEPARATION_PROCESSES)
    from nskinetics.biorefineries.isobutanol import scenarios
    bundle = scenarios.load_scenario(SCENARIO)  # runs one baseline simulation
    results = bundle['solve_TEA'](stream_IDs=('ethanol', 'isobutanol'))
    scenarios.assert_MPSPs_pinned(bundle['expected'], results, 1)
    scenarios.assert_spike_feed_residual(1)
    V406 = bundle['V406']
    d = V406.nsk_results_dict
    time = np.asarray(d['time'], dtype=float)
    columns = {'time': time}
    for key, value in d.items():
        if key == 'time':
            continue
        try:
            arr = np.asarray(value, dtype=float)
        except (TypeError, ValueError):
            continue
        if arr.shape == time.shape:
            columns[key] = arr
    df = pd.DataFrame(columns)
    df.to_csv(TRAJECTORY_STEM + '.csv', index=False)
    specific = V406.nsk_results_specific_tau_dict
    meta = {
        'scenario': SCENARIO,
        'separation_processes': list(SEPARATION_PROCESSES),
        'feeding_kwargs': bundle['feeding_kwargs'],
        'max_n_spikes': int(bundle['fbs_spec'].max_n_spikes),
        'burden_on': bool(bundle['burden_on']),
        'n_glu_spikes': float(specific['curr_n_glu_spikes']),
        'tau': float(V406.tau),
        'tau_stop_aeration': float(V406.tau_stop_aeration),
        'EtOH_titer': float(specific['[s_EtOH]']),
        'IBO_titer': float(specific['[s_IBO]']),
        'cell_density': float(specific['[x]']),
        'MPSPs': {k: float(v) for k, v in results['MPSPs'].items()},
        'IRR': float(results['IRR']),
    }
    with open(TRAJECTORY_STEM + '.json', 'w') as file:
        json.dump(meta, file, indent=2)
    print(f'Cached the trajectory: {TRAJECTORY_STEM}.csv/.json')
    return df, meta


def load_trajectory(resimulate=False):
    if resimulate or not (os.path.exists(TRAJECTORY_STEM + '.csv')
                          and os.path.exists(TRAJECTORY_STEM + '.json')):
        return simulate_trajectory()
    with open(TRAJECTORY_STEM + '.json') as file:
        meta = json.load(file)
    print(f'Using the cached trajectory {TRAJECTORY_STEM}.csv '
          '(--resimulate to refresh).')
    return pd.read_csv(TRAJECTORY_STEM + '.csv'), meta

#%% Figure

def spike_times(df):
    """Times at which the cumulative glucose-spike count steps up."""
    if SPIKE_KEY not in df:
        return np.array([])
    n = df[SPIKE_KEY].to_numpy()
    return df['time'].to_numpy()[1:][np.diff(n) > 0.5]


def nice_upper(value, step):
    return step * math.ceil(value / step + 1e-9)


def draw_series(ax, df, series, t_max):
    shown = df['time'] <= t_max
    for key, label in series:
        if key in df:
            ax.plot(df['time'][shown], df[key][shown], label=label,
                    linewidth=1.6)
    ax.set_xlim(0., t_max)
    ax.set_ylim(bottom=0.)
    ax.xaxis.set_major_locator(MaxNLocator(6, steps=[1, 2, 2.5, 5, 10]))
    ax.yaxis.set_major_locator(MaxNLocator(6, steps=[1, 2, 2.5, 5, 10]))
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_minor_locator(AutoMinorLocator(2))
    ax.set_xlabel(bold_title('Time', 'h'), fontsize=FONTS['axis_title'])
    ax.set_ylabel(bold_title('Concentration', G_PER_L),
                  fontsize=FONTS['axis_title'])
    ax.legend(loc='upper right', fontsize=FONTS['legend'] - 1, frameon=True,
              fancybox=False, edgecolor='black', framealpha=1.,
              handlelength=2.2, borderaxespad=0.6)


def mark_events(ax, df, meta, label=True):
    top = ax.get_ylim()[1]
    # (time, text, colour, line style, label side: +1 right / -1 left of
    # the line); the spike precedes the end of aeration, so their labels
    # sit on opposite sides
    events = [(t, 'glucose spike', 'tab:blue', ':', -1)
              for t in spike_times(df)]
    events += [(meta['tau_stop_aeration'], 'end of aeration', 'grey', '--', 1),
               (meta['tau'], 'end of fermentation', 'grey', '--', -1)]
    for t, text, color, style, side in events:
        ax.axvline(t, color=color, linestyle=style, linewidth=1.0, zorder=0)
        if label:
            ax.annotate(text, xy=(t, top), xytext=(4 * side, 10),
                        textcoords='offset points',
                        ha='left' if side > 0 else 'right', va='bottom',
                        fontsize=FONTS['annotation'] - 1,
                        color='tab:blue' if color == 'tab:blue' else 'black',
                        annotation_clip=False,
                        arrowprops=dict(arrowstyle='-|>', mutation_scale=9,
                                        color=color if color == 'tab:blue' else 'black',
                                        lw=0.9, shrinkA=1, shrinkB=1,
                                        connectionstyle=f'arc3,rad={0.4 * side}'))


def make_figure(df, meta):
    _style.apply_font_rcparams()
    t_max = nice_upper(meta['tau'] * 1.05, 10.)
    fig = plt.figure(figsize=(11.0, 4.4))
    gs = GridSpec(1, 2, figure=fig, wspace=0.28, left=0.07, right=0.98,
                  bottom=0.14, top=0.84)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    draw_series(axA, df, MAIN_SERIES, t_max)
    draw_series(axB, df, INTERMEDIATE_SERIES, t_max)
    mark_events(axA, df, meta, label=True)
    mark_events(axB, df, meta, label=False)
    for ax, letter in ((axA, 'A'), (axB, 'B')):
        ax.text(-0.14, 1.09, letter, transform=ax.transAxes, ha='left',
                va='bottom', fontsize=FONTS['panel_letter'], fontweight='bold')
    fig.canvas.draw()
    _style.style_ticks(axA)
    _style.style_ticks(axB)
    return fig


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--resimulate', action='store_true',
                        help='re-run the scenario simulation even if the '
                             'trajectory cache exists')
    args = parser.parse_args()

    df, meta = load_trajectory(args.resimulate)
    fig = make_figure(df, meta)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ('png', 'pdf'):
        path = os.path.join(OUTPUT_DIR, f'{OUTPUT_STEM}.{ext}')
        fig.savefig(path, dpi=600, bbox_inches='tight', facecolor='white')
        print(f'Saved {path}')
    plt.close(fig)

    print(f"\n{SCENARIO}: feeding {meta['feeding_kwargs']}, "
          f"{meta['n_glu_spikes']:g} spike(s), tau {meta['tau']:.2f} h, end "
          f"of aeration {meta['tau_stop_aeration']:.2f} h, ethanol "
          f"{meta['EtOH_titer']:.1f} g/L, isobutanol {meta['IBO_titer']:.1f} "
          f"g/L, cells {meta['cell_density']:.2f} g/L, MPSPs {meta['MPSPs']}, "
          f"IRR {meta['IRR']:.4f}")


if __name__ == '__main__':
    main()
