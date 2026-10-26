# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""k_13 x isobutanol-inhibition sweep figure for the opt_PI_TRY_informed scenario.

IRR at default prices over the k_13 (ALS capacity) x grouped inhib_isobutanol
multiplier grid of analyses/kinetic_sweeps/evaluate_EtOH_k13_inhib_isobutanol.py
(opt_PI_TRY_informed baseline = trial #1912 of the 2026-09-24 TRY-informed PI
(log-tail) relay campaign _rl15c111dc; enzyme burden ON, the scenario's
feeding strategy (1-spike cap / 189.19 / 194.19), no feeding-strategy
optimization), with the grid optimum of IRR, TCI, cell density and the ethanol
and isobutanol titer / yield / productivity marked and annotated. Format and
style follow plots/plot_glycolysis_inhib_ethanol_A.py (its IRR twin); the
drawing lives in plots/_kinetic_sweep_figure.py (loaded by file path).

Renamed 2026-10-08 from plot_k13_inhib_isobutanol_opt_IRR.py when the opt_IRR
scenario (split_12d trial 1602) was removed in favour of opt_PI_TRY_informed.
Inputs: the 2026-10-26 run of the sweep on opt_PI_TRY_informed under the
campaigns' hensmith 2b5b27d (max_n=1 CSV prefix), so the grid's TEA metrics
are comparable with the campaign record. The old opt_IRR figure
(k13_inhib_isobutanol_opt_IRR_IRR.*, max_n=50 CSVs) is superseded.

Sim-safe: reads the sweep CSVs only and never imports the package.

Output: analyses/results/publication/Kinetic-sweeps/k13_inhib_isobutanol_opt_PI_TRY_informed_IRR.{png,pdf}
"""

import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    '_kinetic_sweep_figure', os.path.join(HERE, '_kinetic_sweep_figure.py'))
ksf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ksf)

# Grid of the sweep to read: the 40 x 40 run of
# analyses/kinetic_sweeps/evaluate_EtOH_k13_inhib_isobutanol.py on
# opt_PI_TRY_informed (2026-10-26, hensmith 2b5b27d), k_13 over 1e-3x-4x the
# optimum's 4.0 g/g/h and the inhib_isobutanol multiplier over 0.75x-1.5x the
# wild-type (scenario A) coefficients (the campaigns' group band; the optimum
# sits at 0.75x). SPEC_1 / SPEC_2 must match that script's linspaces.
SWEEP_STEPS = (40, 40, 1)
BASELINE_K13 = 4.0
SPEC_1 = np.linspace(1e-3 * BASELINE_K13, 4.0 * BASELINE_K13, SWEEP_STEPS[0])  # k_13 [g/g/h]
SPEC_2 = np.linspace(0.75, 1.5, SWEEP_STEPS[1])  # inhib_isobutanol multiplier [-]

# Optimum markers: shape = metric (as in the glycolysis figure), face colour =
# product (ethanol white, isobutanol grey; economic optima cyan). Offsets set
# by eye on the 2026-10-26 opt_PI_TRY_informed grid; retune them if the optima
# move.
IBO_FACE = '#9a9a9a'
IBO_TITER = 'isobutanol\ntiter'
IBO_YIELD = 'isobutanol\nyield'
IBO_PRODUCTIVITY = 'isobutanol\nproductivity'
ETOH_PRODUCTIVITY = 'ethanol\nproductivity'
OPTIMA = [
    ('EtOH Titer',        'max', 'ethanol titer',            '^', 'white',   10, (38, -40),   0.),
    ('EtOH Yield',        'max', 'ethanol yield',            'p', 'white',   10, (26, -16),   0.),
    ('EtOH Productivity', 'max', ETOH_PRODUCTIVITY,          's', 'white',    9, (11, 76),    0.),
    ('IBO Titer',         'max', IBO_TITER,                  '^', IBO_FACE,  10, (4, 35),     0.),
    ('IBO Yield',         'max', IBO_YIELD,                  'p', IBO_FACE,  10, (18, 70),    0.),
    ('IBO Productivity',  'max', IBO_PRODUCTIVITY,           's', IBO_FACE,   9, (24, 116),   0.),
    ('Cell loading',      'max', 'cell density',             'o', 'white',   10, (25, 125),   0.),
    ('TCI',               'min', 'TCI',                      'p', '#33ccff', 10, (-30, -16), -0.2),
    ('IRR',               'max', 'profitability',            '*', '#33ccff', 14, (-20, 40),   0.),
]
# the TCI label sits on the dark (IRR < 0) cells: it and its arrow in white
LABEL_COLOR = {'TCI': 'white'}

FIGURE = ksf.SweepFigure(
    output_stem='k13_inhib_isobutanol_opt_PI_TRY_informed_IRR',
    # max_n=1: opt_PI_TRY_informed's 1-spike cap (the old opt_IRR CSVs carry
    # max_n=50), so stale opt_IRR CSVs are never read by mistake
    csv_prefix=f'ibo_{SWEEP_STEPS}_k_13_inhib_Spike_opt=False_max_n=1_',
    spec_1=SPEC_1, spec_2=SPEC_2,
    produced_by='the k_13 x inhib_isobutanol sweep on opt_PI_TRY_informed (analyses/kinetic_sweeps/evaluate_EtOH_k13_inhib_isobutanol.py)',
    # r13 = k_13·pyr/(K_13 + pyr)·a: k_13 is per gram of active cells
    xlabel=r'$\bf{ALS\ capacity}$, $k_{\mathrm{13}}$' + f' [{ksf.G_PER_G_PER_H}]',
    ylabel=r'$\bf{Isobutanol\ inhibition}$ [× baseline]',
    optima=OPTIMA,
    xlim=(0., 16.), ylim=(0.75, 1.5),
    x_major=4., x_minor_div=4, y_major=0.25, y_minor_div=5,
    # co-located optima fanned apart [pt]: ethanol titer and yield share the
    # top cell next to the left edge, isobutanol titer and yield one cell
    # above the bottom row; cell density and ethanol productivity sit in the
    # two bottom-left cells
    marker_nudge={'ethanol yield': (10., 0.),
                  IBO_TITER: (7., 0.), IBO_YIELD: (-7., 0.),
                  'cell density': (-3., 0.), ETOH_PRODUCTIVITY: (4., 0.)},
    label_color=LABEL_COLOR,
    # 'profitability' centred over its star
    label_ha={'profitability': 'center'},
    # each arrow leaves its label's left edge, ordered so no two cross: the
    # isobutanol column right of the bottom cluster (productivity on top), the
    # stacked left-hand labels to the bottom-left and top-left corners
    arrow_relpos={IBO_TITER: (0., 0.5), IBO_YIELD: (0., 0.5), IBO_PRODUCTIVITY: (0., 0.5),
                  'ethanol titer': (0., 0.5), 'ethanol yield': (0., 0.5),
                  ETOH_PRODUCTIVITY: (0., 0.5), 'cell density': (0., 0.5),
                  'TCI': (1., 0.5)},
    color_axis='IRR',
    # IRR axis -15 to 30 %: the 2026-10-26 grid tops out at 27.2 % (the cell
    # next to the optimum; the campaign's 27.26 %, hensmith 2b5b27d)
    levels=np.arange(-15., 30.00001, 0.5),
    cbar_ticks=np.arange(-15., 30.00001, 5.),
    cbar_minor_step=1.,
)

if __name__ == '__main__':
    ksf.main(FIGURE)
