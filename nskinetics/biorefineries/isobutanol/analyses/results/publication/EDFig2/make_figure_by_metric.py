# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""The Sobol'/Shapley bars of Extended Data Figure 2 one metric per figure,
with the decision variables grouped on the x axis, into ./output/by_metric.

Runs code/by_metric/plot_sobol_bars_by_metric.py with the results root pointed
at ./data (as make_figure.py does). It draws, from the same saved Sobol' index
table and with the bars of the snapshot
code/analyses/sensitivity/analyze_sobol_split12d.py, one figure per metric:
the six enzyme expression levels (glycolysis, ethanol production, isobutanol
production), the three inhibition levels and the three substrate feeding
variables, in that fixed order under group headers. Every figure has the same
canvas, so the metrics overlay exactly on consecutive slides. The grouping is
the GROUPS table at the top of that script. Works from any cwd.

    python make_figure_by_metric.py                       # the five metrics of ED Fig. 2
    python make_figure_by_metric.py --metric "EtOH yield" # one metric (or several names)
    python make_figure_by_metric.py --list                # the metrics that can be plotted

Arguments are passed on (also --dpi N, --formats png pdf). Sim-safe: no
surrogate refit, no simulation. The child runs with PYTHONDONTWRITEBYTECODE=1,
so no __pycache__ appears under code/.
"""
import os
import subprocess
import sys

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
STUDY = ('kin_sobol_ethanol_isobutanol_metabolic_split_12d_custom9f1df84b'
         '_aA_burden_seed20260922')   # as in make_figure.py
DATA = os.path.join(HERE, 'data')


def main():
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1')
    script = os.path.join(HERE, 'code', 'by_metric',
                          'plot_sobol_bars_by_metric.py')
    # the child's exit code is passed on (no second traceback on its errors);
    # --study-name comes last, so it cannot be overridden from the command line
    sys.exit(subprocess.run([PY, script, *sys.argv[1:], '--study-name', STUDY],
                            env=env).returncode)


if __name__ == '__main__':
    main()
