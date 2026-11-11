# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Animated version of make_figure.py: build GIFs of Figure 6 (optimization
trajectories + final proteome allocation), for presentations, with the four
panel-A annotations of the manuscript figure.

One GIF per step towards the full figure (the baseline, the uninformed
profitability campaign, the process-level campaigns' incumbents and their
highest-IRR trials, the scout trials, the TRY-informed campaign, then the
proteome rows), each starting where the one before ended, so they can be
shown one per slide or click:

    output/gifs/kinBO_parameter_sets_irr_step00_start.png
    output/gifs/kinBO_parameter_sets_irr_stepNN_<what it adds>.gif
    output/gifs/kinBO_parameter_sets_irr_stepNN_<what it adds>.png
    output/gifs/kinBO_parameter_sets_irr_all_steps.gif
    output/gifs/kinBO_parameter_sets_irr_steps.txt

The step GIFs play once and rest on their last frame; the PNGs are the
finished state of each step; all_steps.gif runs the whole build in a loop;
steps.txt lists what each step shows. The last step ends on the figure
make_figure.py draws, plus the annotations (with headroom above panel A for
them, the panel title "Optimization trajectories") and less blank space
below panel B.

Runs code/animate_figure.py (which draws the figure with the unmodified
snapshot plot_kin_opt_parameter_sets.py, --variants irr, and animates it with
code/gif_engine.py) through code/run_plot.py, the shim make_figure.py uses
for the relay manifest, with the results root pointed at ./data. Sim-safe,
works from any cwd, writes no __pycache__. Options are passed through, e.g.

    python make_gifs.py --width-px 1920 --speed 1.5 --loop

(python make_gifs.py --help lists them). Needs Pillow besides the packages
make_figure.py needs.
"""
import os
import subprocess
import sys

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
CODE = os.path.join(HERE, 'code')
DATA = os.path.join(HERE, 'data')
SHIM = os.path.join(CODE, 'run_plot.py')
SCRIPT = os.path.join(CODE, 'animate_figure.py')


def main():
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1')
    subprocess.run([PY, SHIM, SCRIPT] + sys.argv[1:], env=env, check=True)


if __name__ == '__main__':
    main()
