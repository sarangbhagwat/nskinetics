# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Animated version of make_figure.py: build GIFs of Figure 4 (Monte Carlo
MESP against TCI, ethanol production, DDGS revenue, ethanol yield, titer and
productivity), for presentations.

One GIF per step towards the full figure (panel A piece by piece: the
baseline, the sample density, the box plots, the Pareto frontier, the
conventional corn ethanol box and the gasoline price band; then the samples
and binned-median line of panels B to F, one panel per step; last their grey
ranges), each starting where the one before ended, so they can be shown one
per slide or click:

    output/gifs/uncertainty_bivariate_A_2026.9.29-23.45_6000sims_step00_start.png
    output/gifs/uncertainty_bivariate_A_2026.9.29-23.45_6000sims_stepNN_<what it adds>.gif
    output/gifs/uncertainty_bivariate_A_2026.9.29-23.45_6000sims_stepNN_<what it adds>.png
    output/gifs/uncertainty_bivariate_A_2026.9.29-23.45_6000sims_all_steps.gif
    output/gifs/uncertainty_bivariate_A_2026.9.29-23.45_6000sims_steps.txt

The step GIFs play once and rest on their last frame; the PNGs are the
finished state of each step; all_steps.gif runs the whole build in a loop;
steps.txt lists what each step shows. The last step ends on the complete
figure, as make_figure.py draws it in Arial (where Arial is not installed the
GIFs use a metric-compatible face, make_figure.py the wider DejaVu Sans).

Runs code/animate_figure.py (which draws the figure with the unmodified
snapshot code/plots/plot_uncertainty_MPSP_vs_TCI.py on the 6,000-sim
scenario-A workbooks and animates it with code/gif_engine.py) with the
results root pointed at ./data. Sim-safe, works from any cwd, writes no
__pycache__. Options are passed through, e.g.

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
DATA = os.path.join(HERE, 'data')
SCRIPT = os.path.join(HERE, 'code', 'animate_figure.py')


def main():
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1')
    done = subprocess.run([PY, SCRIPT] + sys.argv[1:], env=env)
    sys.exit(done.returncode)   # the child has already reported any error


if __name__ == '__main__':
    main()
