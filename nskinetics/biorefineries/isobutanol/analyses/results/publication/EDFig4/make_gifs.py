# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Animated version of make_figure.py: build GIFs of the plot panel of
Extended Data Figure 4, for presentations.

One GIF per step towards the full figure (IRR, then the profitability
index, then ethanol yield, then the burden-derating marker), each starting
where the one before ended, so they can be shown one per slide or click:

    output/gifs/k_3_try_informed_sweep_etoh_yield_step00_start.png
    output/gifs/k_3_try_informed_sweep_etoh_yield_stepNN_<what it adds>.gif
    output/gifs/k_3_try_informed_sweep_etoh_yield_stepNN_<what it adds>.png
    output/gifs/k_3_try_informed_sweep_etoh_yield_all_steps.gif
    output/gifs/k_3_try_informed_sweep_etoh_yield_steps.txt

The step GIFs play once and rest on their last frame; the PNGs are the
finished state of each step; all_steps.gif runs the whole build in a loop;
steps.txt lists what each step shows. The last step ends on the complete
figure, as make_figure.py draws it in Arial (where Arial is not installed the
GIFs use a metric-compatible face, make_figure.py the wider DejaVu Sans).

Runs code/animate_figure.py (which draws the figure with the unmodified
snapshot code/plots/plot_axis_try_informed_sweep.py and animates it with
code/gif_engine.py) with the results root pointed at ./data. Sim-safe, works
from any cwd, writes no __pycache__. Options are passed through, e.g.

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
    # drop the empty output folders the plot script creates under ./data
    src_dir = os.path.join(DATA, 'publication', 'Objective-landscape')
    for d in (src_dir, os.path.dirname(src_dir)):
        if os.path.isdir(d) and not os.listdir(d):
            os.rmdir(d)
    sys.exit(done.returncode)   # the child has already reported any error


if __name__ == '__main__':
    main()
