# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Progressive presentation frames of Extended Data Figure 4, into
./output/frames.

Runs code/frames/plot_k_3_try_informed_sweep_frames.py with the results root
pointed at ./data (as make_figure.py does). It draws the figure once with the
snapshot code/plots/plot_axis_try_informed_sweep.py --param k_3, applies the
manuscript's edits (legend removed, rotated "growth derated at higher k3
values" label added) and saves one PNG per frame on the same canvas: the x
axis, then the ethanol-yield, profitability-index and IRR axes with their
curves, then the burden-derating line and its label. The last frame is the
manuscript figure. The frame sequence is the FRAMES table at the top of that
script. Works from any cwd.

Arguments are passed on, e.g. python make_figure_frames.py --dpi 150
--formats png pdf. The child runs with PYTHONDONTWRITEBYTECODE=1 (no
__pycache__ under code/) and the fixed SOURCE_DATE_EPOCH of make_figure.py;
the plot script's empty output folders under ./data are removed again.
"""
import os
import sys
import subprocess

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
DATA = os.path.join(HERE, 'data')

def main():
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1',
               # fixed PDF CreationDate (--formats pdf), as in make_figure.py
               SOURCE_DATE_EPOCH='1790294400')
    script = os.path.join(HERE, 'code', 'frames',
                          'plot_k_3_try_informed_sweep_frames.py')
    # the child's exit code is passed on (no second traceback on its errors)
    sys.exit(subprocess.run([PY, script, *sys.argv[1:]], env=env).returncode)

if __name__ == '__main__':
    main()
