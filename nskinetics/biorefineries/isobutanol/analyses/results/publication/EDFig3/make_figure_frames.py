# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Progressive presentation frames of Extended Data Figure 3, into
./output/frames.

Runs code/frames/plot_kin_opt_parameter_sets_frames.py --variant process
through code/run_plot.py (the relay-manifest MAX_PATH shim), with
NSKINETICS_IBO_RESULTS_DIR = data/, as make_figure.py runs the plot script.
It draws the figure once with the unmodified snapshot, titles panel A
"Optimization trajectories" as in the manuscript, then saves one PNG per
frame: the empty process-metric cells, each cell's own isobutanol and then
ethanol campaign, the other process-level campaigns' incumbents, the
uninformed and the TRY-informed profitability campaigns' incumbents, then
panel B row by row. The last frame is the complete figure. The frame
sequence is the FRAMES table in that script. Works from any cwd.

Arguments are passed on, e.g. python make_figure_frames.py --dpi 150
--formats png pdf. No __pycache__ is written (PYTHONDONTWRITEBYTECODE=1 in
the child).
"""
import os
import subprocess
import sys

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
PYTHON = sys.executable
CODE = os.path.join(HERE, 'code')
SHIM = os.path.join(CODE, 'run_plot.py')
FRAMES_SCRIPT = os.path.join(CODE, 'frames',
                             'plot_kin_opt_parameter_sets_frames.py')
DATA = os.path.join(HERE, 'data')
VARIANT = 'process'


def main():
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1')
    # the child's exit code is passed on (no second traceback on its errors)
    sys.exit(subprocess.run([PYTHON, SHIM, FRAMES_SCRIPT, '--variant',
                             VARIANT, *sys.argv[1:]], env=env).returncode)


if __name__ == '__main__':
    main()
