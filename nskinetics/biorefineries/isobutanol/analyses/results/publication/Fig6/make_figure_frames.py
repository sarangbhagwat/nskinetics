# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Progressive presentation frames of Figure 6, into ./output/frames.

Runs code/frames/plot_kin_opt_parameter_sets_frames.py --variant irr through
code/run_plot.py (the relay-manifest MAX_PATH shim), with
NSKINETICS_IBO_RESULTS_DIR = data/, as make_figure.py runs the plot script.
It draws the figure once with the unmodified snapshot, titles panel A
"Optimization trajectories" and adds the manuscript's four panel-A
annotations on a canvas extended upwards to hold them, then saves one PNG
per frame: the uninformed profitability campaign, the process-level (TRY)
campaigns, the TRY-informed campaign's preload zone and its own search, each
followed by its annotation, then panel B row by row. The last frame is the
complete figure. The frame sequence is the FRAMES table in that script.
Works from any cwd.

Arguments are passed on, e.g. python make_figure_frames.py --dpi 150
--formats png pdf. No __pycache__ is written (PYTHONDONTWRITEBYTECODE=1 in
the child); the PDF CreationDate is pinned via SOURCE_DATE_EPOCH as in
make_figure.py.
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
VARIANT = 'irr'
SOURCE_DATE_EPOCH = '1792286460'   # as in make_figure.py


def main():
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1',
               SOURCE_DATE_EPOCH=SOURCE_DATE_EPOCH)
    subprocess.run([PYTHON, SHIM, FRAMES_SCRIPT, '--variant', VARIANT,
                    *sys.argv[1:]], env=env, check=True)


if __name__ == '__main__':
    main()
