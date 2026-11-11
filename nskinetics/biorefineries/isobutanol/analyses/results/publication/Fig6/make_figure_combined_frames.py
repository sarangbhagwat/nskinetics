# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Progressive presentation frames of the optimization trajectories of Figure
6 and Extended Data Figure 3 side by side, into ./output/frames.

Runs code/frames/plot_kin_opt_parameter_sets_frames.py --variant combined
through code/run_plot.py (the relay-manifest MAX_PATH shim), with
NSKINETICS_IBO_RESULTS_DIR = data/, as make_figure_frames.py does for Figure
6 alone. The combined figure is Figure 6's panel A (the two profitability
cells with the manuscript's four annotations, and the campaign legend under
them) with Extended Data Figure 3's panel A (the six process-metric cells) to
its right, both at the size they have in their own figures; no panel title,
no panel B. This folder's snapshot and campaign CSVs draw both (they are the
ones Extended Data Figure 3 uses). One PNG per frame, all on the same canvas:
the baseline, the uninformed profitability campaign, the six process-level
(TRY) campaigns in the metrics they optimize and then in profitability, their
trials as the TRY-informed campaign's preload, and that campaign's own
search, with each annotation after what it describes. The last frame is the
complete figure. The frame sequence is FRAMES['combined'] in that script.
Works from any cwd.

Arguments are passed on, e.g. python make_figure_combined_frames.py --dpi 150
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
VARIANT = 'combined'
SOURCE_DATE_EPOCH = '1792286460'   # as in make_figure.py


def main():
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1',
               SOURCE_DATE_EPOCH=SOURCE_DATE_EPOCH)
    # the child's exit code is passed on (no second traceback on its errors);
    # --variant comes last, so it cannot be overridden from the command line
    sys.exit(subprocess.run([PYTHON, SHIM, FRAMES_SCRIPT, *sys.argv[1:],
                             '--variant', VARIANT], env=env).returncode)


if __name__ == '__main__':
    main()
