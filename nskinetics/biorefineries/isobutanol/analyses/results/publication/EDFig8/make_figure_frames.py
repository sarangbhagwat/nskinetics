# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Regenerate Extended Data Figure 8 (k13 x isobutanol-inhibition IRR contour, opt_PI_TRY_informed) as progressive presentation frames, from the snapshot in code/ and the sweep CSVs in data/.

Writes output/frames/k13_inhib_isobutanol_opt_PI_TRY_informed_IRR_frameNN_<name>.png:
cumulative frames of one drawing of the figure (axes; IRR map + colour bar;
break-even line; then the optima, a few at a time), all with the same size
and crop, so they overlay exactly on consecutive slides; the last frame is
the figure make_figure.py saves. The sequence is the FRAMES table of
code/frames/plot_k13_inhib_isobutanol_opt_PI_TRY_informed_frames.py.
Arguments are passed on to that script: --dpi N (default 600) and
--formats png [pdf svg ...] (default png).

Sim-safe: the figure script reads CSVs only and never imports nskinetics,
biosteam, thermosteam or biorefineries. Works from any cwd. Bytecode writing is
disabled (PYTHONDONTWRITEBYTECODE=1
in the child), so no __pycache__ appears under code/.
"""

import os
import subprocess
import sys

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
PYTHON = sys.executable
SCRIPT = os.path.join(HERE, 'code', 'frames', 'plot_k13_inhib_isobutanol_opt_PI_TRY_informed_frames.py')
DATA = os.path.join(HERE, 'data')


def main():
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1')
    # the child's exit code is passed on (no second traceback on its errors)
    sys.exit(subprocess.run([PYTHON, SCRIPT, *sys.argv[1:]], env=env,
                            cwd=os.path.dirname(SCRIPT)).returncode)


if __name__ == '__main__':
    main()
