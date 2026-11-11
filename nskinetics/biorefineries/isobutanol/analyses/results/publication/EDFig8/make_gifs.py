# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Animated version of make_figure.py: build GIFs of Extended Data Figure 8
(k13 x isobutanol-inhibition IRR contour, opt_PI_TRY_informed), for
presentations.

One GIF per step towards the full figure (the IRR map, the break-even
contour, then the optima of the ethanol metrics, the isobutanol metrics, cell
density and TCI, and last profitability), each starting where the one before
ended, so they can be shown one per slide or click:

    output/gifs/k13_inhib_isobutanol_opt_PI_TRY_informed_IRR_step00_start.png
    output/gifs/k13_inhib_isobutanol_opt_PI_TRY_informed_IRR_stepNN_<what it adds>.gif
    output/gifs/k13_inhib_isobutanol_opt_PI_TRY_informed_IRR_stepNN_<what it adds>.png
    output/gifs/k13_inhib_isobutanol_opt_PI_TRY_informed_IRR_all_steps.gif
    output/gifs/k13_inhib_isobutanol_opt_PI_TRY_informed_IRR_steps.txt

The step GIFs play once and rest on their last frame; the PNGs are the
finished state of each step; all_steps.gif runs the whole build in a loop;
steps.txt lists what each step shows. The last step ends on the complete
figure, as make_figure.py draws it (in Arial, or a metric-compatible face
where Arial is not installed).

Runs code/animate_figure.py (which draws the figure with the unmodified
snapshot code/plots/plot_k13_inhib_isobutanol_opt_PI_TRY_informed.py and
animates it with code/gif_engine.py) with the results root pointed at ./data.
Sim-safe, works from any cwd, writes no __pycache__. Options are passed
through, e.g.

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
    subprocess.run([PY, SCRIPT] + sys.argv[1:], env=env, check=True)
    # drop the empty output folders the plot script creates under ./data
    src_dir = os.path.join(DATA, 'publication', 'Kinetic-sweeps')
    for d in (src_dir, os.path.dirname(src_dir)):
        if os.path.isdir(d) and not os.listdir(d):
            os.rmdir(d)


if __name__ == '__main__':
    main()
