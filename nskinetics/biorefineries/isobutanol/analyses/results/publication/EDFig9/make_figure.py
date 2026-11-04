# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Regenerate EDFig9 (trajectory of the TRY-informed optimum) from the snapshot in code/ and the inputs in data/.

Sim-safe: the script runs from its trajectory cache in data/ (stage 2 only);
it never imports nskinetics, biosteam, thermosteam or biorefineries. It refuses
to run if the cache is missing (that would trigger the simulating stage 1).
Works from any working directory.
"""
import os
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True  # keep __pycache__ out of this folder

HERE = os.path.dirname(os.path.abspath(__file__))
# Hardcoded on purpose (the IBO_2026 env); edit this line to port elsewhere.
PYTHON = sys.executable
SCRIPT = os.path.join(HERE, 'code', 'plots', 'plot_trajectory_opt_PI_TRY_informed.py')
DATA = os.path.join(HERE, 'data')
OUTPUT = os.path.join(HERE, 'output')
THEME = 'Kinetic-trajectories'
STEM = 'trajectory_opt_PI_TRY_informed'
CACHE = os.path.join(DATA, 'trajectory_opt_PI_TRY_informed')


def main():
    for ext in ('.csv', '.json'):
        if not os.path.exists(CACHE + ext):
            sys.exit(f'missing trajectory cache {CACHE}{ext}: refusing to run '
                     '(the script would simulate)')
    # PYTHONDONTWRITEBYTECODE: the child loads code/ siblings by path; write no __pycache__.
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1')
    subprocess.run([PYTHON, SCRIPT], env=env, cwd=HERE, check=True)
    os.makedirs(OUTPUT, exist_ok=True)
    for ext in ('png', 'pdf'):
        src = os.path.join(DATA, 'publication', THEME, f'{STEM}.{ext}')
        shutil.copy2(src, os.path.join(OUTPUT, f'{STEM}.{ext}'))
        print(f'Copied {src} -> output/')


if __name__ == '__main__':
    main()
