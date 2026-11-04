# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Regenerate Extended Data Figure 8 (k13 x isobutanol-inhibition IRR contour, opt_PI_TRY_informed) from the snapshot in code/ and the sweep CSVs in data/.

Sim-safe: the figure script reads CSVs only and never imports nskinetics,
biosteam, thermosteam or biorefineries. Works from any cwd. Bytecode writing is
disabled (PYTHONDONTWRITEBYTECODE=1
in the child), so no __pycache__ appears under code/.
"""

import os
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
PYTHON = r'C:\Users\saran\anaconda3\envs\IBO_2026\python.exe'
if not os.path.exists(PYTHON):
    PYTHON = sys.executable
SCRIPT = os.path.join(HERE, 'code', 'plots', 'plot_k13_inhib_isobutanol_opt_PI_TRY_informed.py')
DATA = os.path.join(HERE, 'data')
PRODUCED = os.path.join(DATA, 'publication', 'Kinetic-sweeps')
OUTPUT = os.path.join(HERE, 'output')
STEMS = ['k13_inhib_isobutanol_opt_PI_TRY_informed_IRR']


def main():
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1')
    subprocess.run([PYTHON, SCRIPT], env=env, cwd=os.path.dirname(SCRIPT), check=True)
    os.makedirs(OUTPUT, exist_ok=True)
    for stem in STEMS:
        for ext in ('png', 'pdf'):
            src = os.path.join(PRODUCED, f'{stem}.{ext}')
            shutil.copy2(src, os.path.join(OUTPUT, f'{stem}.{ext}'))
            print('copied', os.path.join(OUTPUT, f'{stem}.{ext}'))


if __name__ == '__main__':
    main()
