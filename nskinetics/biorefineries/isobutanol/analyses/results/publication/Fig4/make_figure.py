# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Regenerate manuscript Figure 4 (Monte Carlo bivariate KDE panels A-F, MESP vs
TCI, ethanol production, DDGS revenue, ethanol yield, titer, productivity).

Sim-safe: runs the snapshot code/plots/plot_uncertainty_MPSP_vs_TCI.py (which
loads only code/paths.py by file path) on the 6,000-sim scenario-A workbooks in
data/, and writes the PNG + PDF into output/. Works from any cwd.
"""
import os
import sys
import glob
import shutil
import subprocess

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
PY = r"C:\Users\saran\anaconda3\envs\IBO_2026\python.exe"
if not os.path.exists(PY):
    PY = sys.executable
TAG = "_IBO_2026.9.29-23.45_['A']_6000sims_A"
SCRIPT = os.path.join(HERE, 'code', 'plots', 'plot_uncertainty_MPSP_vs_TCI.py')
DATA = os.path.join(HERE, 'data')
OUT = os.path.join(HERE, 'output')


def main():
    # PYTHONDONTWRITEBYTECODE: the child loads code/paths.py by file path, which
    # would otherwise leave code/__pycache__/ inside this folder.
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1')
    os.makedirs(OUT, exist_ok=True)
    tmp = os.path.join(HERE, 'output', '_tmp')
    os.makedirs(tmp, exist_ok=True)
    subprocess.run([PY, SCRIPT,
                    os.path.join(DATA, TAG + '_1_full_evaluation.xlsx'),
                    '--baseline-file', os.path.join(DATA, TAG + '_0_baseline.xlsx'),
                    '--out-dir', tmp],
                   env=env, cwd=HERE, check=True)
    for f in glob.glob(os.path.join(tmp, '*')):
        shutil.copy2(f, OUT)
        print(f'-> {os.path.join(OUT, os.path.basename(f))}')
    shutil.rmtree(tmp)


if __name__ == '__main__':
    main()
