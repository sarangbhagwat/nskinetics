# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Regenerate manuscript Extended Data Figure 6 (Spearman rho tornado for MESP, TCI and
ethanol titer, productivity, yield; Monte Carlo, 6,000 sims).

Sim-safe: runs the snapshot code/plots/plot_uncertainty_spearman_tornado.py (which
loads only code/paths.py by file path) on the 6,000-sim scenario-A workbook in
data/, and writes the PNG + PDF + CSV into output/. Works from any cwd.
"""
import os
import sys
import glob
import shutil
import subprocess

sys.dont_write_bytecode = True  # keep this folder free of __pycache__

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
TAG = "_IBO_2026.9.29-23.45_['A']_6000sims_A"
SCRIPT = os.path.join(HERE, 'code', 'plots', 'plot_uncertainty_spearman_tornado.py')
DATA = os.path.join(HERE, 'data')
OUT = os.path.join(HERE, 'output')


def main():
    # PYTHONDONTWRITEBYTECODE: the plot script loads code/paths.py by file path,
    # which would otherwise drop a __pycache__ into code/.
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1')
    os.makedirs(OUT, exist_ok=True)
    tmp = os.path.join(HERE, 'output', '_tmp')
    os.makedirs(tmp, exist_ok=True)
    subprocess.run([PY, SCRIPT,
                    os.path.join(DATA, TAG + '_1_full_evaluation.xlsx'),
                    '--out-dir', tmp],
                   env=env, cwd=HERE, check=True)
    for f in glob.glob(os.path.join(tmp, '*')):
        shutil.copy2(f, OUT)
        print(f'-> {os.path.join(OUT, os.path.basename(f))}')
    shutil.rmtree(tmp)


if __name__ == '__main__':
    main()
