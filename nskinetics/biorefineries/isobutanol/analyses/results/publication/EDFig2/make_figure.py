# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Regenerate Extended Data Figure 2 (Sobol'/Shapley bars, five outputs).

Runs the snapshot code/analyses/sensitivity/analyze_sobol_split12d.py in its
``--replot`` mode (sim-safe: redraws from the saved Sobol' index table, no
surrogate refit, no simulation) with the results root pointed at ./data, then
moves every figure it writes into ./output (the script writes them next to
its input, i.e. transiently into ./data). Works from any cwd. Bytecode writing
is disabled so no __pycache__ appears in this folder.
"""
import os
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True  # keep this folder free of __pycache__

HERE = os.path.dirname(os.path.abspath(__file__))
PY = r'C:\Users\saran\anaconda3\envs\IBO_2026\python.exe'
STUDY = ('kin_sobol_ethanol_isobutanol_metabolic_split_12d_custom9f1df84b'
         '_aA_burden_seed20260922')
DATA = os.path.join(HERE, 'data')
OUT = os.path.join(HERE, 'output')

def main():
    # PYTHONDONTWRITEBYTECODE: the by-path sibling loads would otherwise write
    # code/__pycache__/*.pyc into this folder
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1')
    script = os.path.join(HERE, 'code', 'analyses', 'sensitivity', 'analyze_sobol_split12d.py')
    subprocess.run([PY, script, '--study-name', STUDY, '--replot'], env=env, check=True)
    os.makedirs(OUT, exist_ok=True)
    # --replot writes every figure next to its input; move them all to ./output
    # (the stacked bars are the manuscript figure; PI bars / heatmap are by-products)
    for name in sorted(os.listdir(DATA)):
        if name.startswith(STUDY + '_sobol_') and name.endswith(('.png', '.pdf')):
            dst = os.path.join(OUT, 'EDFig2_' + name[len(STUDY + '_sobol_'):])
            shutil.move(os.path.join(DATA, name), dst)
            print('wrote', dst)

if __name__ == '__main__':
    main()
