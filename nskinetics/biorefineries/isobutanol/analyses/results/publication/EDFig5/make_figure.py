# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Regenerate Extended Data Figure 5 (3x3 TEA cost / utility breakdowns of
the scenario-A baseline and the eight split_12d campaign optima).

Runs stage 2 only, the snapshot code/plots/plot_tea_breakdowns_split12d.py
(sim-safe: reads the stage-1 JSON in ./data), writing into ./output. The
script always appends a render timestamp to the figure stem, so the newest
render is renamed to output/EDFig5.png / .pdf. Stage 1
(analyses/tea_breakdowns/collect_tea_breakdowns_split12d.py) simulates and is
NOT run here. Works from any cwd.

The child gets PYTHONDONTWRITEBYTECODE=1 (no __pycache__ under code/) and a
fixed SOURCE_DATE_EPOCH so the PDF's CreationDate, and with it the PDF bytes,
are the same on every run. Neither touches the plotted content or the PNG.
"""
import glob
import os
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
PY = r'C:\Users\saran\anaconda3\envs\IBO_2026\python.exe'
DATA = os.path.join(HERE, 'data')
OUT = os.path.join(HERE, 'output')
JSON = os.path.join(DATA, 'tea_breakdowns_split12d_2026.09.24-17.48.json')
STEM = 'tea_breakdowns_split12d_figure'
# 2026-10-02 12:09 PDT, the manuscript render (BP 126fb46f); PDF metadata only.
SOURCE_DATE_EPOCH = '1790968140'

def main():
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1',
               SOURCE_DATE_EPOCH=SOURCE_DATE_EPOCH)
    script = os.path.join(HERE, 'code', 'plots', 'plot_tea_breakdowns_split12d.py')
    os.makedirs(OUT, exist_ok=True)
    subprocess.run([PY, script, '--data', JSON, '--out-dir', OUT, '--stem', STEM,
                    '--dpi', '300'], env=env, check=True)
    for ext in ('png', 'pdf'):
        newest = max(glob.glob(os.path.join(glob.escape(OUT), f'{STEM}_*.{ext}')),
                     key=os.path.getmtime)
        dst = os.path.join(OUT, f'EDFig5.{ext}')
        shutil.move(newest, dst)   # the timestamped name is render time only
        print('wrote', dst, '<-', os.path.basename(newest))

if __name__ == '__main__':
    main()
