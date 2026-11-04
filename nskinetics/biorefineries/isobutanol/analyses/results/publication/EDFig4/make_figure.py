# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Regenerate the plot panel of Extended Data Figure 4 (PI, ethanol yield and
IRR across the 1-D k_3 slice through the TRY-informed optimum).

Runs the snapshot code/plots/plot_axis_try_informed_sweep.py --param k_3
(sim-safe: reads the saved sweep CSV + anchor JSON only) with the results
root pointed at ./data, then copies the figure into ./output. The manuscript
image is this render after PowerPoint edits (legend removed, rotated
"growth derated at higher k3 values" label added). Works from any cwd.

The plot script writes its figure under its results root, i.e. at
data/publication/Objective-landscape/; this wrapper moves it to ./output and
removes those temporary folders. The child runs with PYTHONDONTWRITEBYTECODE=1
(no __pycache__ under code/) and a fixed SOURCE_DATE_EPOCH, so the PDF's
CreationDate is reproducible (the PNG carries no timestamp).
"""
import os
import sys
import shutil
import subprocess

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
PY = r'C:\Users\saran\anaconda3\envs\IBO_2026\python.exe'
DATA = os.path.join(HERE, 'data')
OUT = os.path.join(HERE, 'output')
STEM = 'k_3_try_informed_sweep_etoh_yield'

def main():
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1',
               # fixed PDF CreationDate: 2026-09-25 00:00 UTC (the day of the
               # original raw-panel render; any constant would do)
               SOURCE_DATE_EPOCH='1790294400')
    script = os.path.join(HERE, 'code', 'plots', 'plot_axis_try_informed_sweep.py')
    subprocess.run([PY, script, '--param', 'k_3', '--dpi', '300'], env=env, check=True)
    src_dir = os.path.join(DATA, 'publication', 'Objective-landscape')
    os.makedirs(OUT, exist_ok=True)
    for ext in ('png', 'pdf'):
        dst = os.path.join(OUT, f'{STEM}.{ext}')
        shutil.move(os.path.join(src_dir, f'{STEM}.{ext}'), dst)
        print('wrote', dst)
    # drop the now-empty output folders the script created under ./data
    for d in (src_dir, os.path.dirname(src_dir)):
        if os.path.isdir(d) and not os.listdir(d):
            os.rmdir(d)

if __name__ == '__main__':
    main()
