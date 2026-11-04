# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Regenerate Figure 6 (optimization trajectories + final proteome allocation) from the snapshot in code/ and the campaign CSVs in data/.

Runs code/nskinetics/biorefineries/isobutanol/plots/plot_kin_opt_parameter_sets.py
(the repo script at nskinetics commit aca102e == HEAD 15ac19d) with
``--variants irr`` through code/run_plot.py (a MAX_PATH shim for the
relay manifest; see that file), with NSKINETICS_IBO_RESULTS_DIR = data/.
Writes output/kinBO_parameter_sets_irr.{png,pdf}. No __pycache__ is written
(PYTHONDONTWRITEBYTECODE=1 in the child, sys.dont_write_bytecode here and in
the shim); the PDF CreationDate is pinned via SOURCE_DATE_EPOCH.

Sim-safe: the plot script loads paths.py / kinetic_optimization.py /
enzyme_burden.py by file path and reads CSVs, two workbooks and the antimony
text; it never imports nskinetics, biosteam, thermosteam or biorefineries.
Works from any cwd.
"""
import glob
import os
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
PYTHON = sys.executable
CODE = os.path.join(HERE, 'code')
SCRIPT = os.path.join(CODE, 'nskinetics', 'biorefineries', 'isobutanol',
                      'plots', 'plot_kin_opt_parameter_sets.py')
SHIM = os.path.join(CODE, 'run_plot.py')
DATA = os.path.join(HERE, 'data')
OUTPUT = os.path.join(HERE, 'output')
VARIANT = 'irr'
# Fixed PDF CreationDate (matplotlib honours SOURCE_DATE_EPOCH) so the PDF is
# reproducible run to run: 2026-10-17 18:21 -07:00, the timestamp of the
# original kinBO_parameter_sets_irr_2026.10.17-18.21 render. PNG unaffected.
SOURCE_DATE_EPOCH = '1792286460'


def main():
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1',
               SOURCE_DATE_EPOCH=SOURCE_DATE_EPOCH)
    render = os.path.join(OUTPUT, '_render')
    shutil.rmtree(render, ignore_errors=True)
    os.makedirs(render)
    subprocess.run([PYTHON, SHIM, SCRIPT, '--variants', VARIANT,
                    '--out-dir', render],
                   env=env, cwd=os.path.dirname(SCRIPT), check=True)
    for ext in ('png', 'pdf'):
        (src,) = glob.glob(os.path.join(
            render, f'kinBO_parameter_sets_{VARIANT}_*.{ext}'))
        dst = os.path.join(OUTPUT, f'kinBO_parameter_sets_{VARIANT}.{ext}')
        shutil.move(src, dst)
        print('wrote', dst)
    shutil.rmtree(render)


if __name__ == '__main__':
    main()
