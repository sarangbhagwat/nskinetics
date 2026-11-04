# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Regenerate manuscript Figure 3 (simplified block flow diagram of the
biorefinery, panels A and B plus legend).

Sim-safe (matplotlib only, no data inputs): runs the snapshot
code/plots/plot_block_flow_diagram.py with --out-dir output/, which writes
block_flow_diagram.{pdf,svg,png} (PNG at 600 dpi). Works from any cwd.
"""
import os
import sys
import subprocess

# Keep the snapshot folder free of __pycache__ (also for the child below).
sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
SCRIPT = os.path.join(HERE, 'code', 'plots', 'plot_block_flow_diagram.py')
DATA = os.path.join(HERE, 'data')
OUT = os.path.join(HERE, 'output')


def main():
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1')
    os.makedirs(OUT, exist_ok=True)
    subprocess.run([PY, SCRIPT, '--out-dir', OUT], env=env, cwd=HERE,
                   check=True)


if __name__ == '__main__':
    main()
