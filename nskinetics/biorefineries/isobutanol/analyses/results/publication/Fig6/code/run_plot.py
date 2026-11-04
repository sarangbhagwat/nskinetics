# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""MAX_PATH shim for running the snapshot of plot_kin_opt_parameter_sets.py
from this figure folder.

The relay campaign's preload manifest is looked up by the plot script as
``<results root>/<relay study>_relay_manifest.csv``. With the results root at
``<figure folder>/data`` that path is 260-262 characters, over Windows'
MAX_PATH (260 incl. NUL; LongPathsEnabled = 0), so the file cannot be stored
under its real name here. It is stored as ``data/rl15c111dc_relay_manifest.csv``
(byte-identical copy) and this shim redirects the one long path to it, by
wrapping ``os.path.isfile`` and ``pandas.read_csv`` in this process only. The
snapshot plot script itself is unmodified.

Usage: python run_plot.py <plot script> [plot-script args...]
Sim-safe: imports only stdlib + pandas, then runs the plot script as __main__.
"""
import os
import runpy
import sys

sys.dont_write_bytecode = True  # keep the figure folder free of __pycache__

import pandas as pd

RELAY_STUDY = ('kin_opt_ethanol_isobutanol_metabolic_split_12d_pi_log-tail_gp'
               '_rb0.001-4_ib0.75-1.5_aA_rl15c111dc_burden')
SHORT_MANIFEST = 'rl15c111dc_relay_manifest.csv'


def _key(path):
    return os.path.normcase(os.path.abspath(str(path)))


def main():
    root = os.path.abspath(os.environ['NSKINETICS_IBO_RESULTS_DIR'])
    long_path = os.path.join(root, f'{RELAY_STUDY}_relay_manifest.csv')
    short_path = os.path.join(root, SHORT_MANIFEST)
    redirect = {_key(long_path): short_path}

    def remap(path):
        if isinstance(path, (str, os.PathLike)):
            return redirect.get(_key(path), path)
        return path

    _isfile, _read_csv = os.path.isfile, pd.read_csv
    os.path.isfile = lambda p: _isfile(remap(p))
    pd.read_csv = lambda p, *a, **k: _read_csv(remap(p), *a, **k)

    script = os.path.abspath(sys.argv[1])
    sys.argv = [script] + sys.argv[2:]
    runpy.run_path(script, run_name='__main__')


if __name__ == '__main__':
    main()
