# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Regenerate manuscript Figure 2 (biochemically structured kinetic model:
reaction network, controls and strain design variables).

Sim-safe (matplotlib + numpy only): runs the snapshot
code/models/s_cerevisiae_ferm_fb_inhib_mod_ibo/conceptual_diagram.py, calling
draw_conceptual_diagram(show_process_controls=False,
filename='conceptual_diagram_network_only') -- the network-only variant -- and
writes PNG (600 dpi) + PDF into output/. Works from any cwd.
"""
import os
import sys
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
SCRIPT_DIR = os.path.join(HERE, 'code', 'models',
                          's_cerevisiae_ferm_fb_inhib_mod_ibo')
DATA = os.path.join(HERE, 'data')
OUT = os.path.join(HERE, 'output')

SNIPPET = (
    "import importlib.util, sys\n"
    "sys.dont_write_bytecode = True\n"
    "spec = importlib.util.spec_from_file_location('_cd', sys.argv[1])\n"
    "m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)\n"
    "m.draw_conceptual_diagram(save_dir=sys.argv[2], formats=('png', 'pdf'),\n"
    "    dpi=600, show_process_controls=False,\n"
    "    filename='conceptual_diagram_network_only')\n"
)


def main():
    # NSKINETICS_IBO_RESULTS_DIR is set only to keep the per-figure template
    # uniform; conceptual_diagram.py never reads it (data/ stays empty).
    # PYTHONDONTWRITEBYTECODE keeps __pycache__ out of code/.
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1')
    os.makedirs(OUT, exist_ok=True)
    subprocess.run([PY, '-c', SNIPPET,
                    os.path.join(SCRIPT_DIR, 'conceptual_diagram.py'), OUT],
                   env=env, cwd=HERE, check=True)


if __name__ == '__main__':
    main()
