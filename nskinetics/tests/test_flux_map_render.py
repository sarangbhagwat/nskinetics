# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Rendering tests for :func:`nskinetics.visualization.draw_flux_map`.

Uses hand-built :class:`~nskinetics.FluxSummary` objects (no simulation, no
shipped kinetic model), so the module stays in the default, non-slow run.
``matplotlib.use('Agg')`` is called at import time; this module is therefore
deliberately not registered in ``nskinetics/tests/__init__.py``, so that a
plain ``import nskinetics`` never changes a user's matplotlib backend.
"""

import os
import matplotlib
matplotlib.use('Agg')
from nskinetics.engine.flux_analysis import FluxSummary
from nskinetics.visualization import FluxMapSpec, draw_flux_map


def _spec():
    return FluxMapSpec(
        nodes={'S': (10., 40., 'S'), 'P': (10., 20., 'P'), 'Q': (40., 20., 'Q')},
        edges={'r1': ('S', 'P'), 'r2': ('P', 'Q')},
        inhibitors={'P': '#D55E00', 'Q': '#0072B2'},
        inhibition_map={'ki': ('r1', 'P'), 'Ke': ('r2', 'Q')},
        enhancement_reactions=set(),
        products=['Q'],
        strip_offsets={}, size_mm=(120., 70.))


def _summary(label, r1flux):
    return FluxSummary(
        label=label, reaction_ids=['r1', 'r2'],
        cumulative_mass={'r1': r1flux, 'r2': 5.0},
        cumulative_flux={'r1': r1flux, 'r2': 5.0},
        fraction_lost={'r1': {'P': 0.3}, 'r2': {'Q': 0.1}},
        fraction_lost_all={'r1': 0.35, 'r2': 0.1},
        final_volume=1.0, final_concentrations={'Q': 42.0}, t_end=40.0,
        inhibitors=['P', 'Q'])


def test_draw_writes_files(tmp_path):
    fig, axes = draw_flux_map([_summary('A', 10.), _summary('B', 6.)], _spec(),
                              save_dir=str(tmp_path), formats=('png', 'pdf'))
    assert len(axes) == 2
    assert os.path.exists(os.path.join(tmp_path, 'flux_map.png'))
    assert os.path.exists(os.path.join(tmp_path, 'flux_map.pdf'))
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_draw_rejects_mismatched_reaction_sets():
    a = _summary('A', 10.)
    b = _summary('B', 6.)
    b.reaction_ids = ['r1']            # mismatch
    try:
        draw_flux_map([a, b], _spec())
        assert False, 'expected ValueError'
    except ValueError:
        pass
