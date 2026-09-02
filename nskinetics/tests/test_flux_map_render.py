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
Assertions inspect the artists actually added to each axes, so a panel that
silently drew nothing would fail rather than pass.

``matplotlib.use('Agg')`` is called at import time; this module is therefore
deliberately not registered in ``nskinetics/tests/__init__.py``, so that a
plain ``import nskinetics`` never changes a user's matplotlib backend.
"""

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pytest
from matplotlib.colors import to_rgba
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

from nskinetics.engine.flux_analysis import FluxSummary
from nskinetics.visualization import FluxMapSpec, draw_flux_map
from nskinetics.visualization.flux_map import C_FLUX, C_JOINT

C_ZERO = '#B4B2A9'   # color of a zero-flux (dashed) edge
C_P = '#D55E00'      # inhibitor colors used by _spec()
C_Q = '#0072B2'
STRIP_LEN = 9.0      # full strip-row length in mm


def _spec(**kwargs):
    """The synthetic two-edge layout; ``kwargs`` override any field."""
    fields = dict(
        nodes={'S': (10., 40., 'S'), 'P': (10., 20., 'P'), 'Q': (40., 20., 'Q')},
        edges={'r1': ('S', 'P'), 'r2': ('P', 'Q')},
        inhibitors={'P': C_P, 'Q': C_Q},
        inhibition_map={'ki': ('r1', 'P'), 'Ke': ('r2', 'Q')},
        enhancement_reactions=set(),
        products=['Q'],
        strip_offsets={}, size_mm=(120., 70.))
    fields.update(kwargs)
    return FluxMapSpec(**fields)


def _summary(label, r1flux, **kwargs):
    """A synthetic summary; ``kwargs`` override any field."""
    fields = dict(
        label=label, reaction_ids=['r1', 'r2'],
        cumulative_mass={'r1': r1flux, 'r2': 5.0},
        cumulative_flux={'r1': r1flux, 'r2': 5.0},
        fraction_lost={'r1': {'P': 0.3}, 'r2': {'Q': 0.1}},
        fraction_lost_all={'r1': 0.35, 'r2': 0.1},
        final_volume=1.0, final_concentrations={'Q': 42.0}, t_end=40.0,
        inhibitors=['P', 'Q'])
    fields.update(kwargs)
    return FluxSummary(**fields)


def _arrows(ax):
    return [p for p in ax.patches if isinstance(p, FancyArrowPatch)]


def _nodes(ax):
    return [p for p in ax.patches if isinstance(p, FancyBboxPatch)]


def _strip_rects(ax):
    return [p for p in ax.patches if isinstance(p, Rectangle)]


def _fills(ax):
    """Strip rectangles that are actually filled (alpha > 0)."""
    return [r for r in _strip_rects(ax) if r.get_facecolor()[3] > 0]


def _expected_strip_rects(spec, summary):
    """Strip rectangles the spec/summary pair should produce.

    One outline per inhibitor row plus one proportional fill per nonzero row,
    then the same pair for the grey joint row, for every reaction that appears
    in ``inhibition_map``.
    """
    mapped = {rid for rid, _inh in spec.inhibition_map.values()}
    n = 0
    for rid in spec.edges:
        if rid not in mapped:
            continue
        fl = summary.fraction_lost.get(rid, {})
        for inh in spec.inhibitors:
            n += 1                                        # row outline
            if fl.get(inh, 0.0):
                n += 1                                    # proportional fill
        n += 1                                            # joint outline
        if summary.fraction_lost_all.get(rid, 0.0):
            n += 1                                        # joint fill
    return n


def test_draw_writes_files(tmp_path):
    spec = _spec()
    summaries = [_summary('A', 10.), _summary('B', 6.)]
    fig, axes = draw_flux_map(summaries, spec, save_dir=str(tmp_path),
                              formats=('png', 'pdf'))
    try:
        assert len(axes) == 2
        assert os.path.exists(os.path.join(tmp_path, 'flux_map.png'))
        assert os.path.exists(os.path.join(tmp_path, 'flux_map.pdf'))
        for ax, summary in zip(axes, summaries):
            n_strip = _expected_strip_rects(spec, summary)
            # 2 reactions x (2 row outlines + 1 nonzero fill
            #                + joint outline + joint fill)
            assert n_strip == 10
            assert len(_arrows(ax)) == len(spec.edges) == 2
            assert len(_strip_rects(ax)) == n_strip
            assert len(_nodes(ax)) == len(spec.nodes) == 3
            # nothing counted twice and nothing else on the panel
            assert len(ax.patches) == 2 + n_strip + 3
    finally:
        plt.close(fig)


def test_draw_rejects_mismatched_reaction_sets():
    a = _summary('A', 10.)
    b = _summary('B', 6.)
    b.reaction_ids = ['r1']            # mismatch
    with pytest.raises(ValueError):
        draw_flux_map([a, b], _spec())


def test_empty_summaries_raises():
    with pytest.raises(ValueError):
        draw_flux_map([], _spec())


def test_single_summary_returns_a_one_element_axes_list():
    fig, axes = draw_flux_map([_summary('A', 10.)], _spec())
    try:
        assert isinstance(axes, list)
        assert len(axes) == 1
    finally:
        plt.close(fig)


def test_zero_flux_edge_is_dashed_grey_and_unlabeled():
    spec = _spec()
    summary = _summary('A', 10., cumulative_flux={'r1': 10., 'r2': 0.0})
    fig, (ax,) = draw_flux_map([summary], spec)
    try:
        arrows = _arrows(ax)
        assert len(arrows) == 2
        dashed = [a for a in arrows
                  if tuple(a.get_edgecolor()) == to_rgba(C_ZERO)]
        solid = [a for a in arrows
                 if tuple(a.get_edgecolor()) == to_rgba(C_FLUX)]
        assert len(dashed) == 1 and len(solid) == 1
        assert dashed[0].get_linestyle() not in ('solid', '-')
        # only the nonzero edge is labeled; header and node labels remain
        assert [t.get_text() for t in ax.texts] == [
            'A', 'harvest 40 h; Q 42 g/L', '10', 'S', 'P', 'Q']
    finally:
        plt.close(fig)


def test_enhancement_rows_are_hatched_and_marked():
    spec = _spec(enhancement_reactions={'r1'})
    summary = _summary('A', 10.,
                       fraction_lost={'r1': {'P': -0.4}, 'r2': {'Q': 0.1}},
                       fraction_lost_all={'r1': -0.45, 'r2': 0.1})
    fig, (ax,) = draw_flux_map([summary], spec)
    try:
        hatched = [r for r in _strip_rects(ax) if r.get_hatch()]
        # r1's P row and r1's grey joint row; r2 is not an enhancement
        assert len(hatched) == 2
        assert {tuple(r.get_facecolor()) for r in hatched} == {
            to_rgba(C_P), to_rgba(C_JOINT)}
        assert [t.get_text() for t in ax.texts].count('+') == 2
    finally:
        plt.close(fig)


def test_fractions_above_one_clamp_to_the_full_row():
    spec = _spec()
    summary = _summary('A', 10.,
                       fraction_lost={'r1': {'P': 1.8}, 'r2': {'Q': 0.1}},
                       fraction_lost_all={'r1': 1.9, 'r2': 0.1})
    fig, (ax,) = draw_flux_map([summary], spec)
    try:
        widths = sorted(round(r.get_width(), 6) for r in _fills(ax))
        # r2 keeps its 0.1 rows; both of r1's over-unity rows stop at the row
        assert widths == [0.9, 0.9, STRIP_LEN, STRIP_LEN]
        assert max(widths) == STRIP_LEN
    finally:
        plt.close(fig)


def test_draw_does_not_mutate_global_rcparams(tmp_path):
    def snapshot():
        return (plt.rcParams['pdf.fonttype'],
                plt.rcParams['ps.fonttype'],
                plt.rcParams['svg.fonttype'],
                list(plt.rcParams['font.family']),
                list(plt.rcParams['font.sans-serif']),
                plt.rcParams['text.color'])

    before = snapshot()
    fig, _axes = draw_flux_map([_summary('A', 10.)], _spec(),
                               save_dir=str(tmp_path), formats=('png', 'pdf'))
    plt.close(fig)
    assert snapshot() == before
