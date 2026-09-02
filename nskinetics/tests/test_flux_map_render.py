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
from nskinetics.visualization.flux_map import (
    C_FLUX, C_JOINT, C_ZERO, LEADER_LW, STRIP_LEN, STRIP_ROW_H,
    STRIP_ROW_GAP, LABEL_OFFSET, _fmt_value, _panel_letter)

C_P = '#D55E00'      # inhibitor colors used by _spec()
C_Q = '#0072B2'


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
        # only the nonzero edge is labeled; letter, header and node labels stay
        assert [t.get_text() for t in ax.texts] == [
            'a', 'A', 'harvest 40 h · Q 42 g/L', '10', 'S', 'P', 'Q']
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


# --- value formatting -------------------------------------------------------

def test_fmt_value_rules():
    assert _fmt_value(0) == '0'          # only an exact zero is bare '0'
    assert _fmt_value(0.0) == '0'
    assert _fmt_value(1e-9) == '<0.1'    # nonzero never prints as '0.0'
    assert _fmt_value(0.049) == '<0.1'
    assert _fmt_value(-1e-3) == '≈0'     # tiny negative reads as ~nothing
    assert _fmt_value(-0.049) == '≈0'
    assert _fmt_value(-0.06) == '-0.1'   # a real negative still prints
    assert _fmt_value(0.05) == '0.1'     # from here it rounds honestly
    assert _fmt_value(0.16) == '0.2'
    assert _fmt_value(1.44) == '1.4'
    assert _fmt_value(9.0) == '9.0'
    assert _fmt_value(12.4) == '12'      # >= 10 -> integer
    assert _fmt_value(236.04) == '236'


def test_edge_labels_use_the_value_formatter():
    spec = _spec()
    summary = _summary('A', 0.24, cumulative_flux={'r1': 0.24, 'r2': 15.6})
    fig, (ax,) = draw_flux_map([summary], spec)
    try:
        texts = [t.get_text() for t in ax.texts]
        assert '0.2' in texts and '16' in texts
    finally:
        plt.close(fig)


# --- panel letters, product labels, header ---------------------------------

def test_panel_letters_are_bold_lowercase_and_ordered():
    summaries = [_summary('A', 10.), _summary('B', 6.)]
    fig, axes = draw_flux_map(summaries, _spec())
    try:
        for ax, letter in zip(axes, 'ab'):
            first = ax.texts[0]
            assert first.get_text() == letter
            assert first.get_fontweight() == 'bold'
            assert first.get_fontsize() == 8.0
            # the scenario label follows on the same baseline
            label = ax.texts[1]
            assert label.get_text() in ('A', 'B')
            assert label.get_position()[1] == first.get_position()[1]
            assert label.get_position()[0] > first.get_position()[0]
    finally:
        plt.close(fig)


def test_product_labels_rename_header_titers():
    spec = _spec(product_labels={'Q': 'quinone'})
    fig, (ax,) = draw_flux_map([_summary('A', 10.)], spec)
    try:
        assert ax.texts[2].get_text() == 'harvest 40 h · quinone 42 g/L'
    finally:
        plt.close(fig)
    # with no product_labels entry the raw species id is used
    fig, (ax,) = draw_flux_map([_summary('A', 10.)], _spec())
    try:
        assert ax.texts[2].get_text() == 'harvest 40 h · Q 42 g/L'
    finally:
        plt.close(fig)


# --- perpendicular value labels --------------------------------------------

def test_value_labels_sit_off_the_line_perpendicular_to_the_edge():
    # r1 S(10,40) -> P(10,20) is vertical: the label goes to the RIGHT.
    # r2 P(10,20) -> Q(40,20) is horizontal: the label goes ABOVE.
    spec = _spec()
    fig, (ax,) = draw_flux_map([_summary('A', 10.)], spec)
    try:
        by_text = {t.get_text(): t for t in ax.texts}
        vert = by_text['10']
        x, y = vert.get_position()
        assert y == 30.0                       # the edge midpoint, unmoved
        assert x == 10.0 + LABEL_OFFSET        # pushed off the line, rightward
        assert vert.get_ha() == 'left' and vert.get_va() == 'center'
        horiz = by_text['5.0']
        x, y = horiz.get_position()
        assert x == 25.0                       # the edge midpoint, unmoved
        assert y == 20.0 + LABEL_OFFSET        # pushed off the line, upward
        assert horiz.get_ha() == 'center' and horiz.get_va() == 'bottom'
    finally:
        plt.close(fig)


def test_value_label_normal_is_flipped_upward_for_a_reversed_edge():
    # Q(40,20) -> P(10,20) runs leftward; the normal must still point up, so
    # the label never lands under the line.
    spec = _spec(edges={'r1': ('S', 'P'), 'r2': ('Q', 'P')})
    fig, (ax,) = draw_flux_map([_summary('A', 10.)], spec)
    try:
        horiz = [t for t in ax.texts if t.get_text() == '5.0'][0]
        assert horiz.get_position()[1] == 20.0 + LABEL_OFFSET
        assert horiz.get_va() == 'bottom'
    finally:
        plt.close(fig)


# --- inactive-branch strips -------------------------------------------------

def test_zero_flux_strip_outlines_are_faint():
    spec = _spec()
    summary = _summary('A', 10., cumulative_flux={'r1': 10., 'r2': 0.0},
                       fraction_lost={'r1': {'P': 0.3}, 'r2': {'Q': 0.0}},
                       fraction_lost_all={'r1': 0.35, 'r2': 0.0})
    fig, (ax,) = draw_flux_map([summary], spec)
    try:
        outlines = [r for r in _strip_rects(ax) if not r.get_fill()]
        faint = [r for r in outlines
                 if tuple(r.get_edgecolor()) == to_rgba(C_ZERO)]
        # r2 carried no flux: all three of its rows (2 inhibitors + joint) are
        # drawn faint; r1's three stay in the normal outline color.
        assert len(outlines) == 6
        assert len(faint) == 3
    finally:
        plt.close(fig)


# --- legend -----------------------------------------------------------------

def test_legend_axes_is_present_and_keyed():
    spec = _spec()
    summaries = [_summary('A', 10.), _summary('B', 6.)]
    fig, axes = draw_flux_map(summaries, spec)
    try:
        # panel axes only in the return value; the legend is the extra axes
        assert len(axes) == 2
        assert len(fig.axes) == 3
        legend = [a for a in fig.axes if a not in axes][0]
        texts = [t.get_text() for t in legend.texts]
        assert 'pathway inactive' in texts
        assert 'all inhibitors' in texts
        for inh in spec.inhibitors:
            assert inh in texts
        assert any('cumulative flux' in t for t in texts)
        assert any('fraction of potential flux lost' in t for t in texts)
        # three flux samples, all at or below the largest drawn flux
        samples = [float(t) for t in texts
                   if t.replace('.', '', 1).isdigit()]
        assert len(samples) == 3
        assert samples == sorted(samples)
        assert max(samples) <= 10.0
        # one line per flux sample plus the dashed "inactive" sample
        assert len(legend.lines) == 4
        assert legend.lines[-1].get_linestyle() not in ('solid', '-')
        # one swatch per inhibitor + the joint one
        assert len(legend.patches) == len(spec.inhibitors) + 1
        # the figure stays within the Nature height budget
        assert fig.get_size_inches()[1] * 25.4 < 170.
    finally:
        plt.close(fig)


def test_legend_gains_an_enhancement_swatch_only_for_a_drawn_reaction():
    fig, axes = draw_flux_map([_summary('A', 10.)],
                              _spec(enhancement_reactions={'r1'}))
    try:
        legend = [a for a in fig.axes if a not in axes][0]
        assert 'enhancement' in [t.get_text() for t in legend.texts]
        assert [p for p in legend.patches if p.get_hatch()]
    finally:
        plt.close(fig)
    # r9 is mapped but not drawn, so no enhancement key is added
    fig, axes = draw_flux_map([_summary('A', 10.)],
                              _spec(enhancement_reactions={'r9'}))
    try:
        legend = [a for a in fig.axes if a not in axes][0]
        assert 'enhancement' not in [t.get_text() for t in legend.texts]
        assert not [p for p in legend.patches if p.get_hatch()]
    finally:
        plt.close(fig)


# --- legend layout ----------------------------------------------------------

def _shipped_like_spec(width_mm):
    """A spec with as many legend items as the shipped one, at ``width_mm``.

    Three inhibitors plus an enhancement reaction that IS drawn, so the legend
    carries its widest possible item set (three flux samples, the inactive
    key, three inhibitor swatches, the joint swatch and the enhancement key).
    """
    return _spec(
        nodes={'S': (10., 50., 'Substrate'), 'P': (10., 25., 'Product'),
               'Q': (40., 25., 'Byproduct')},
        inhibitors={'ethanol': C_P, 'acetate': '#CC79A7',
                    'isobutanol': C_Q},
        enhancement_reactions={'r2'},
        size_mm=(width_mm, 92.))


def _legend_axes(fig, axes):
    return [a for a in fig.axes if a not in axes][0]


def _assert_legend_is_clean(fig, axes):
    """Every legend text inside the legend axes, and no two overlapping."""
    legend = _legend_axes(fig, axes)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    box = legend.get_window_extent(renderer)
    extents = [t.get_window_extent(renderer) for t in legend.texts]
    assert extents, 'the legend drew no text at all'
    for text, e in zip(legend.texts, extents):
        assert box.x0 <= e.x0 and e.x1 <= box.x1, \
            f'{text.get_text()!r} overflows the legend horizontally'
        assert box.y0 <= e.y0 and e.y1 <= box.y1, \
            f'{text.get_text()!r} overflows the legend vertically'
    for i, a in enumerate(extents):
        for b, other in zip(extents[i + 1:], legend.texts[i + 1:]):
            assert not (a.x0 < b.x1 and b.x0 < a.x1
                        and a.y0 < b.y1 and b.y0 < a.y1), (
                f'{legend.texts[i].get_text()!r} overlaps '
                f'{other.get_text()!r}')


def test_legend_is_clean_at_one_panel_width():
    spec = _shipped_like_spec(88.)
    fig, axes = draw_flux_map([_summary('A', 10.)], spec)
    try:
        legend = _legend_axes(fig, axes)
        # the enhancement key is present: this is the crowded case
        assert 'enhancement' in [t.get_text() for t in legend.texts]
        _assert_legend_is_clean(fig, axes)
    finally:
        plt.close(fig)


def test_legend_is_clean_at_two_panel_width():
    spec = _shipped_like_spec(88.)
    fig, axes = draw_flux_map([_summary('A', 10.), _summary('B', 6.)], spec)
    try:
        assert fig.get_size_inches()[0] * 25.4 == pytest.approx(176.)
        _assert_legend_is_clean(fig, axes)
    finally:
        plt.close(fig)


def test_legend_items_wrap_rather_than_collide_in_a_narrow_band():
    # A band far too narrow for one row must wrap, not overprint.
    spec = _shipped_like_spec(60.)
    fig, axes = draw_flux_map([_summary('A', 10.)], spec)
    try:
        legend = _legend_axes(fig, axes)
        ys = {round(t.get_position()[1], 6) for t in legend.texts}
        assert len(ys) > 2                    # items wrapped onto more rows
        _assert_legend_is_clean(fig, axes)
    finally:
        plt.close(fig)


def test_legend_caption_names_the_flux_unit():
    fig, axes = draw_flux_map([_summary('A', 10.)], _spec())
    try:
        texts = [t.get_text() for t in _legend_axes(fig, axes).texts]
        assert ("cumulative flux, g of each step's substrate per L "
                'final broth') in texts
    finally:
        plt.close(fig)


# --- panel letters beyond z -------------------------------------------------

def test_panel_letters_continue_past_z():
    assert [_panel_letter(i) for i in (0, 1, 25)] == ['a', 'b', 'z']
    assert [_panel_letter(i) for i in (26, 27, 51, 52)] == \
        ['aa', 'ab', 'az', 'ba']


# --- enhancement is a property of the value, not of the spec ----------------

def test_negative_fraction_is_hatched_without_an_enhancement_entry():
    # `enhancement_reactions` keys the LEGEND; a negative fraction is drawn as
    # a gain wherever it occurs.
    spec = _spec(enhancement_reactions=set())
    summary = _summary('A', 10.,
                       fraction_lost={'r1': {'P': -0.4}, 'r2': {'Q': 0.1}},
                       fraction_lost_all={'r1': -0.45, 'r2': 0.1})
    fig, (ax,) = draw_flux_map([summary], spec)
    try:
        hatched = [r for r in _strip_rects(ax) if r.get_hatch()]
        assert len(hatched) == 2                       # r1's P row and joint
        assert [t.get_text() for t in ax.texts].count('+') == 2
    finally:
        plt.close(fig)


# --- strip-to-edge leaders --------------------------------------------------

def test_each_offset_strip_gets_a_leader_to_its_edge_midpoint():
    spec = _spec(strip_offsets={'r1': (6., 4.), 'r2': (0., -8.)})
    fig, (ax,) = draw_flux_map([_summary('A', 10.)], spec)
    try:
        leaders = [ln for ln in ax.lines if ln.get_linewidth() == LEADER_LW]
        assert len(leaders) == 2                    # one per inhibited edge
        n_rows = len(spec.inhibitors) + 1
        for rid, (a, b) in spec.edges.items():
            xa, ya, _ = spec.nodes[a]
            xb, yb, _ = spec.nodes[b]
            mx, my = (xa + xb) / 2, (ya + yb) / 2
            hit = [ln for ln in leaders
                   if ln.get_xdata()[-1] == mx and ln.get_ydata()[-1] == my]
            assert len(hit) == 1, f'no leader points at {rid}'
            # ... and it starts on a corner of that reaction's strip
            dx, dy = spec.strip_offsets[rid]
            sx, sy = mx + dx, my + dy
            y_bot = sy - (n_rows - 1) * (STRIP_ROW_H + STRIP_ROW_GAP)
            corners = {(sx, y_bot), (sx, sy + STRIP_ROW_H),
                       (sx + STRIP_LEN, y_bot),
                       (sx + STRIP_LEN, sy + STRIP_ROW_H)}
            start = (hit[0].get_xdata()[0], hit[0].get_ydata()[0])
            assert start in corners
            assert to_rgba(hit[0].get_color()) == to_rgba(C_ZERO)
    finally:
        plt.close(fig)


def test_no_leader_when_the_strip_covers_its_own_edge_midpoint():
    spec = _spec(strip_offsets={'r1': (-STRIP_LEN / 2, 0.),
                                'r2': (-STRIP_LEN / 2, 0.)})
    fig, (ax,) = draw_flux_map([_summary('A', 10.)], spec)
    try:
        assert not [ln for ln in ax.lines
                    if ln.get_linewidth() == LEADER_LW]
    finally:
        plt.close(fig)


# --- spec.reactions ---------------------------------------------------------

def test_spec_reactions_appends_mapped_but_undrawn_reactions():
    spec = _spec(inhibition_map={'ki': ('r1', 'P'), 'Ke': ('r2', 'Q'),
                                 'kd': ('r9', 'P'), 'kd2': ('r9', 'Q')})
    # drawn edges first, in edge order; then undrawn mapped ones, first-seen
    assert spec.reactions == ['r1', 'r2', 'r9']
    # a spec whose map names only drawn reactions is just the edge list
    assert _spec().reactions == list(_spec().edges)
