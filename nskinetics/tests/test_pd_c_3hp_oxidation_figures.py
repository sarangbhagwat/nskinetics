# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Rendering tests for ``nskinetics.models.pd_c_3hp_oxidation.figures``
(Fig. 9, SI Fig. S12 and Fig. 11 of Kim et al., Appl. Catal. B: Environ.
Energy 386 (2026) 126403, drawn from the NSKinetics models, and the
verification figures drawn next to them).

The figures are built in memory and their artists inspected, so a panel that
silently drew nothing would fail; one small PNG is written to ``tmp_path``.
The models are small, so the module stays in the default, non-slow run.

``matplotlib.use('Agg')`` is called at import time; this module is therefore
deliberately not registered in ``nskinetics/tests/__init__.py``.
"""
import os
import struct

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.colors import to_hex
from matplotlib.contour import ContourSet
from matplotlib.lines import Line2D
from matplotlib.text import Text

import nskinetics.models.pd_c_3hp_oxidation as pkg
from nskinetics.models.pd_c_3hp_oxidation import figures


def _lines_with_data(ax):
    return [line for line in ax.get_lines() if len(line.get_xdata()) > 1]


def _species_key(ax):
    """The two-column species key of a profile panel."""
    (key,) = [artist for artist in ax.artists
              if artist.get_gid() == 'species_key']
    return key


def _key_labels(key):
    """Text artists of a species key: the 'mod' column, then the 'exp'
    column, five rows each."""
    labels = [text for text in key.findobj(Text) if text.get_text()]
    assert [text.get_text() for text in labels] == (
        ['%s$^{mod}$' % name for name in ('3HP', 'MA', 'AAc', 'OXA',
                                          'Others')]
        + ['%s$^{exp}$' % name for name in ('3HP', 'MA', 'AAc', 'OXA',
                                            'Others')])
    return labels[:5], labels[5:]


def _marker_only_lines(ax):
    """The series drawn as markers without a line (experimental points)."""
    return [line for line in ax.get_lines()
            if line.get_linestyle() == 'None' and len(line.get_xdata()) > 1]


def _model_lines(ax, n_points):
    """The solid model lines of a profile panel, in the order of
    ``NORMALIZED``."""
    return [line for line in ax.get_lines()
            if line.get_linestyle() == '-' and line.get_marker() == 'None'
            and len(line.get_xdata()) == n_points]


def _assert_lines_are_the_model(lines, km, selections):
    """Each line holds the trajectory of ``km`` for its selection."""
    assert len(lines) == len(selections)
    times = lines[0].get_xdata()
    sim = pkg.simulate_at_times(km, times, list(selections))
    for line, name in zip(lines, selections):
        assert np.array_equal(line.get_xdata(), times)
        assert np.allclose(line.get_ydata(), sim[name], rtol=0., atol=1e-12)


def _assert_species_series(ax, sim, experiment):
    """Each species is drawn with its own data in its own colour.

    For every name of ``NORMALIZED``: the one solid line of colour
    ``figures._COLORS[name]`` holds ``sim[name]``, and the one series of
    circles with that edge colour holds ``experiment[name]`` at the
    samples on the time axis.
    """
    drawn = experiment[experiment['time'] <= figures._T_END]
    for name in pkg.NORMALIZED:
        color = figures._COLORS[name]
        (model,) = [line for line in ax.get_lines()
                    if line.get_linestyle() == '-'
                    and line.get_marker() == 'None'
                    and to_hex(line.get_color()) == color]
        assert np.array_equal(model.get_xdata(), sim['time']), name
        assert np.allclose(model.get_ydata(), sim[name], rtol=0.,
                           atol=1e-12), name
        (circles,) = [line for line in _marker_only_lines(ax)
                      if line.get_marker() == 'o'
                      and to_hex(line.get_markeredgecolor()) == color]
        assert circles.get_markerfacecolor() == 'none', name
        assert np.array_equal(circles.get_xdata(), drawn['time']), name
        assert np.array_equal(circles.get_ydata(), drawn[name]), name


def _texts(ax):
    return [text.get_text() for text in ax.texts]


# Oxidant tags: text, box face and box edge (the face colours are those of
# the published Fig. 9).
_TAG_BOXES = {'H2O2': ('H$_2$O$_2$', '#fbb4b4', '#e98686'),
              'O2': ('O$_2$', '#c3c3dc', '#8f8fb8')}


def _assert_tag(ax, oxidant):
    """The oxidant tag of ``ax`` has the text and box colours of
    ``oxidant``."""
    text, face, edge = _TAG_BOXES[oxidant]
    (tag,) = [t for t in ax.texts if t.get_text() == text]
    box = tag.get_bbox_patch()
    assert to_hex(box.get_facecolor()) == face, oxidant
    assert to_hex(box.get_edgecolor()) == edge, oxidant


def _png_size(path):
    with open(path, 'rb') as f:
        header = f.read(24)
    return struct.unpack('>II', header[16:24])


def test_series_colours_and_labels_are_those_of_the_published_figures():
    # The colours of the species series as published (Fig. 9, SI Fig. S12;
    # the solver-emulation figure uses the same), one per species, and the
    # contour colours of Fig. 11. The figure tests below look the series
    # up by these colours.
    assert figures._COLORS == {
        'X_HP3': '#4d4d4d', 'X_MA': '#a66bd6', 'X_AAc': '#1f6fdb',
        'X_OXA': '#f26a1b', 'X_D': '#5b9bd5'}
    assert figures._LABELS == {'X_HP3': '3HP', 'X_MA': 'MA', 'X_AAc': 'AAc',
                               'X_OXA': 'OXA', 'X_D': 'Others'}
    assert (figures._C_H2O2, figures._C_OPA) == ('#f28e8e', '#2aa35a')
    assert figures._CONTOUR_COLORS == {0.1: (0.30, 0.50, 0.82),
                                       0.3: (0.35, 0.57, 0.40),
                                       0.5: (0.82, 0.56, 0.09)}
    assert (figures._C_RIDGE, figures._C_RED_STAR,
            figures._C_PURPLE_STAR) == ('#e41a1c', '#f0525a', '#7b6fc4')


# Reaction lists of the published Fig. 9, (reactants, products) per
# reaction, as drawn.
FIG9_REACTIONS = {
    'H2O2': [('3-HP + H$_2$O$_2$', 'OPA + 2H$_2$O'),
             ('OPA + H$_2$O$_2$', 'MA + H$_2$O'),
             ('MA', 'AAc + CO$_2$'),
             ('MA + H$_2$O$_2$', 'OXA + 2H$_2$ + CO$_2$'),
             ('3-HP', 'Others'), ('MA', 'Others'),
             ('H$_2$O$_2$', 'H$_2$O + 0.5 O$_2$')],
    'O2': [('3-HP + 0.5 O$_2$', 'OPA + H$_2$O'),
           ('OPA + 0.5 O$_2$', 'MA'),
           ('MA', 'AAc + CO$_2$'),
           ('MA + O$_2$', 'OXA + H$_2$ + CO$_2$'),
           ('3-HP', 'Others'), ('MA', 'Others')]}


def test_fig9_draws_model_experiment_and_deposited_curve():
    fig = figures.plot_fig9()
    try:
        panel_a, reactions_a, panel_b, reactions_b = fig.axes
        for panel, letter, tag, oxidant, read, load in (
                (panel_a, 'a', 'H$_2$O$_2$', 'H2O2', pkg.read_fig9a_h2o2,
                 pkg.load_h2o2_model),
                (panel_b, 'b', 'O$_2$', 'O2', pkg.read_fig9b_o2,
                 pkg.load_o2_model)):
            # 5 model lines + 5 experimental series + 5 deposited polylines
            assert len(_lines_with_data(panel)) == 15
            # the model lines are the Table 3 model of the panel's oxidant
            _assert_lines_are_the_model(
                _model_lines(panel, len(figures._T_DENSE)), load(),
                pkg.NORMALIZED)
            # each species: its model line and its experimental circles, in
            # its colour. Markers are drawn whole on the frame, so no sample
            # taken after the end of the time axis may be drawn (Fig. 9b has
            # samples at 1470 and 1680 min).
            _assert_species_series(
                panel, pkg.simulate_at_times(load(), figures._T_DENSE),
                read()['experiment'])
            assert len(_marker_only_lines(panel)) == 5
            # the five deposited polylines, in the order of NORMALIZED,
            # black and dotted, with their key
            deposited = read()['model']
            drawn = [line for line in panel.get_lines()
                     if len(line.get_xdata()) == len(deposited)]
            assert len(drawn) == 5
            for line, name in zip(drawn, pkg.NORMALIZED):
                assert np.array_equal(line.get_xdata(), deposited['time'])
                assert np.array_equal(line.get_ydata(), deposited[name]), name
                assert to_hex(line.get_color()) == '#000000'
                # matplotlib has no public getter of the dash sequence;
                # this attribute holds it unscaled (matplotlib 3.10, 3.11)
                assert line.is_dashed()
                assert line._unscaled_dash_pattern == (0, (1, 2))
            legend = panel.get_legend()
            assert [text.get_text() for text in legend.get_texts()] == [
                'Deposited model curve']
            (handle,) = legend.legend_handles
            assert to_hex(handle.get_color()) == '#000000'
            assert handle.is_dashed()
            assert handle._unscaled_dash_pattern == (0, (1, 2))
            # species key: a line ('mod') and a circle ('exp') per species,
            # in the colour of its series
            key = _species_key(panel)
            _key_labels(key)
            handles = key.findobj(Line2D)
            colors = [figures._COLORS[name] for name in pkg.NORMALIZED]
            assert [to_hex(line.get_color()) for line in handles] == 2*colors
            assert [(line.get_linestyle(), line.get_marker())
                    for line in handles] == (5*[('-', 'None')]
                                             + 5*[('None', 'o')])
            # oxidant tag and panel letter
            assert _texts(panel) == [tag, letter]
            _assert_tag(panel, oxidant)
            assert panel.get_xlim() == (-30., 1440.)
            assert panel.get_ylim() == (-0.05, 1.)
            assert panel.get_xlabel() == 'Time (min)'
            assert panel.get_ylabel() == 'Normalized\nConcentration'
            # inward ticks with minor ticks half-way, as published
            assert list(panel.get_xticks()) == [0, 500, 1000]
            assert 250. in panel.xaxis.get_minorticklocs()
            assert panel.xaxis.get_tick_params()['direction'] == 'in'
            assert panel.yaxis.get_tick_params(which='minor')[
                'direction'] == 'in'
        # five pieces of text per reaction: "Rxn n" (navy, as published),
        # reactants, arrow, k over the arrow (red), products
        for reactions, oxidant in ((reactions_a, 'H2O2'),
                                   (reactions_b, 'O2')):
            texts = reactions.texts
            n = len(FIG9_REACTIONS[oxidant])
            assert len(texts) == 5*n
            assert [text.get_text() for text in texts[::5]] == [
                'Rxn %d' % i for i in range(1, n + 1)]
            assert [text.get_text() for text in texts[3::5]] == [
                'k$_%d$' % i for i in range(1, n + 1)]
            assert list(zip([text.get_text() for text in texts[1::5]],
                            [text.get_text() for text in texts[4::5]])) == (
                FIG9_REACTIONS[oxidant])
            assert {to_hex(text.get_color()) for text in texts[::5]} == {
                '#3653a0'}
            assert {to_hex(text.get_color()) for text in texts[3::5]} == {
                '#ee2c32'}
    finally:
        plt.close(fig)


@pytest.mark.parametrize('font', ['DejaVu Sans', None])
def test_species_key_rows_are_aligned_and_inside_the_panel(font,
                                                           monkeypatch):
    # With the last fallback font (DejaVu Sans, the widest) as well as with
    # the default font list: the rows of the two key columns share their
    # baselines, and the key stays inside the axes frame. The figures are
    # drawn inside the rc context of the module, as save_figures draws
    # them: the font list is looked up when a text is drawn.
    if font is not None:
        monkeypatch.setitem(figures._RC_PARAMS, 'font.sans-serif', [font])
    with plt.rc_context(figures._RC_PARAMS):
        fig = figures.plot_fig9()
        try:
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            panels = [ax for ax in fig.axes
                      if ax.get_xlabel() == 'Time (min)']
            assert len(panels) == 2
            for panel in panels:
                key = _species_key(panel)
                frame = panel.get_window_extent(renderer)
                box = key.get_window_extent(renderer)
                assert frame.x0 < box.x0 and box.x1 < frame.x1 - 2.
                assert frame.y0 < box.y0 and box.y1 < frame.y1
                model, experiment = _key_labels(key)
                rows = []
                for left, right in zip(model, experiment):
                    y_left = left.get_window_extent(renderer).y0
                    y_right = right.get_window_extent(renderer).y0
                    assert abs(y_left - y_right) < 0.5       # pixels
                    rows.append(y_left)
                # equal row pitch
                assert np.ptp(np.diff(rows)) < 0.5
        finally:
            plt.close(fig)


def test_figS12_draws_the_insets(tmp_path):
    fig = figures.plot_figS12()
    try:
        panel_a, panel_b = fig.axes[:2]
        (inset_a,), (inset_b,) = panel_a.child_axes, panel_b.child_axes
        assert inset_a.get_ylim() == (-0.3, 15.)
        assert inset_b.get_xlim() == (-50., 600.)
        n = len(figures._T_DENSE)
        # species lines of the main axes and of inset a: the Table 3 models
        for ax, load in ((panel_a, pkg.load_h2o2_model),
                         (inset_a, pkg.load_h2o2_model),
                         (panel_b, pkg.load_o2_model)):
            _assert_lines_are_the_model(_model_lines(ax, n), load(),
                                        pkg.NORMALIZED)
        # each species: its model line and its experimental circles, in its
        # colour, on the main axes and in inset a (inset b holds no species)
        for ax, load, read in ((panel_a, pkg.load_h2o2_model,
                                pkg.read_fig9a_h2o2),
                               (inset_a, pkg.load_h2o2_model,
                                pkg.read_fig9a_h2o2),
                               (panel_b, pkg.load_o2_model,
                                pkg.read_fig9b_o2)):
            _assert_species_series(
                ax, pkg.simulate_at_times(load(), figures._T_DENSE),
                read()['experiment'])
        assert not _marker_only_lines(inset_b)
        # the dashed lines are the predicted H2O2 and OPA of those models, in
        # the colours of the published panels
        for ax, load, name, scale, color in (
                (panel_a, pkg.load_h2o2_model, 'eq_H2O2', 1., '#f28e8e'),
                (inset_a, pkg.load_h2o2_model, 'eq_H2O2', 1., '#f28e8e'),
                (panel_b, pkg.load_o2_model, 'x_OPA', 1., '#2aa35a'),
                (inset_b, pkg.load_o2_model, 'x_OPA', 1e5, '#2aa35a')):
            (dashed,) = [line for line in ax.get_lines()
                         if line.get_linestyle() == '--']
            sim = pkg.simulate_at_times(load(), dashed.get_xdata(), [name])
            assert np.allclose(dashed.get_ydata(), scale*sim[name],
                               rtol=0., atol=1e-12*scale)
            assert to_hex(dashed.get_color()) == color
        # panel letters and oxidant tags
        assert (_texts(panel_a), _texts(panel_b)) == (['a'], ['b'])
        assert _texts(inset_a) == ['H$_2$O$_2$']
        assert _texts(inset_b) == [r'$\times10^{-5}$', 'O$_2$']
        _assert_tag(inset_a, 'H2O2')
        _assert_tag(inset_b, 'O2')
        # key of panel b: line and circle of each species in its colour,
        # labelled as published (D for "Others")
        legend = panel_b.get_legend()
        assert [text.get_text() for text in legend.get_texts()] == [
            '3HP', 'MA', 'AAc', 'OXA', 'D']
        assert [to_hex(handle.get_color())
                for handle in legend.legend_handles] == [
            figures._COLORS[name] for name in pkg.NORMALIZED]
        # keys of the insets: each label on the handle of its series
        for ax, keys in (
                (inset_a, [('H$_2$O$_2$$^{predicted}$', '#f28e8e', '--',
                            'None'),
                           ('Published (digitized)', '#1a1a1a', '-', 'None'),
                           ('Deposited trace', '#1a1a1a', 'None', 's')]),
                (inset_b, [('OPA$^{predicted}$', '#2aa35a', '--', 'None'),
                           ('Published (digitized)', '#1a1a1a', '-',
                            'D')])):
            legend = ax.get_legend()
            assert [(text.get_text(), to_hex(handle.get_color()),
                     handle.get_linestyle(), handle.get_marker())
                    for text, handle in zip(legend.get_texts(),
                                            legend.legend_handles)] == keys
        h2o2 = inset_a.get_lines()[-1]
        assert h2o2.get_ydata()[0] == pytest.approx(15.)
        # (a) the 'Published (digitized)' line is the digitized inset of SI
        # Fig. S12a, and the squares are the deposited trace
        published = pkg.read_figS12a_h2o2_inset_digitized()
        trace = pkg.read_fig9a_h2o2()['h2o2_trace']
        (line,) = [line for line in inset_a.get_lines()
                   if line.get_linestyle() == '-'
                   and to_hex(line.get_color()) == figures._C_INK]
        assert np.array_equal(line.get_xdata(), published['time'])
        assert np.array_equal(line.get_ydata(), published['eq_H2O2'])
        (squares,) = [line for line in inset_a.get_lines()
                      if line.get_marker() == 's']
        assert np.array_equal(squares.get_xdata(), trace['time'])
        assert np.array_equal(squares.get_ydata(), trace['eq_H2O2'])
        # (b) published polyline and predicted OPA, on the 1e-5 axis
        published, opa = inset_b.get_lines()
        digitized = pkg.read_figS12b_opa_inset_digitized()
        assert np.array_equal(published.get_xdata(), digitized['time'])
        assert np.array_equal(published.get_ydata(), digitized['x_OPA']*1e5)
        assert len(published.get_xdata()) == 13
        assert 1.40 <= max(opa.get_ydata()) <= 1.47
        path = os.path.join(str(tmp_path), 'figS12.png')
        fig.savefig(path, dpi=50)
        assert os.path.getsize(path) > 0
    finally:
        plt.close(fig)


def test_figures_do_not_mutate_global_rcparams():
    def snapshot():
        return (plt.rcParams['pdf.fonttype'],
                plt.rcParams['font.size'],
                list(plt.rcParams['font.family']),
                list(plt.rcParams['font.sans-serif']),
                plt.rcParams['axes.linewidth'])

    before = snapshot()
    plt.close(figures.plot_figS12())
    plt.close(figures.plot_fig11())
    assert snapshot() == before


# --- Fig. 11 (heat maps) -----------------------------------------------------

def _heat_map_axes(fig):
    """Main axes of the heat-map panels (those that hold an inset or a mesh
    and are not themselves insets or color bars)."""
    return [ax for ax in fig.axes if ax.get_xlabel() == 'Time (min)']


def _mesh(ax):
    (mesh,) = [c for c in ax.collections if c.get_array() is not None
               and np.size(c.get_array()) == 71*300]
    return np.asarray(mesh.get_array()).reshape(71, 300)


def _star_positions(ax):
    return [(line.get_xdata()[0], line.get_ydata()[0])
            for line in ax.get_lines()
            if line.get_marker() == '*' and line.get_markerfacecolor() != 'none']


def _assert_heat_map_lines(ax, T_C, times, matrix):
    """Contour lines at 0.1, 0.3 and 0.5 in their colours, the red dashed
    line through the row maxima, the red star on the largest cell and the
    purple star on the maximum of the 50 C row."""
    (contours,) = [c for c in ax.collections if isinstance(c, ContourSet)]
    assert list(contours.levels) == [0.1, 0.3, 0.5]
    assert [to_hex(color) for color in contours.get_edgecolor()] == [
        to_hex(figures._CONTOUR_COLORS[level]) for level in (0.1, 0.3, 0.5)]
    (ridge,) = [line for line in ax.get_lines()
                if len(line.get_xdata()) == len(T_C)]
    assert np.array_equal(ridge.get_xdata(), pkg.heat_map_ridge(times, matrix))
    assert np.array_equal(ridge.get_ydata(), T_C)
    assert (to_hex(ridge.get_color()), ridge.get_linestyle()) == (
        '#e41a1c', '--')
    stars = {to_hex(line.get_markerfacecolor()): (line.get_xdata()[0],
                                                  line.get_ydata()[0])
             for line in ax.get_lines() if line.get_marker() == '*'
             and line.get_markerfacecolor() != 'none'}
    _, T_max, t_max = pkg.heat_map_maximum(T_C, times, matrix)
    row = list(T_C).index(50.)
    assert stars == {'#f0525a': (t_max, T_max),
                     '#7b6fc4': (times[np.argmax(matrix[row])], 50.)}


def test_fig11_draws_both_heat_maps_with_the_default_constant_set():
    fig = figures.plot_fig11()
    try:
        panel_a, panel_b = _heat_map_axes(fig)
        T_C, times, selectivity, MA_yield = pkg.fig11_heat_maps()
        # contour labels of the main panels as published: (level, side of
        # the red dashed line) and how many labels it carries
        published_labels = {
            'a': {(0.1, True): 1, (0.3, True): 1, (0.5, True): 1,
                  (0.1, False): 1, (0.3, False): 1, (0.5, False): 1},
            'b': {(0.1, True): 1, (0.3, True): 1, (0.5, True): 2,
                  (0.1, False): 2, (0.3, False): 2, (0.5, False): 2}}
        for panel, letter, matrix, label, value in (
                (panel_a, 'a', selectivity, 'MA Selectivity', 0.5931),
                (panel_b, 'b', MA_yield, 'MA Yield', 0.5484)):
            assert np.array_equal(_mesh(panel), matrix)
            _assert_heat_map_lines(panel, T_C, times, matrix)
            assert _texts(panel)[-1] == letter
            ridge_times = pkg.heat_map_ridge(times, matrix)
            labels = {}
            for text in panel.texts:
                if text.get_text() in ('0.1', '0.3', '0.5'):
                    t_label, T_label = text.get_position()
                    key = (float(text.get_text()),
                           t_label < np.interp(T_label, T_C, ridge_times))
                    labels[key] = labels.get(key, 0) + 1
            assert labels == published_labels[letter]
            assert panel.get_xlim() == (7.2, 1682.8)
            assert panel.get_ylim() == (29.5, 100.5)
            # inward ticks on all four sides, as published
            for axis, far_side in ((panel.xaxis, 'top'),
                                   (panel.yaxis, 'right')):
                ticks = axis.get_tick_params()
                assert ticks['direction'] == 'in' and ticks[far_side]
            # red star on the largest cell, purple star on the 50 C row maximum
            red, purple = _star_positions(panel)
            maximum = pkg.heat_map_maximum(T_C, times, matrix)
            assert maximum[0] == pytest.approx(value, abs=2e-4)
            _, T_max, t_max = maximum
            assert red == (t_max, T_max)
            assert purple == (times[np.argmax(matrix[20])], 50.)
            assert [text.get_text() for text in panel.texts
                    if text.get_text().startswith('Max')] == [
                'Max %.3f\n%g \N{DEGREE SIGN}C, %.0f min' % maximum]
            # the red dashed line follows the row maxima
            (ridge,) = [line for line in panel.get_lines()
                        if len(line.get_xdata()) == 71]
            assert np.array_equal(ridge.get_xdata(),
                                  pkg.heat_map_ridge(times, matrix))
            # inset: a zoom of the same matrix on 47 to 57 C
            colorbar, inset = panel.child_axes
            assert np.array_equal(_mesh(inset), matrix)
            _assert_heat_map_lines(inset, T_C, times, matrix)
            assert inset.get_ylim() == (46.5, 57.5)
            # inset contour labels as published: every level on both sides
            # of the red dashed line, turned along the steep lines of the
            # zoom (within 20 degrees of upright before the red line)
            labelled = set()
            for text in inset.texts:
                t_label, T_label = text.get_position()
                early = t_label < np.interp(T_label, T_C, ridge.get_xdata())
                labelled.add((float(text.get_text()), early))
                if early:
                    assert abs(text.get_rotation() % 180. - 90.) < 20.
            assert labelled == {(level, early) for level in (0.1, 0.3, 0.5)
                                for early in (True, False)}
            assert colorbar.get_ylabel() == label
            assert colorbar.get_ylim() == (0., 1.)
    finally:
        plt.close(fig)


def test_fig11_verification_overlays_the_digitized_lines():
    fig = figures.plot_fig11_verification()
    try:
        tops = _heat_map_axes(fig)
        bottoms = [ax for ax in fig.axes
                   if ax.get_xlabel().startswith('Temperature of the sample')]
        assert len(tops) == len(bottoms) == 2
        maps = pkg.fig11_heat_maps()
        for top, bottom, read, matrix, rms, matched in zip(
                tops, bottoms, (pkg.read_fig11a_selectivity_digitized,
                                pkg.read_fig11b_yield_digitized), maps[2:],
                (1.05, 0.30), ('not matched', 'matched')):
            digitized = read()
            # both panels are the maps of the one default constant set
            assert np.array_equal(_mesh(top), matrix)
            T_C, times = maps[:2]
            _assert_heat_map_lines(top, T_C, times, matrix)
            # the digitized contour lines and red dashed line as black dots,
            # the digitized stars as open black stars
            for samples in (digitized['contours'], digitized['ridge']):
                (dots,) = [line for line in top.get_lines()
                           if len(line.get_xdata()) == len(samples)]
                assert np.array_equal(dots.get_xdata(), samples['time'])
                assert np.array_equal(dots.get_ydata(), samples['T_C'])
                assert to_hex(dots.get_color()) == '#000000'
            (stars,) = [line for line in top.get_lines()
                        if line.get_marker() == '*'
                        and line.get_markerfacecolor() == 'none']
            assert np.array_equal(stars.get_xdata(),
                                  digitized['stars']['time'])
            assert np.array_equal(stars.get_ydata(),
                                  digitized['stars']['T_C'])
            drawn = [line for line in top.get_lines()
                     if len(line.get_xdata()) == len(digitized['contours'])]
            # the digitized dots (a '.' marker is half its size wide) are
            # narrower than the lines of the map, which stay visible
            (contours,) = [c for c in top.collections
                           if isinstance(c, ContourSet)]
            assert drawn[0].get_marker() == '.'
            dot_width = (0.5*drawn[0].get_markersize()
                         + drawn[0].get_markeredgewidth())
            assert dot_width <= 0.5*min(contours.get_linewidths())
            # six groups of offsets (three levels, two sides). The misfit
            # of panel a (1.05 px rms, single samples up to 3.0 px) is
            # drawn and named: in the title, and in the offsets panel.
            groups = [line for line in bottom.get_lines()
                      if len(line.get_xdata()) > 2]
            assert len(groups) == 6
            # each group: the offsets of its own samples from the contour of
            # its level, in the colour of that level, filled on the late side
            # of the red dashed line and open on the early side
            for line, ((level, branch), samples) in zip(
                    groups, digitized['contours'].groupby(['level',
                                                           'branch'])):
                assert line.get_label() == '%g, %s' % (level, branch)
                inside = samples[samples['time'].between(times[0], times[-1])
                                 & samples['T_C'].between(T_C[0], T_C[-1])]
                assert np.array_equal(line.get_xdata(), inside['T_C'])
                assert np.array_equal(line.get_ydata(), pkg.contour_offsets(
                    T_C, times, matrix, level, inside['time'],
                    inside['T_C']))
                color = to_hex(figures._CONTOUR_COLORS[level])
                assert to_hex(line.get_markeredgecolor()) == color
                assert (line.get_marker(), line.get_markerfacecolor()) == (
                    ('o', 'none') if branch == 'early'
                    else ('s', line.get_markeredgecolor()))
            offsets = np.concatenate([line.get_ydata() for line in groups])
            assert np.sqrt(np.mean(offsets**2)) == pytest.approx(rms,
                                                                 abs=0.005)
            low, high = bottom.get_ylim()
            assert low < offsets.min() and offsets.max() < high
            assert np.abs(offsets).max() == pytest.approx(
                3.02 if matched == 'not matched' else 1.64, abs=0.01)
            assert top.get_title() == (
                'Published contours %s by the model: rms offset %.2f px'
                % (matched, rms))
            assert [text.get_text() for text in bottom.texts][0].startswith(
                'rms %.2f px' % rms)
            # key of the overlay: the model's lines and stars in their
            # colours, the published ones in black
            legend = top.get_legend()
            assert [text.get_text() for text in legend.get_texts()] == [
                'Model: contours 0.1, 0.3, 0.5', 'Model: row maxima',
                'Model: maximum, 50 \N{DEGREE SIGN}C maximum',
                'Published: lines (digitized)',
                'Published: stars (digitized)']
            assert [(to_hex(line.get_color()), line.get_marker())
                    for line in legend.findobj(Line2D)] == [
                (to_hex(figures._CONTOUR_COLORS[level]), 'None')
                for level in (0.1, 0.3, 0.5)] + [
                ('#e41a1c', 'None'), ('#f0525a', '*'), ('#7b6fc4', '*'),
                ('#000000', 'None'), ('#000000', '*')]
            assert bottom.get_legend() is not None
    finally:
        plt.close(fig)


def test_solver_emulation_figure_marks_the_emulated_steps_on_the_nodes():
    fig = figures.plot_solver_emulation()
    try:
        tops = [ax for ax in fig.axes if ax.get_title()]
        bottoms = [ax for ax in fig.axes if ax.get_yscale() == 'log']
        assert len(tops) == len(bottoms) == 3 and len(fig.axes) == 6
        assert [ax.get_title().split('\n')[1].split(';')[0]
                for ax in tops] == ['23 steps', '19 steps', '24 steps']
        colors = dict(figures._COLORS, eq_H2O2='#f28e8e')
        for top, bottom, nodes, km in zip(
                tops, bottoms,
                (pkg.read_fig9a_h2o2()['model'], pkg.read_fig9b_o2()['model'],
                 pkg.read_fig9a_h2o2()['h2o2_trace']),
                (pkg.load_h2o2_model('fig9a_nodes'),
                 pkg.load_o2_model('fig9b_nodes'),
                 pkg.load_h2o2_model('figS12a_h2o2_nodes'))):
            columns = list(nodes.columns[1:])
            circles = [line for line in top.get_lines()
                       if line.get_marker() == 'o']
            crosses = [line for line in top.get_lines()
                       if line.get_marker() == '+']
            assert len(circles) == len(crosses) == nodes.shape[1] - 1
            # per plotted quantity, in its colour: the deposited nodes and
            # the accurate solution of the model that carries the constants
            # inferred from them
            lines = [line for line in top.get_lines()
                     if line.get_marker() == 'None']
            assert len(lines) == len(columns)
            for circle, line, name in zip(circles, lines, columns):
                assert np.array_equal(circle.get_ydata(), nodes[name]), name
                assert to_hex(circle.get_markeredgecolor()) == colors[name]
                assert to_hex(line.get_color()) == colors[name]
                accurate = pkg.simulate_at_times(km, line.get_xdata(),
                                                 [name])
                assert np.allclose(line.get_ydata(), accurate[name],
                                   rtol=0., atol=1e-12), name
            # every deposited node carries a cross
            n = len(nodes)
            for circle, cross in zip(circles, crosses):
                assert np.array_equal(circle.get_xdata(), nodes['time'])
                assert np.abs(cross.get_xdata()[:n]
                              - nodes['time']).max() <= 2e-5
                assert np.abs(cross.get_ydata()[:n]
                              - circle.get_ydata()).max() <= 1e-5
            # below: the emulated steps lie under the line of the last
            # printed decimal, the accurate solution rises above it
            accurate, emulated = [line for line in bottom.get_lines()
                                  if line.get_marker() in ('o', '+')]
            assert (accurate.get_marker(), emulated.get_marker()) == ('o', '+')
            assert emulated.get_ydata().max() <= 1e-5
            assert accurate.get_ydata().max() > 1e-4
            low, high = bottom.get_ylim()
            assert low < emulated.get_ydata().min()
            assert accurate.get_ydata().max() < high
        (legend,) = fig.legends
        assert [(text.get_text(), to_hex(handle.get_color()))
                for text, handle in zip(legend.get_texts()[:6],
                                        legend.legend_handles)] == [
            (figures._LABELS[name], figures._COLORS[name])
            for name in pkg.NORMALIZED] + [('H$_2$O$_2$', '#f28e8e')]
    finally:
        plt.close(fig)


@pytest.mark.parametrize('font', ['DejaVu Sans', None])
def test_solver_emulation_titles_do_not_run_into_each_other(font,
                                                            monkeypatch):
    # With the last fallback font (DejaVu Sans, the widest) as well as with
    # the default font list, the two-line titles of neighbouring panels
    # stay apart and inside the figure (drawn inside the rc context of the
    # module, as save_figures draws them).
    if font is not None:
        monkeypatch.setitem(figures._RC_PARAMS, 'font.sans-serif', [font])
    with plt.rc_context(figures._RC_PARAMS):
        fig = figures.plot_solver_emulation()
        try:
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            boxes = [ax.title.get_window_extent(renderer)
                     for ax in fig.axes if ax.get_title()]
            assert len(boxes) == 3
            em = 11.5*fig.dpi/72.       # font size of the titles in pixels
            for left, right in zip(boxes, boxes[1:]):
                assert right.x0 - left.x1 > 2.*em
            assert boxes[0].x0 > 0. and boxes[-1].x1 < fig.bbox.width
        finally:
            plt.close(fig)


def test_save_figures_writes_the_figures_named(tmp_path):
    # the names the command line offers (--figures), each mapped to the
    # plot function of that name
    assert list(figures._FIGURES) == [
        'fig9', 'figS12', 'fig11', 'fig11_verification', 'solver_emulation']
    for name, plot in figures._FIGURES.items():
        assert plot is getattr(figures, 'plot_' + name), name
    paths = figures.save_figures(str(tmp_path), formats=('pdf',),
                                 figures=['fig9'])
    assert [os.path.basename(path) for path in paths] == [
        'pd_c_3hp_oxidation_fig9.pdf']
    assert os.path.getsize(paths[0]) > 0
    with pytest.raises(ValueError, match='unknown figures'):
        figures.save_figures(str(tmp_path), figures=['fig9c'])


def test_command_line_writes_the_figures_named(tmp_path, monkeypatch,
                                              capsys):
    # python -m nskinetics.models.pd_c_3hp_oxidation.figures OUTPUT_DIR ...
    monkeypatch.setattr('sys.argv', [
        'figures', str(tmp_path), '--figures', 'figS12', '--formats', 'png',
        '--dpi', '20'])
    figures.main()
    path = os.path.join(str(tmp_path), 'pd_c_3hp_oxidation_figS12.png')
    assert capsys.readouterr().out == f'Saved {path}\n'
    assert os.listdir(str(tmp_path)) == ['pd_c_3hp_oxidation_figS12.png']
    monkeypatch.setattr('sys.argv', ['figures', str(tmp_path), '--figures',
                                     'fig9c'])
    with pytest.raises(SystemExit):
        figures.main()


def test_save_figures_dpi_sets_the_size_of_the_png(tmp_path):
    # The images of the documentation page are written at dpi=100; the
    # default is 600.
    sizes = []
    for dpi in (40, 80):
        (path,) = figures.save_figures(
            os.path.join(str(tmp_path), str(dpi)), formats=('png',),
            figures=['figS12'], dpi=dpi)
        assert os.path.basename(path) == 'pd_c_3hp_oxidation_figS12.png'
        sizes.append(_png_size(path))
    assert sizes[1] == pytest.approx([2*size for size in sizes[0]], abs=2)
    assert figures.save_figures.__defaults__[-1] == 600


def test_figS12_has_outward_major_ticks_as_published():
    fig = figures.plot_figS12()
    try:
        for panel in fig.axes[:2]:
            assert panel.xaxis.get_major_ticks()[0].get_tickdir() == 'out'
            assert len(panel.xaxis.get_minorticklocs()) == 0
    finally:
        plt.close(fig)
