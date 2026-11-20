# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Rendering tests for ``nskinetics.models.pd_c_3hp_oxidation.figures``
(Fig. 9 and SI Fig. S12 of Kim et al., Appl. Catal. B: Environ. Energy 386
(2026) 126403, regenerated from the NSKinetics models, and the comparison
and verification figures drawn next to them).

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


def _png_size(path):
    with open(path, 'rb') as f:
        header = f.read(24)
    return struct.unpack('>II', header[16:24])


def test_fig9_draws_model_experiment_and_deposited_curve():
    fig = figures.plot_fig9()
    try:
        panel_a, reactions_a, panel_b, reactions_b = fig.axes
        for panel, read in ((panel_a, pkg.read_fig9a_h2o2),
                            (panel_b, pkg.read_fig9b_o2)):
            # 5 model lines + 5 experimental series + 5 deposited polylines
            assert len(_lines_with_data(panel)) == 15
            assert panel.get_xlim() == (-30., 1440.)
            assert panel.get_xlabel() == 'Time (min)'
            # markers are drawn whole on the frame, so no sample taken
            # after the end of the time axis may be drawn (Fig. 9b has
            # samples at 1470 and 1680 min)
            experiment = read()['experiment']
            drawn = experiment[experiment['time'] <= figures._T_END]
            markers = _marker_only_lines(panel)
            assert len(markers) == 5
            for line in markers:
                assert np.array_equal(line.get_xdata(), drawn['time'])
            deposited = read()['model']
            drawn = [line for line in panel.get_lines()
                     if len(line.get_xdata()) == len(deposited)]
            assert len(drawn) == 5
            assert np.array_equal(drawn[0].get_ydata(), deposited['X_HP3'])
            # inward ticks with minor ticks half-way, as published
            assert list(panel.get_xticks()) == [0, 500, 1000]
            assert 250. in panel.xaxis.get_minorticklocs()
            assert panel.xaxis.get_tick_params()['direction'] == 'in'
            assert panel.yaxis.get_tick_params(which='minor')[
                'direction'] == 'in'
        # five pieces of text per reaction: "Rxn n" (navy, as published),
        # reactants, arrow, k over the arrow (red), products
        for reactions, n in ((reactions_a, 7), (reactions_b, 6)):
            texts = reactions.texts
            assert len(texts) == 5*n
            assert [text.get_text() for text in texts[::5]] == [
                'Rxn %d' % i for i in range(1, n + 1)]
            assert [text.get_text() for text in texts[3::5]] == [
                'k$_%d$' % i for i in range(1, n + 1)]
            assert {to_hex(text.get_color()) for text in texts[::5]} == {
                '#3653a0'}
            assert {to_hex(text.get_color()) for text in texts[3::5]} == {
                '#ee2c32'}
        assert reactions_a.texts[1].get_text() == '3-HP + H$_2$O$_2$'
        assert reactions_b.texts[4].get_text() == 'OPA + H$_2$O'
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
        # (a) species, published curve, deposited trace, predicted H2O2
        h2o2 = inset_a.get_lines()[-1]
        assert h2o2.get_ydata()[0] == pytest.approx(15.)
        # (b) published polyline and predicted OPA, on the 1e-5 axis
        published, opa = inset_b.get_lines()
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
    assert snapshot() == before


# --- Table 3 O2 model against the experiments of Fig. 10a-c; Fig. 10d ---------

def _two_point_lines(ax):
    return [line for line in ax.get_lines() if len(line.get_xdata()) == 2]


def _slope_and_value_at_100C(line):
    """Slope in K and value at 1/373.15 K of a line drawn against 1000/T."""
    x, y = line.get_xdata()/1000., line.get_ydata()
    slope = (y[1] - y[0])/(x[1] - x[0])
    return slope, y[0] + slope*(1./373.15 - x[0])


def test_fig10_table3_o2_draws_the_comparison(tmp_path):
    fig = figures.plot_fig10_table3_o2()
    try:
        panels = fig.axes
        assert len(panels) == 4
        # pooled R2 of the Table 3 model
        for panel, temperature_C, r2 in zip(
                panels, pkg.FIG10_TEMPERATURES_C,
                ('0.579', '0.644', '0.746')):
            # 5 Table 3 lines + 5 experimental series
            lines = _lines_with_data(panel)
            assert len(lines) == 10
            # the samples after the end of the time axis (1740 and 1890 min
            # at 25 C, 1560 min at 100 C) are not drawn
            experiment = pkg.read_fig10_o2(temperature_C)['experiment']
            drawn = experiment[experiment['time'] <= figures._T_END]
            assert len(experiment) - len(drawn) == {25: 2, 75: 0,
                                                    100: 1}[temperature_C]
            markers = _marker_only_lines(panel)
            assert len(markers) == 5
            for line in markers:
                assert np.array_equal(line.get_xdata(), drawn['time'])
            times = lines[0].get_xdata()[[0, 480, 2880]]
            expected = pkg.simulate_at_times(
                pkg.load_o2_model(T=temperature_C + 273.15), times,
                ['X_HP3'])
            assert np.allclose(lines[0].get_ydata()[[0, 480, 2880]],
                               expected['X_HP3'], rtol=0., atol=1e-7)
            assert [text.get_text() for text in panel.texts][:2] == [
                '%d \N{DEGREE SIGN}C' % temperature_C, 'R$^2$ = ' + r2]
        # the keys of panels a to c are above the panels
        assert len(fig.legends) == 2
        # panel d: the three series of four points, as deposited
        arrhenius = panels[3]
        assert arrhenius.get_xlim() == (2.6, 3.45)
        assert arrhenius.get_ylim() == (-13., -3.)
        assert arrhenius.get_xlabel() == '1000/T (K$^{-1}$)'
        points = pkg.fig10d_points()
        series = [line for line in arrhenius.get_lines()
                  if len(line.get_xdata()) == 4]
        assert len(series) == 3
        for line, name in zip(series, ('k_1', 'k_2', 'k_5')):
            assert np.allclose(line.get_xdata(), 1000.*points['inv_T'])
            assert np.array_equal(line.get_ydata(), points['ln_' + name])
        # in the colors of the published series: those of MA, AAc and Others
        assert [to_hex(line.get_markeredgecolor()) for line in series] == [
            color for _, color in figures._FIG10D_SERIES.values()] == [
            '#a66bd6', '#1f6fdb', '#6f9bc4']
        # per constant, the regression line (from 100 C to the end of the
        # published line) and the line of the activation energy quoted in
        # the text through the 50 C point (from 100 C to 25 C)
        lines = _two_point_lines(arrhenius)
        assert len(lines) == 6
        for (regression, text), name, Ea, (slope, at_100C, T_end) in zip(
                zip(lines[::2], lines[1::2]), ('k_1', 'k_2', 'k_5'),
                (59781.2, 91594.6, 59549.8),
                ((-7208.69, -4.61102, 25.), (-8293.97, -6.32951, 25.13),
                 (-6265.49, -4.17048, 23.50))):
            assert regression.get_linestyle() == '--'
            assert _slope_and_value_at_100C(regression) == pytest.approx(
                (slope, at_100C), rel=1e-4, abs=2e-5)
            assert list(1000./regression.get_xdata() - 273.15) == (
                pytest.approx([100., T_end], abs=1e-9))
            assert list(1000./text.get_xdata() - 273.15) == (
                pytest.approx([100., 25.], abs=1e-9))
            slope, at_100C = _slope_and_value_at_100C(text)
            assert -8.314*slope == pytest.approx(Ea, rel=1e-9)
            assert at_100C + slope*(1./323.15 - 1./373.15) == pytest.approx(
                points['ln_' + name][1], abs=1e-9)
        assert [text.get_text()
                for text in arrhenius.get_legend().get_texts()] == [
            '3-HP $\\rightarrow$ MA', 'MA $\\rightarrow$ AAc',
            'MA $\\rightarrow$ Others', 'Regression line (published)',
            'E$_a$ of the text, from 50 \N{DEGREE SIGN}C']
        path = os.path.join(str(tmp_path), 'fig10_table3_o2.png')
        fig.savefig(path, dpi=50)
        assert os.path.getsize(path) > 0
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
        for top, bottom, nodes in zip(
                tops, bottoms,
                (pkg.read_fig9a_h2o2()['model'], pkg.read_fig9b_o2()['model'],
                 pkg.read_fig9a_h2o2()['h2o2_trace'])):
            circles = [line for line in top.get_lines()
                       if line.get_marker() == 'o']
            crosses = [line for line in top.get_lines()
                       if line.get_marker() == '+']
            assert len(circles) == len(crosses) == nodes.shape[1] - 1
            # every deposited node carries a cross
            n = len(nodes)
            for circle, cross in zip(circles, crosses):
                assert np.array_equal(circle.get_xdata(), nodes['time'])
                assert np.abs(cross.get_xdata()[:n]
                              - nodes['time']).max() <= 0.01
                assert np.abs(cross.get_ydata()[:n]
                              - circle.get_ydata()).max() <= 1e-5
            # below: the emulated steps lie under the line of the last
            # printed decimal, the accurate solution rises above it
            accurate, emulated = [line for line in bottom.get_lines()
                                  if line.get_marker() in ('o', '+')]
            assert emulated.get_ydata().max() <= 1e-5
            assert accurate.get_ydata().max() > 1e-4
            low, high = bottom.get_ylim()
            assert low < emulated.get_ydata().min()
            assert accurate.get_ydata().max() < high
        assert len(fig.legends) == 1
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
    # the names the command line offers (--figures)
    assert list(figures._FIGURES) == ['fig9', 'fig10_table3_o2', 'figS12',
                                      'solver_emulation']
    paths = figures.save_figures(str(tmp_path), formats=('pdf',),
                                 figures=['fig9'])
    assert [os.path.basename(path) for path in paths] == [
        'pd_c_3hp_oxidation_fig9.pdf']
    assert os.path.getsize(paths[0]) > 0
    with pytest.raises(ValueError, match='unknown figures'):
        figures.save_figures(str(tmp_path), figures=['fig9c'])


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
