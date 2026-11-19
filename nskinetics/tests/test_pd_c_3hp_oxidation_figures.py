# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Rendering tests for ``nskinetics.models.pd_c_3hp_oxidation.figures``
(Fig. 9, SI Fig. S12 and Fig. 11 of Kim et al., Appl. Catal. B: Environ.
Energy 386 (2026) 126403, regenerated from the NSKinetics models).

The figures are built in memory and their artists inspected, so a panel that
silently drew nothing would fail; one small PNG is written to ``tmp_path``.
The models are small, so the module stays in the default, non-slow run.

``matplotlib.use('Agg')`` is called at import time; this module is therefore
deliberately not registered in ``nskinetics/tests/__init__.py``.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.contour import ContourSet

import nskinetics.models.pd_c_3hp_oxidation as pkg
from nskinetics.models.pd_c_3hp_oxidation import figures


def _lines_with_data(ax):
    return [line for line in ax.get_lines() if len(line.get_xdata()) > 1]


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
            deposited = read()['model']
            drawn = [line for line in panel.get_lines()
                     if len(line.get_xdata()) == len(deposited)]
            assert len(drawn) == 5
            assert np.array_equal(drawn[0].get_ydata(), deposited['X_HP3'])
        # one line of text per reaction (label + equation)
        assert len(reactions_a.texts) == 2*7
        assert len(reactions_b.texts) == 2*6
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


def test_fig11_draws_both_heat_maps_with_their_own_constant_sets():
    fig = figures.plot_fig11()
    try:
        panel_a, panel_b = _heat_map_axes(fig)
        T_C, times, selectivity, MA_yield = pkg.fig11_heat_maps()
        for panel, matrix, label, value in (
                (panel_a, selectivity, 'MA Selectivity', 0.5906),
                (panel_b, MA_yield, 'MA Yield', 0.5483)):
            assert np.array_equal(_mesh(panel), matrix)
            assert panel.get_xlim() == (7.2, 1682.8)
            assert panel.get_ylim() == (29.5, 100.5)
            # red star on the largest cell, purple star on the 50 C row maximum
            red, purple = _star_positions(panel)
            maximum = pkg.heat_map_maximum(T_C, times, matrix)
            assert maximum[0] == pytest.approx(value, abs=0.003)
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
        for top, bottom, read in zip(
                tops, bottoms, (pkg.read_fig11a_selectivity_digitized,
                                pkg.read_fig11b_yield_digitized)):
            digitized = read()
            drawn = [line for line in top.get_lines()
                     if len(line.get_xdata()) == len(digitized['contours'])]
            assert len(drawn) == 1
            assert np.array_equal(drawn[0].get_ydata(),
                                  digitized['contours']['T_C'])
            # the digitized dots (a '.' marker is half its size wide) are
            # narrower than the lines of the map, which stay visible
            (contours,) = [c for c in top.collections
                           if isinstance(c, ContourSet)]
            assert drawn[0].get_marker() == '.'
            dot_width = (0.5*drawn[0].get_markersize()
                         + drawn[0].get_markeredgewidth())
            assert dot_width <= 0.5*min(contours.get_linewidths())
            # six groups of offsets (three levels, two sides), all small
            groups = [line for line in bottom.get_lines()
                      if len(line.get_xdata()) > 2]
            assert len(groups) == 6
            offsets = np.concatenate([line.get_ydata() for line in groups])
            assert np.sqrt(np.mean(offsets**2)) <= 0.45
            assert np.abs(offsets).max() < 1.5      # all inside the axes
            assert top.get_legend() is not None
            assert bottom.get_legend() is not None
    finally:
        plt.close(fig)


def test_fig11_table3_o2_draws_the_table3_maps(tmp_path):
    fig = figures.plot_fig11_table3_o2()
    try:
        panel_a, panel_b = _heat_map_axes(fig)
        T_C, times, selectivity, MA_yield = pkg.table3_o2_heat_maps()
        for panel, matrix in ((panel_a, selectivity), (panel_b, MA_yield)):
            assert np.array_equal(_mesh(panel), matrix)
            assert panel.get_legend() is not None
            _, T_max, t_max = pkg.heat_map_maximum(T_C, times, matrix)
            assert _star_positions(panel)[0] == (t_max, T_max)
        path = os.path.join(str(tmp_path), 'fig11_table3_o2.png')
        fig.savefig(path, dpi=50)
        assert os.path.getsize(path) > 0
    finally:
        plt.close(fig)
