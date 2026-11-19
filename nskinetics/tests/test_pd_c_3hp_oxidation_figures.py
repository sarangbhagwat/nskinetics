# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Rendering tests for ``nskinetics.models.pd_c_3hp_oxidation.figures``
(Fig. 9 and SI Fig. S12 of Kim et al., Appl. Catal. B: Environ. Energy 386
(2026) 126403, regenerated from the NSKinetics models).

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
    assert snapshot() == before
