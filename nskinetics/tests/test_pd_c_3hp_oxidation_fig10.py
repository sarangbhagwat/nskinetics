# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests of the Fig. 10 data (Kim et al., Appl. Catal. B: Environ. Energy
386 (2026) 126403) in ``nskinetics.models.pd_c_3hp_oxidation``: the tables
deposited for the O2 experiments at 25, 75 and 100 C (Fig. 10a-c) and for
the Arrhenius plot (Fig. 10d) with their readers, and the numbers
``NOTES.md`` reports for the Table 3 O2 model against those experiments.

The straight lines of the Arrhenius plot (Fig. 10d) are recomputed from the
deposited points and compared with the deposited lines and with the
activation energies and the 50 C constants of the O2 model.

The models are small, so this module is not marked slow. It imports pytest
at top level and is therefore not registered in
``nskinetics/tests/__init__.py``.
"""
import hashlib
import itertools
import os

import numpy as np
import pandas as pd
import pytest

import nskinetics.models.pd_c_3hp_oxidation as pkg
from nskinetics.models.pd_c_3hp_oxidation import (
    NORMALIZED, load_o2_model, pooled_r2, rms_difference, simulate_at_times)

_NORMALIZED = list(NORMALIZED)
TEMPERATURES_C = pkg.FIG10_TEMPERATURES_C


# --- 1. reference data of Fig. 10a-c ---------------------------------------------

def test_deposited_fig10_files_are_verbatim_copies():
    # SHA-256 of the files in the Dryad deposit (doi:10.5061/dryad.jsxksn0q4)
    for filename, sha256 in (
            ('11-1_25C_Kinetic.csv',
             'f6037edf98df13acbbc74455cb9009bdf376cb08e5203ca7b7a4738c4f8e2892'),
            ('11-2_75C_Kinetic.csv',
             'c3fa90a544c958e71ea3f006aa53f2bfbf0367524fbff0939cc9105479557c70'),
            ('11-3_100C_Kinetic.csv',
             '78b790344240f007237d8cb8f89ebd15c081d9018147647783f3d64bf09b3a76')):
        with open(os.path.join(pkg.reference_data_dir, filename), 'rb') as f:
            assert hashlib.sha256(f.read()).hexdigest() == sha256, filename


def test_fig10_reader_returns_experiment_and_model_blocks():
    assert TEMPERATURES_C == (25, 75, 100)
    rows = {25: (14, 19), 75: (12, 20), 100: (11, 25)}
    last_sample = {25: 1890., 75: 1440., 100: 1560.}
    for temperature_C in TEMPERATURES_C:
        tables = pkg.read_fig10_o2(temperature_C)
        assert list(tables) == ['experiment', 'model']
        experiment, model = tables['experiment'], tables['model']
        assert (len(experiment), len(model)) == rows[temperature_C]
        for table in (experiment, model):
            assert list(table.columns) == ['time'] + _NORMALIZED
            assert list(table.iloc[0]) == [0., 1., 0., 0., 0., 0.]
            assert (np.diff(table['time']) > 0.).all()
            assert np.allclose(table[_NORMALIZED].sum(axis=1), 1., atol=5e-5)
        # sampling times are whole minutes; the curve ends 10 % later
        assert (experiment['time'] % 10 == 0).all()
        assert experiment['time'].iloc[-1] == last_sample[temperature_C]
        assert model['time'].iloc[-1] == pytest.approx(
            1.1*last_sample[temperature_C])
    assert pkg.read_fig10_o2(25)['model']['time'].iloc[1] == 0.38
    assert pkg.read_fig10_o2(100)['experiment']['X_MA'].max() == 0.56794


def test_fig10_reader_rejects_other_temperatures():
    # the argument is a temperature in degrees Celsius; 348.15 is 75 C in K.
    # An argument that cannot be looked up, such as a list, is refused with
    # the same error.
    for temperature_C in (50, 20, '25', '50C', '75 C', 348.15, None, True,
                          [25], np.array([75])):
        with pytest.raises(ValueError, match='no Fig. 10 data') as error:
            pkg.read_fig10_o2(temperature_C)
        assert 'read_fig9b_o2' in str(error.value)


def test_experiments_of_fig10_as_the_text_and_the_notes_quote_them():
    # Largest measured MA fraction and its time. The text gives maximum MA
    # yields of 52.6 % (75 C) and 56.8 % (100 C) and 2.3 % after 31.5 h at
    # 25 C: the 75 C value is not in the deposited data.
    for temperature_C, value, time in ((25, 0.02335, 1890.),
                                       (75, 0.43581, 180.),
                                       (100, 0.56794, 40.)):
        experiment = pkg.read_fig10_o2(temperature_C)['experiment']
        largest = experiment['X_MA'].idxmax()
        assert (experiment['X_MA'][largest], experiment['time'][largest]) == (
            value, time)
    # The published panels end at 1440 min; these samples are later.
    late = {T: [time for time in pkg.read_fig10_o2(T)['experiment']['time']
                if time > 1440.] for T in TEMPERATURES_C}
    assert late == {25: [1740., 1890.], 75: [], 100: [1560.]}
    # Measured X_HP3 at 60, 90 and 120 min: it falls to 0.5 between 60 and
    # 90 min at 75 C as at 50 C (Fig. 9b).
    for experiment, left in (
            (pkg.read_fig10_o2(75)['experiment'], [0.533, 0.493, 0.283]),
            (pkg.read_fig9b_o2()['experiment'], [0.522, 0.443, 0.255])):
        rows = experiment[experiment['time'].isin([60., 90., 120.])]
        assert list(rows['X_HP3']) == pytest.approx(left, abs=5e-4)


# --- 2. Arrhenius plot (Fig. 10d) ------------------------------------------------

# The plotted rate constants, taken to be in 1/min (NOTES.md) and named as
# read_fig10d_arrhenius names them (3-HP -> MA, MA -> AAc, MA -> D): the
# shortest decimals whose logarithms are the deposited points.
FIG10D_K = {'k_1': (7.75e-5, 2.14e-3, 2.43685e-3, 1.00763e-2),
            'k_2': (6.64e-4, 5.2e-5, 4.44096e-4, 1.59645e-3),
            'k_5': (1.676e-3, 1.0625e-3, 5.4797e-3, 1.41038e-2)}
# Slope (K), intercept (NOTES.md gives it to four decimals), value at
# 1/373.15 K (the first row of the deposited line) and activation energy
# (kJ/mol) of the three lines.
FIG10D_LINES = {'k_1': (-7208.69, 14.707468, -4.61102, 59.93),
                'k_2': (-8293.97, 15.897396, -6.32951, 68.96),
                'k_5': (-6265.49, 12.620339, -4.17048, 52.09)}
# The recomputed line minus the first row of the deposited line: the two
# agree to 7e-6 at 100 C.
FIG10D_DIFFERENCE_AT_100C = {'k_1': 7.2e-6, 'k_2': -2.2e-6, 'k_5': 2.7e-6}
# The steps of the O2 model that correspond to the three plotted constants.
O2_MODEL_STEPS = {'k_1': ('1', '2'), 'k_2': ('3',), 'k_5': ('6',)}


def test_arrhenius_table_reader():
    tables = pkg.read_fig10d_arrhenius()
    assert list(tables) == ['points', 'line_k_1', 'line_k_2', 'line_k_5']
    points = tables['points']
    assert list(points.columns) == ['inv_T', 'ln_k_1', 'ln_k_2', 'ln_k_5']
    assert list(points['inv_T']) == [0.00335, 0.00309, 0.00287, 0.00268]
    assert list(points.iloc[1]) == [0.00309, -6.14695, -9.86427, -6.84713]
    assert [len(tables[name]) for name in ('line_k_1', 'line_k_2',
                                           'line_k_5')] == [10, 18, 20]
    assert list(tables['line_k_5'].iloc[-1]) == [0.00337, -8.50046]
    # The file gives a unit for the abscissa only; that of the constants is
    # inferred (NOTES.md).
    with open(os.path.join(pkg.reference_data_dir, '11-4_Arrhenius.csv'),
              encoding='utf-8-sig') as f:
        labels, units = (f.readline().strip().split(',') for _ in range(2))
    assert labels[:4] == ['1/T', 'ln k1', 'ln k2', 'ln k3']
    assert set(units) == {'/K', 'exp', 'fitting', ''}


def test_fig10d_points_are_the_deposited_ones():
    points = pkg.fig10d_points()
    assert list(points['T_C']) == [25., 50., 75., 100.]
    assert np.allclose(points['inv_T'], 1./(points['T_C'] + 273.15),
                       rtol=1e-14)
    deposited = pkg.read_fig10d_arrhenius()['points']
    # the deposited abscissa is this one, printed to five decimals
    assert list(np.round(points['inv_T'], 5)) == list(deposited['inv_T'])
    for name, k in FIG10D_K.items():
        assert list(points['ln_' + name]) == list(deposited['ln_' + name])
        # the largest difference is 4.9e-6 (k_5)
        assert np.abs(np.log(k) - points['ln_' + name]).max() <= 5e-6


def test_fig10d_lines_are_least_squares_through_three_points():
    assert pkg.FIG10D_REGRESSION_TEMPERATURES_C == {
        'k_1': (25, 75, 100), 'k_2': (50, 75, 100), 'k_5': (50, 75, 100)}
    lines = pkg.fig10d_regression_lines()
    deposited = pkg.read_fig10d_arrhenius()
    for name, (slope, intercept, at_100C, Ea) in FIG10D_LINES.items():
        line = lines[name]
        assert line['slope'] == pytest.approx(slope, rel=1e-4)
        assert line['intercept'] == pytest.approx(intercept, abs=2e-6)
        assert line['Ea'] == -8.314*line['slope']
        assert line['Ea']/1000. == pytest.approx(Ea, abs=0.005)
        assert line['temperatures_C'] == (
            pkg.FIG10D_REGRESSION_TEMPERATURES_C[name])
        drawn = deposited['line_' + name]['ln_' + name]
        assert drawn.iloc[0] == at_100C
        assert (line['slope']/373.15 + line['intercept'] - at_100C
                == pytest.approx(FIG10D_DIFFERENCE_AT_100C[name], abs=2e-7))
        # The published line ends where the deposited line ends: at the
        # temperature at which the recomputed line takes the ln k of the
        # last row, to two decimals.
        T_end = pkg.FIG10D_LINE_ENDS_C[name][1]
        assert (line['slope']/(T_end + 273.15) + line['intercept']
                == pytest.approx(drawn.iloc[-1], abs=2e-4))
        assert (line['slope']/(drawn.iloc[-1] - line['intercept']) - 273.15
                == pytest.approx(T_end, abs=0.005))
    # The 1/T printed in the last rows reads 25.36, 25.36 and 23.59 C.
    assert [deposited['line_' + name]['inv_T'].iloc[-1]
            for name in FIG10D_LINES] == [0.00335, 0.00335, 0.00337]
    assert (1./0.00335 - 273.15, 1./0.00337 - 273.15) == pytest.approx(
        (25.3575, 23.5859), abs=1e-3)
    # at 25 C the line of k_1 is 4.5e-7 below the last row of the deposited
    # one
    k_1 = lines['k_1']
    assert k_1['slope']/298.15 + k_1['intercept'] + 9.47060 == pytest.approx(
        -4.5e-7, abs=5e-8)
    # At the abscissae as printed in the file the recomputed lines differ
    # from the deposited ones by up to 0.04 (0.035, 0.041 and 0.042).
    at_printed = [np.abs(
        lines[name]['slope']*deposited['line_' + name]['inv_T']
        + lines[name]['intercept']
        - deposited['line_' + name]['ln_' + name]).max()
        for name in FIG10D_LINES]
    assert at_printed == pytest.approx([0.0345, 0.0406, 0.0419], abs=5e-4)
    # The file prints 1/T to five decimals. A regression of each deposited
    # line on the abscissa as printed therefore gives other slopes:
    printed = [-8.314e-3*np.polyfit(deposited['line_' + name]['inv_T'],
                                    deposited['line_' + name]['ln_' + name],
                                    1)[0] for name in FIG10D_LINES]
    assert printed == pytest.approx([59.97, 69.32, 52.38], abs=0.005)


def test_each_fig10d_line_leaves_out_one_point():
    points, lines = pkg.fig10d_points(), pkg.fig10d_regression_lines()
    # ln k of the point minus the line, at 25, 50, 75 and 100 C
    for name, T_out, offset in (('k_1', 50., 1.45), ('k_2', 25., 4.60),
                                ('k_5', 25., 2.00)):
        residual = points['ln_' + name] - (
            lines[name]['slope']*points['inv_T'] + lines[name]['intercept'])
        left_out = points['T_C'] == T_out
        assert residual[left_out].iloc[0] == pytest.approx(offset, abs=0.01)
        assert np.abs(residual[~left_out]).max() < 0.21
        if name == 'k_1':
            # the 50 C point of 3-HP -> MA: a factor of 4.3 above the line
            assert np.exp(residual[left_out].iloc[0]) == pytest.approx(
                4.28, abs=0.005)
    # No other subset of two, three or four points gives the first row of
    # a deposited line (its value at 100 C); the nearest is 0.013 off.
    deposited = pkg.read_fig10d_arrhenius()
    inv_T = points['inv_T'].to_numpy()
    nearest_other = np.inf
    for name, temperatures_C in pkg.FIG10D_REGRESSION_TEMPERATURES_C.items():
        ln_k = points['ln_' + name].to_numpy()
        at_100C = deposited['line_' + name]['ln_' + name].iloc[0]
        for size in (2, 3, 4):
            for rows in itertools.combinations(range(4), size):
                slope, intercept = np.polyfit(inv_T[list(rows)],
                                              ln_k[list(rows)], 1)
                fitted = tuple(int(points['T_C'][i]) for i in rows)
                difference = abs(slope/373.15 + intercept - at_100C)
                assert (difference < 0.01) == (fitted == temperatures_C), (
                    name, fitted)
                if fitted != temperatures_C:
                    nearest_other = min(nearest_other, difference)
    assert nearest_other == pytest.approx(0.01345, abs=1e-5)


def test_activation_energies_of_the_text_are_not_those_of_the_lines():
    assert pkg.TEXT_ACTIVATION_ENERGIES == {
        'k_1': 59781.2, 'k_2': 91594.6, 'k_5': 59549.8}
    points = pkg.fig10d_points()
    regression, text = pkg.fig10d_regression_lines(), pkg.fig10d_text_lines()
    for name, Ea in pkg.TEXT_ACTIVATION_ENERGIES.items():
        line = text[name]
        assert line['Ea'] == pytest.approx(Ea, rel=1e-12)
        # through the 50 C point
        assert line['slope']/323.15 + line['intercept'] == pytest.approx(
            points['ln_' + name][1], abs=1e-12)
        assert line['temperatures_C'] == (50,)
    difference = [(regression[name]['Ea'] - text[name]['Ea'])/1000.
                  for name in ('k_1', 'k_2', 'k_5')]
    assert difference == pytest.approx([0.15, -22.64, -7.46], abs=0.01)
    # ln k of the other three points minus the text line through 50 C
    for name, offsets in (('k_1', (-1.45, -1.47, -1.43)),
                          ('k_2', (5.41, -0.30, -1.14)),
                          ('k_5', (2.31, 0.05, -0.38))):
        residual = points['ln_' + name] - (
            text[name]['slope']*points['inv_T'] + text[name]['intercept'])
        assert list(residual[[0, 2, 3]]) == pytest.approx(offsets, abs=0.01)


def test_no_subset_of_the_points_gives_the_activation_energies_of_the_text():
    points = pkg.fig10d_points()
    nearest = {}
    for name, Ea in pkg.TEXT_ACTIVATION_ENERGIES.items():
        subsets = {}
        for size in (2, 3, 4):
            for rows in itertools.combinations(range(4), size):
                slope = np.polyfit(points['inv_T'].to_numpy()[list(rows)],
                                   points['ln_' + name].to_numpy()[list(rows)],
                                   1)[0]
                subsets[tuple(int(points['T_C'][i]) for i in rows)] = (
                    -8.314*slope/1000.)
        best = min(subsets, key=lambda rows: abs(subsets[rows] - Ea/1000.))
        nearest[name] = (best, round(subsets[best], 1))
        if name == 'k_5':
            # close to the heat-map script's 51788.73 J/mol for 3-HP -> D
            assert subsets[50, 100] == pytest.approx(51.85, abs=0.01)
    assert nearest == {'k_1': ((25, 75, 100), 59.9), 'k_2': ((50, 75), 80.2),
                       'k_5': ((50, 75), 61.4)}


def test_text_lines_carry_the_arrhenius_factors_of_the_o2_model():
    # The O2 model has the activation energies of the text for R1 and R2
    # (3-HP -> MA by way of OPA), R3 (MA -> AAc) and R6 (MA -> D), and its
    # 3-HP -> D step the value the 50 and 100 C points of MA -> D come
    # close to (previous test).
    km = load_o2_model()
    for name, steps in O2_MODEL_STEPS.items():
        for step in steps:
            assert km.get_value('Ea_' + step) == (
                pkg.TEXT_ACTIVATION_ENERGIES[name])
    assert km.get_value('Ea_5') == 51788.73
    # Each line of fig10d_text_lines is therefore the plotted 50 C point
    # multiplied by the Arrhenius factor k(T)/k_ref of the model.
    points, text = pkg.fig10d_points(), pkg.fig10d_text_lines()
    for temperature_C in TEMPERATURES_C:
        km = load_o2_model(T=temperature_C + 273.15)
        for name, steps in O2_MODEL_STEPS.items():
            rise = (text[name]['slope']/(temperature_C + 273.15)
                    + text[name]['intercept'] - points['ln_' + name][1])
            for step in steps:
                factor = (km.get_value('k_' + step)
                          / km.get_value('k_%s_ref' % step))
                assert rise == pytest.approx(np.log(factor), abs=1e-12)


def test_fig10d_50C_points_against_the_table_3_constants():
    # First-order constants of the Table 3 O2 model at 50 C in 1/min, for
    # the 0.2 g batch: 3-HP -> OPA (R2, OPA -> MA, is 19.5 times faster),
    # MA -> AAc and MA -> D.
    km = load_o2_model()
    m_cat = km.get_value('m_cat')
    scale = np.array([m_cat*km.get_value('x_O2_ref'), m_cat, m_cat])
    k = np.array([km.get_value('k_' + step) for step in '136'])
    assert list(k) == [7.14, 3.16e-4, 5.52e-3]
    model = k*scale
    assert model == pytest.approx([7.261e-3, 6.32e-5, 1.104e-3], rel=1e-3)
    assert km.get_value('k_2')/km.get_value('k_1') == pytest.approx(19.5,
                                                                    abs=0.05)
    # 95 % confidence intervals of Table 3, 7.14 (+- 0.4), 3.16 (+- 1.4)e-4
    # and 5.52 (+- 0.3)e-3, read as half-widths and scaled like the
    # constants
    half_width = np.array([0.4, 1.4e-4, 0.3e-3])
    low, high = (k - half_width)*scale, (k + half_width)*scale
    assert low == pytest.approx([6.855e-3, 3.52e-5, 1.044e-3], rel=1e-3)
    assert high == pytest.approx([7.668e-3, 9.12e-5, 1.164e-3], rel=1e-3)
    plotted = np.array([FIG10D_K[name][1] for name in O2_MODEL_STEPS])
    assert list(plotted) == [2.14e-3, 5.2e-5, 1.0625e-3]
    # Taken as they stand, the plotted constants of MA -> AAc and MA -> D
    # lie inside the intervals and that of 3-HP -> MA does not: the model
    # has 3.39, 1.22 and 1.04 times the plotted values, 1.22, 0.20 and 0.04
    # in ln k.
    assert plotted/model == pytest.approx([0.295, 0.823, 0.962], abs=0.0005)
    assert model/plotted == pytest.approx([3.39, 1.22, 1.04], abs=0.005)
    assert np.log(model/plotted) == pytest.approx([1.22, 0.20, 0.04],
                                                  abs=0.005)
    assert list((low < plotted) & (plotted < high)) == [False, True, True]
    assert np.log(low/model) == pytest.approx([-0.06, -0.59, -0.06],
                                              abs=0.005)
    assert np.log(high/model) == pytest.approx([0.05, 0.37, 0.05],
                                               abs=0.005)
    # In the deposited heat-map script the 3-HP -> MA constant multiplies
    # an O2 activity P/f0 = 3.5/1.01325 and the other two constants the MA
    # mole fraction alone. On that basis the plotted value of 3-HP -> MA
    # corresponds to 7.39e-3 1/min, 1.02 times the model's, and all three
    # lie inside the intervals.
    with_activity = plotted*np.array([3.5/1.01325, 1., 1.])
    assert with_activity[0] == pytest.approx(7.39e-3, rel=1e-3)
    assert with_activity/model == pytest.approx([1.018, 0.823, 0.962],
                                                abs=0.0005)
    assert np.all((low < with_activity) & (with_activity < high))
    # read as percentages, the +- values would contain none of them
    percent = np.array([0.4, 1.4, 0.3])/100.
    assert not np.any(np.abs(plotted/model - 1.) < percent)
    assert not np.any(np.abs(with_activity/model - 1.) < percent)
    # The unit of the plotted constants is inferred (NOTES.md): read as 1/s
    # or 1/h instead of 1/min, none would lie inside its interval.
    for factor in (60., 1./60.):
        for values in (plotted, with_activity):
            assert not np.any((low < factor*values) & (factor*values < high))


# --- 3. Table 3 O2 model against the experiments at 25, 75 and 100 C (NOTES.md) --

def test_rms_difference_and_pooled_r2_by_hand():
    # Two rows and five species, written out: residuals 0.1 (twice) and
    # -0.2 (once), all others zero.
    observed = pd.DataFrame([[1., 0., 0., 0., 0.], [0.5, 0.3, 0.1, 0., 0.1]],
                            columns=_NORMALIZED)
    modeled = pd.DataFrame([[0.9, 0.1, 0., 0., 0.], [0.5, 0.3, 0.1, 0., -0.1]],
                           columns=_NORMALIZED)
    ss_res = 0.01 + 0.01 + 0.04
    assert rms_difference(observed, modeled) == pytest.approx(
        np.sqrt(ss_res/10.), rel=1e-12)
    # sum of squares about the mean of the ten observed values, 0.2
    ss_tot = 0.64 + 0.09 + 0.01 + 0.01 + 0.01 + 5*0.04
    assert pooled_r2(observed, modeled) == pytest.approx(1. - ss_res/ss_tot,
                                                         rel=1e-12)


def test_table3_o2_model_against_the_experiments():
    for temperature_C, r2, rms in ((25, 0.579, 0.243), (75, 0.644, 0.174),
                                   (100, 0.746, 0.150)):
        experiment = pkg.read_fig10_o2(temperature_C)['experiment']
        sim = simulate_at_times(load_o2_model(T=temperature_C + 273.15),
                                experiment['time'])
        assert pooled_r2(experiment, sim) == pytest.approx(r2, abs=2e-3)
        assert rms_difference(experiment, sim) == pytest.approx(rms, abs=6e-4)


def test_what_makes_the_table3_o2_model_too_fast():
    # NOTES.md: in Fig. 10d the 50 C point of 3-HP -> MA lies a factor of
    # 4.3 above the regression through the other temperatures. Lowering
    # the 50 C constants of reactions 1 and 2 by that factor recovers part
    # of the misfit only (and nothing at 100 C); lowering the constant of
    # 3-HP -> D tenfold as well recovers most of the rest. Pooled R2 at
    # 25, 75 and 100 C:
    for factors, expected in (
            ({'k_1_ref': 1./4.3, 'k_2_ref': 1./4.3}, (0.818, 0.851, 0.749)),
            ({'k_1_ref': 1./4.3, 'k_2_ref': 1./4.3, 'k_5_ref': 0.1},
             (0.957, 0.982, 0.928))):
        for temperature_C, r2 in zip(TEMPERATURES_C, expected):
            km = load_o2_model(T=temperature_C + 273.15)
            for name, factor in factors.items():
                km.set_value(name, factor*km.get_value(name))
            experiment = pkg.read_fig10_o2(temperature_C)['experiment']
            sim = simulate_at_times(km, experiment['time'])
            assert pooled_r2(experiment, sim) == pytest.approx(r2, abs=2e-3)
    # 3-HP -> D of the Table 3 model: 2.3e-3 1/min at 50 C and, with the
    # activation energy of the script, 4.5e-4 1/min at 25 C.
    km = load_o2_model(T=298.15)
    assert km.get_value('k_5_ref')*km.get_value('m_cat') == pytest.approx(
        2.26e-3, rel=1e-3)
    assert km.get_value('k_5')*km.get_value('m_cat') == pytest.approx(
        4.5e-4, rel=0.01)


def test_no_activation_energy_of_3hp_to_d_fits_all_three_temperatures():
    # NOTES.md: with the 50 C constant of 3-HP -> D kept (and reactions 1
    # and 2 lowered by 4.3), a low Ea_5 fits 75 and 100 C and not 25 C, a
    # high one the reverse. Pooled R2 at 25, 75 and 100 C for Ea_5 from 0
    # to 150 kJ/mol in steps of 1: the lowest of the three is 0.79 at best,
    # between 47 and 48 kJ/mol, where the values at 25 and 100 C cross.
    models, experiments = [], []
    for temperature_C in TEMPERATURES_C:
        km = load_o2_model(T=temperature_C + 273.15)
        for name in ('k_1_ref', 'k_2_ref'):
            km.set_value(name, km.get_value(name)/4.3)
        models.append(km)
        experiments.append(pkg.read_fig10_o2(temperature_C)['experiment'])
    r2 = {}
    for Ea_5 in range(0, 150001, 1000):
        r2[Ea_5] = []
        for km, experiment in zip(models, experiments):
            km.set_value('Ea_5', float(Ea_5))
            r2[Ea_5].append(pooled_r2(
                experiment, simulate_at_times(km, experiment['time'])))
    assert r2[0] == pytest.approx([0.258, 0.985, 0.931], abs=2e-3)
    assert r2[120000] == pytest.approx([0.956, 0.100, 0.202], abs=2e-3)
    assert r2[47000] == pytest.approx([0.787, 0.879, 0.794], abs=2e-3)
    assert r2[48000] == pytest.approx([0.794, 0.874, 0.785], abs=2e-3)
    assert max(r2, key=lambda Ea_5: min(r2[Ea_5])) == 47000


def time_fraction_falls_to_half(times, X_HP3):
    """Time at which ``X_HP3`` falls through 0.5, interpolated linearly
    between the two grid times around the crossing (``inf`` if it never
    does)."""
    below = np.flatnonzero(X_HP3 < 0.5)
    if len(below) == 0:
        return np.inf
    i = below[0]
    return times[i - 1] + (times[i] - times[i - 1])*(
        X_HP3[i - 1] - 0.5)/(X_HP3[i - 1] - X_HP3[i])


def maximum(times, values):
    """(value, time) of the maximum of a curve, from the parabola through
    the largest grid value and its two neighbours."""
    i = int(np.argmax(values))
    a, b, c = np.polyfit(times[i - 1:i + 2] - times[i], values[i - 1:i + 2],
                         2)
    return c - b*b/(4.*a), times[i] - b/(2.*a)


def test_table3_o2_model_is_too_fast_at_25_75_and_100C():
    # Time at which X_HP3 falls to 0.5 and maximum of the MA fraction of
    # the Table 3 O2 model. Both are interpolated on a 0.1 min grid, so the
    # documented values (NOTES.md) are pinned to their last digit.
    half, peak = {}, {}
    times = np.arange(0., 1890.01, 0.1)
    for temperature_C in TEMPERATURES_C:
        sim = simulate_at_times(load_o2_model(T=temperature_C + 273.15),
                                times, ['X_HP3', 'X_MA'])
        half[temperature_C] = time_fraction_falls_to_half(
            times, sim['X_HP3'].to_numpy())
        peak[temperature_C] = maximum(times, sim['X_MA'].to_numpy())
    assert [half[T] for T in TEMPERATURES_C] == pytest.approx(
        [458.1, 16.04, 4.16], abs=0.02)
    for T, (value, time) in ((25, (0.5218, 1479.1)), (75, (0.5842, 51.7)),
                             (100, (0.5866, 12.8))):
        assert peak[T][0] == pytest.approx(value, abs=1e-4)
        assert peak[T][1] == pytest.approx(time, abs=0.06)
    # 25 C, last sample (1890 min): 0.876 of the 3-HP is left; Table 3: 0.051
    experiment = pkg.read_fig10_o2(25)['experiment']
    km = load_o2_model(T=298.15)
    sim = simulate_at_times(km, experiment['time'])
    assert experiment['X_HP3'].iloc[-1] == 0.87637
    assert sim['X_HP3'].iloc[-1] == pytest.approx(0.051, abs=1e-3)
    assert sim['X_MA'].iloc[-1] == pytest.approx(0.51, abs=5e-3)
    # first-order constant of 3-HP loss (R1 + R5), 1/min: model and measured
    loss = (km.get_value('k_1')*km.get_value('m_cat')*km.get_value('x_O2_ref')
            + km.get_value('k_5')*km.get_value('m_cat'))
    measured = -np.log(experiment['X_HP3'].iloc[-1])/1890.
    assert (loss, measured) == pytest.approx((1.57e-3, 7.0e-5), rel=0.01)
    # Arrhenius factor of R1 and R2 from 50 to 75 C
    km = load_o2_model(T=348.15)
    assert km.get_value('k_1')/km.get_value('k_1_ref') == pytest.approx(
        4.94, abs=0.01)
    # 100 C: final acetic acid fraction 0.18; measured 0.09 to 0.11
    sim = simulate_at_times(load_o2_model(T=373.15), [1560.])
    assert sim['X_AAc'][0] == pytest.approx(0.182, abs=1e-3)
    late = pkg.read_fig10_o2(100)['experiment'].iloc[-3:]
    assert late['X_AAc'].between(0.09, 0.11).all()
