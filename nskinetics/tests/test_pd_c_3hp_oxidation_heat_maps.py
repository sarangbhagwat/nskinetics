# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests of the simplified reversible O2 model and of the temperature-time
heat maps (Fig. 11 of Kim et al., Appl. Catal. B: Environ. Energy 386 (2026)
126403) in ``nskinetics.models.pd_c_3hp_oxidation``.

The authors deposited the script that draws Fig. 11 but neither the heat-map
matrices nor the fitted constants, so the reference is the published figure
itself: the contour lines of the two main panels, digitized and shipped with
the package, and the two maxima quoted in the text. The model is checked in
two steps:

* its rate laws, Arrhenius factors and initial conditions against an
  independent evaluation of the equations of the authors' script, written
  out again in this module; and
* its heat maps against the digitized contour lines: an rms offset of at
  most 0.45 px per panel is accepted, where 1 px of the published raster is
  3.07 min or 0.148 C and the digitization is good to about 0.3 px (the
  measured offsets are 0.17 and 0.27 px). Negative controls show what this
  limit excludes: the constant set of the other panel, the starting vector
  printed in the script, and the Table 3 O2 model. The red dashed line of
  the figure (the row maxima) is compared with the digitized one as well:
  it has to sit in the published grid column on average.

The last section pins the numbers ``NOTES.md`` reports for the Table 3 O2
model against the published maps.

A heat map takes about 0.2 s, so this module is not marked slow. It imports
pytest at top level and is therefore not registered in
``nskinetics/tests/__init__.py``.
"""
import os
import re

import numpy as np
import pytest

import nskinetics as nsk
import nskinetics.models.pd_c_3hp_oxidation as pkg
from nskinetics.models.pd_c_3hp_oxidation import (
    FIG11_T_C, FIG11_TIMES, contour_offsets, contour_shifts,
    heat_map_maximum, heat_map_ridge, load_h2o2_model, load_o2_model,
    load_o2_reversible_model, load_o2_simplified_reversible_model, ode23s,
    simulate_at_times, simulate_heat_maps)

_SPECIES = ['HP3', 'MA', 'AAc', 'OXA', 'D']
_CONCENTRATIONS = ['[%s]' % s for s in _SPECIES]
_NAMES = (['k_%d_ref' % i for i in range(1, 6)]
          + ['K_%d_ref' % i for i in range(6, 11)])
C0 = 0.06                       # M
X_3HP_FITTING = 0.000339        # 3-HP mole fraction of the script's fits

# The documented constant sets at 50 C (NOTES.md), in the script's order.
FIG11B = (2.139999e-3, 5.199984e-5, 1.435e-5, 2.3040e-3, 1.062501e-3,
          4.7840, 0.5918, 3.5116, 81.6870, 52.0472)
FIG11A = (2.139999e-3, 5.199984e-5, 1.4912e-5, 2.2397e-3, 1.1110e-3,
          4.7840, 0.5918, 3.5116, 81.6870, 23.28)
K0_PRINTED = (0.0021, 0.0001, 0.0000, 0.0023, 0.0010,
              4.7840, 0.5918, 3.5116, 81.6870, 52.0472)

# Constants of APCATB_D_25_08770_heatmap.m.
P_OVER_F0 = 3.5/1.01325
R_GAS, T_REF = 8.314, 323.15
EA = np.array([59781.2, 91594.6, 10000., 51788.73, 59549.8]*2)
#              r1  r2  r3  r4  r5
NU = np.array([[-1, 0, 0, -1, 0],        # 3HP
               [1, -1, -1, 0, -1],       # MA
               [0, 1, 0, 0, 0],          # AA
               [0, 0, 1, 0, 0],          # OXA
               [0, 0, 0, 1, 1],          # D
               [1, 0, 0, 0, 0],          # H2O
               [0, 0, 1, 0, 0],          # H2
               [-1, 0, -1, 0, 0],        # O2
               [0, 1, 1, 0, 0]])         # CO2


def script_rates(k, y):
    """The function ``rates`` of the authors' script, for its 9-state y."""
    x3HP, xMA, xAA, xOXA, xD = y[:5]
    xH2O = 1. - sum(y[:5])
    yH2, yCO2 = y[6], y[8]
    yO2 = 1. - yH2 - yCO2
    return np.array([
        k[0]*(x3HP*yO2*P_OVER_F0 - xMA*xH2O/k[5]),
        k[1]*(xMA - xAA*yCO2*P_OVER_F0/k[6]),
        k[2]*(xMA*yO2*P_OVER_F0 - xOXA*yH2*yCO2*P_OVER_F0**2/k[7]),
        k[3]*(x3HP - xD/k[8]),
        k[4]*(xMA - xD/k[9])])


def script_k_adj(pars, T):
    """``k_adj`` of the script: every parameter scaled, the K's included."""
    return np.asarray(pars)*np.exp(-EA/R_GAS*(1./T - 1./T_REF))


def contour_scores(maps, digitized, panel, max_offset=30.):
    """Rms offset (px) over all contour samples, and (rms, mean) per branch.

    As for the numbers of NOTES.md, a far-off sample counts as 30 px unless
    ``max_offset`` is None."""
    T_C, times = maps[:2]
    matrix = maps[2] if panel == 'a' else maps[3]
    per_branch, offsets = {}, []
    for (level, branch), samples in digitized['contours'].groupby(
            ['level', 'branch']):
        d = contour_offsets(T_C, times, matrix, level, samples['time'],
                            samples['T_C'], max_offset=max_offset)
        per_branch[level, branch] = (np.sqrt(np.mean(d**2)), d.mean())
        offsets.append(d)
    return np.sqrt(np.mean(np.concatenate(offsets)**2)), per_branch


@pytest.fixture(scope='module')
def digitized():
    return {'a': pkg.read_fig11a_selectivity_digitized(),
            'b': pkg.read_fig11b_yield_digitized()}


@pytest.fixture(scope='module')
def maps():
    """Heat maps of the three constant sets, run as the script runs them."""
    return {name: simulate_heat_maps(
                load_o2_simplified_reversible_model(parameter_set=name),
                time_offset=pkg.FIG11_TIME_OFFSET)
            for name in ('fig11a', 'fig11b', 'k0_printed')}


@pytest.fixture(scope='module')
def table3_maps():
    return pkg.table3_o2_heat_maps()


# --- 1. loader and constant sets ---------------------------------------------

def test_loader_returns_a_configured_kinetic_model():
    km = load_o2_simplified_reversible_model()
    assert isinstance(km, nsk.KineticModel)
    assert km._units == {'time': 'min', 'conc': 'M'}
    km.validate_units()
    assert list(km._te.getFloatingSpeciesIds()) == _SPECIES
    assert km._te.integrator.absolute_tolerance == 1e-10
    assert km._te.integrator.relative_tolerance == 1e-9
    assert km.get_value('[HP3]') == pytest.approx(C0)
    for name in ('[MA]', '[AAc]', '[OXA]', '[D]', 'yH2', 'yCO2'):
        assert km.get_value(name) == 0.
    assert km.get_value('yO2') == 1.
    assert km.get_value('T') == km.get_value('T_ref') == T_REF
    assert km.get_value('PF') == pytest.approx(P_OVER_F0, rel=1e-12)
    assert load_o2_simplified_reversible_model() is not km
    km.simulate(0., 1670., 300, reset=True)
    assert len(km.results_df) == 300


def test_constant_sets_are_the_documented_ones():
    sets = pkg.SIMPLIFIED_REVERSIBLE_PARAMETER_SETS
    assert list(sets) == ['fig11b', 'fig11a', 'k0_printed']
    for name, values in (('fig11b', FIG11B), ('fig11a', FIG11A),
                         ('k0_printed', K0_PRINTED)):
        assert sets[name] == dict(zip(_NAMES, values)), name
        km = load_o2_simplified_reversible_model(parameter_set=name)
        assert [km.get_value(n) for n in _NAMES] == list(values), name
    # The Fig. 11a set differs from the Fig. 11b set in four constants only.
    assert [n for n in _NAMES if sets['fig11a'][n] != sets['fig11b'][n]] == [
        'k_3_ref', 'k_4_ref', 'k_5_ref', 'K_10_ref']


@pytest.mark.parametrize('parameter_set', ['table3', 'fig11c', 'FIG11A',
                                           None, ['fig11a']])
def test_loader_rejects_an_unknown_parameter_set(parameter_set):
    with pytest.raises(ValueError, match='unknown parameter set') as error:
        load_o2_simplified_reversible_model(parameter_set=parameter_set)
    assert "'fig11b', 'fig11a', 'k0_printed'" in str(error.value)


@pytest.mark.parametrize('T', [75., 50, 0., -5., 249.9, 450.1,
                               float('nan')])
def test_loader_rejects_a_temperature_that_is_not_in_kelvin(T):
    with pytest.raises(ValueError, match='T is in kelvin') as error:
        load_o2_simplified_reversible_model(T=T)
    assert '348.15' in str(error.value)


def test_antimony_file_holds_the_fig11b_set_and_the_heat_map_convention():
    as_shipped = pkg.pd_c_3hp_oxidation._load(
        pkg.o2_simplified_reversible_antimony_filename)
    assert [as_shipped.get_value(n) for n in _NAMES] == list(FIG11B)
    assert as_shipped.get_value('x_per_M') == pytest.approx(1./C0, rel=1e-14)
    activation = [as_shipped.get_value('Ea_%d' % i) for i in range(1, 6)]
    assert activation == list(EA[:5])
    assert as_shipped.get_value('R_gas') == R_GAS


def test_fig11b_rate_constants_come_from_the_deposited_arrhenius_table():
    points = pkg.read_fig10d_arrhenius()['points']
    row_50C = points[points['inv_T'] == 0.00309].iloc[0]
    fig11b = pkg.SIMPLIFIED_REVERSIBLE_PARAMETER_SETS['fig11b']
    for name in ('k_1', 'k_2', 'k_5'):
        assert fig11b[name + '_ref'] == pytest.approx(
            np.exp(row_50C['ln_' + name]), rel=1e-6), name


def test_initial_condition_selects_the_mole_fraction_scale():
    assert pkg.SIMPLIFIED_REVERSIBLE_INITIAL_CONDITIONS == {
        'heatmap': 1., 'fitting': X_3HP_FITTING}
    heatmap = load_o2_simplified_reversible_model()
    fitting = load_o2_simplified_reversible_model(initial_condition='fitting')
    assert heatmap.get_value('x_HP3') == pytest.approx(1., rel=1e-14)
    assert heatmap.get_value('x_H2O') == pytest.approx(0., abs=1e-14)
    assert fitting.get_value('x_HP3') == pytest.approx(X_3HP_FITTING,
                                                       rel=1e-12)
    assert fitting.get_value('x_H2O') == pytest.approx(1. - X_3HP_FITTING,
                                                       rel=1e-12)
    for km in (heatmap, fitting):
        assert km.get_value('[HP3]') == pytest.approx(C0)
        assert km.get_value('X_HP3') == pytest.approx(1.)


def test_initial_condition_aliases():
    # 'unscaled' and 'dilute' are other names of 'heatmap' and 'fitting'
    for alias, name in (('unscaled', 'heatmap'), ('dilute', 'fitting')):
        aliased = load_o2_simplified_reversible_model(initial_condition=alias)
        named = load_o2_simplified_reversible_model(initial_condition=name)
        assert aliased.get_value('x_per_M') == named.get_value('x_per_M')
    for unknown in ('concentrated', 'dilut', None, 1.):
        with pytest.raises(ValueError,
                           match='unknown initial condition') as error:
            load_o2_simplified_reversible_model(initial_condition=unknown)
        # the message lists the names and their aliases
        for name in ('heatmap', 'unscaled', 'fitting', 'dilute'):
            assert repr(name) in str(error.value)


# --- 2. rate laws against the equations of the authors' script ----------------

@pytest.mark.parametrize('initial_condition, parameter_set, T', [
    ('heatmap', 'fig11b', 357.15), ('heatmap', 'fig11a', 303.15),
    ('fitting', 'fig11b', 323.15), ('fitting', 'fig11a', 373.15)])
def test_rate_laws_are_those_of_the_script(initial_condition, parameter_set,
                                           T):
    km = load_o2_simplified_reversible_model(
        T=T, parameter_set=parameter_set, initial_condition=initial_condition)
    x_per_M = km.get_value('x_per_M')
    # A state off any trajectory, so that every term of every rate counts:
    # the liquid fractions do not add up to their initial sum (x_H2O is not
    # 0 in the heat-map convention) and the gas fractions are not small.
    concentrations = [0.021, 0.014, 0.004, 0.003, 0.008]
    yH2, yCO2 = 0.07, 0.11
    for name, value in zip(_CONCENTRATIONS, concentrations):
        km.set_value(name, value)
    km.set_value('yH2', yH2)
    km.set_value('yCO2', yCO2)

    x = x_per_M*np.array(concentrations)
    y = np.concatenate((x, [1. - x.sum(), yH2, 1. - yH2 - yCO2, yCO2]))
    pars = pkg.SIMPLIFIED_REVERSIBLE_PARAMETER_SETS[parameter_set].values()
    r = script_rates(script_k_adj(list(pars), T), y)
    dydt = NU @ r
    assert (np.abs(r) > 0.).all()

    assert km.get_value('x_H2O') == pytest.approx(y[5], rel=1e-12, abs=1e-15)
    assert km.get_value('yO2') == pytest.approx(y[7], rel=1e-12)
    for i in range(5):
        # the script's rate (1/min), and the same rate as a reaction in M/min
        assert km.get_value('r_%d' % (i + 1)) == pytest.approx(
            r[i], rel=1e-11), i
        assert km.get_value('R%d' % (i + 1)) == pytest.approx(
            r[i]/x_per_M, rel=1e-11), i
        assert km.get_value(_CONCENTRATIONS[i] + "'")*x_per_M == pytest.approx(
            dydt[i], rel=1e-10), i
    assert km.get_value("yH2'") == pytest.approx(dydt[6], rel=1e-11)
    assert km.get_value("yCO2'") == pytest.approx(dydt[8], rel=1e-11)


@pytest.mark.parametrize('T', [303.15, 348.15, 373.15])
def test_arrhenius_factors_scale_all_ten_parameters(T):
    km = load_o2_simplified_reversible_model(T=T)
    assert km.get_value('T') == T
    adjusted = script_k_adj(FIG11B, T)
    names = (['k_%d' % i for i in range(1, 6)]
             + ['K_%d' % i for i in range(6, 11)])
    for name, value in zip(names, adjusted):
        assert km.get_value(name) == pytest.approx(value, rel=1e-12), name
    # Step i and its equilibrium constant share one factor, so the ratio
    # k_i/K_(i+5) does not depend on temperature.
    for i in range(1, 6):
        ratio = km.get_value('k_%d' % i)/km.get_value('K_%d' % (i + 5))
        assert ratio == pytest.approx(FIG11B[i - 1]/FIG11B[i + 4], rel=1e-12)
    km.reset()                      # restores the state, keeps the temperature
    assert km.get_value('T') == T


def test_reference_temperature_gives_the_reference_constants():
    km = load_o2_simplified_reversible_model()
    for name in _NAMES:
        assert km.get_value(name[:-4]) == km.get_value(name), name


# --- 3. what the two initial conditions imply --------------------------------

def test_heat_map_convention_has_no_water_and_couples_gas_to_liquid():
    km = load_o2_simplified_reversible_model(T=368.15)      # 95 C
    selections = _CONCENTRATIONS + ['x_H2O', 'yH2', 'yCO2', 'yO2', 'x_AAc',
                                    'x_OXA', 'r_1', 'k_1', 'x_HP3', 'PF']
    sim = simulate_at_times(km, FIG11_TIMES - 10., selections)
    assert np.allclose(sim[_CONCENTRATIONS].sum(axis=1), C0, rtol=1e-9)
    assert (sim[_CONCENTRATIONS].to_numpy() > -1e-12).all()
    # No water: the reverse term of R1 is absent and R1 is irreversible.
    assert np.abs(sim['x_H2O']).max() < 1e-9
    assert np.allclose(sim['r_1'], sim['k_1']*sim['x_HP3']*sim['yO2']
                       * sim['PF'], rtol=1e-8, atol=1e-15)
    # The gas fractions are liquid fractions: yO2 = 1 - x_AAc - 2*x_OXA.
    assert np.allclose(sim['yH2'], sim['x_OXA'], atol=1e-12)
    assert np.allclose(sim['yCO2'], sim['x_AAc'] + sim['x_OXA'], atol=1e-12)
    assert np.allclose(sim['yO2'], 1. - sim['x_AAc'] - 2.*sim['x_OXA'],
                       atol=1e-12)
    assert sim['yO2'].min() == pytest.approx(0.850, abs=0.002)


def test_fitting_convention_keeps_water_and_oxygen_at_one():
    km = load_o2_simplified_reversible_model(initial_condition='fitting')
    sim = simulate_at_times(km, np.linspace(0., 1680., 57),
                            _CONCENTRATIONS + ['x_H2O', 'yO2'])
    assert np.allclose(sim[_CONCENTRATIONS].sum(axis=1), C0, rtol=1e-9)
    assert np.allclose(sim['x_H2O'], 1. - X_3HP_FITTING, rtol=1e-12)
    assert sim['yO2'].min() > 1. - 1e-4


def test_the_two_conventions_give_different_curves_at_50C():
    times = np.arange(0., 1680.1, 10.)
    X_MA = [simulate_at_times(load_o2_simplified_reversible_model(
                initial_condition=name), times, ['X_MA'])['X_MA']
            for name in ('heatmap', 'fitting')]
    # peak MA fraction 0.571 (heat-map convention) against 0.549 (fitting)
    assert X_MA[0].max() == pytest.approx(0.571, abs=0.002)
    assert X_MA[1].max() == pytest.approx(0.549, abs=0.002)
    assert np.abs(X_MA[0] - X_MA[1]).max() == pytest.approx(0.029, abs=0.003)


def test_pressure_dependence_is_that_of_the_rate_form():
    # The time courses behind the models (Figs. 9 and 10) are all at 3 bar
    # of O2 and the script uses 3.5 bar, so the dependence on P is that of
    # the script's rate form, which was fitted to no other pressure. It is
    # strong (NOTES.md): X_HP3 at 60 min and X_MA at 240 min, 50 C, for 1,
    # 3.5 (the default of the script) and 10 bar.
    for P, X_HP3, X_MA in ((1., 0.767, 0.268), (3.5, 0.559, 0.571),
                           (10., 0.246, 0.703)):
        km = load_o2_simplified_reversible_model()
        assert km.get_value('P') == 3.5
        km.set_value('P', P)
        sim = simulate_at_times(km, [60., 240.], ['X_HP3', 'X_MA'])
        assert (sim['X_HP3'][0], sim['X_MA'][1]) == pytest.approx(
            (X_HP3, X_MA), abs=1e-3)


# The paper's 8 h screening experiments at 50 C under 1, 3, 5, 10 and 20 bar
# of O2 (Fig. 3a), as deposited in 3-1_O2_experiments.csv: 3-HP conversion
# and MA yield in %. The text gives 72.8 % for 3 bar.
FIG3A_PRESSURES = (1., 3., 5., 10., 20.)                    # bar
FIG3A_CONVERSION = (21.66, 72.77, 63.16, 25.01, 25.87)
FIG3A_MA_YIELD = (2.654, 41.055, 20.531, 1.259, 0.589)


@pytest.mark.parametrize('initial_condition, conversion, X_MA', [
    ('heatmap', (0.8770, 0.9822, 0.9967, 0.9995, 0.9997),
     (0.2998, 0.4758, 0.5140, 0.5271, 0.5116)),
    ('fitting', (0.8551, 0.9588, 0.9786, 0.9892, 0.9946),
     (0.2726, 0.4451, 0.4893, 0.5110, 0.4973))])
def test_pressure_dependence_is_contradicted_by_the_screening_experiments(
        initial_condition, conversion, X_MA):
    # The measured conversion is highest at 3 bar and falls to a third of
    # that at 10 and 20 bar (the paper: "likely due to catalyst poisoning by
    # excess O2").
    # The model, default constant set, 50 C, 480 min: the conversion only
    # rises with P, and the MA fraction stays near 0.5 where the measured
    # yield is about 1 %. P must not be used to explore pressure.
    assert np.argmax(FIG3A_CONVERSION) == 1
    assert max(FIG3A_CONVERSION[3:]) < FIG3A_CONVERSION[1]/2.8
    modelled = []
    for P in FIG3A_PRESSURES:
        km = load_o2_simplified_reversible_model(
            initial_condition=initial_condition)
        km.set_value('P', P)
        sim = simulate_at_times(km, [480.], ['X_HP3', 'X_MA'])
        modelled.append((1. - sim['X_HP3'][0], sim['X_MA'][0]))
    modelled = np.array(modelled)
    assert modelled[:, 0] == pytest.approx(conversion, abs=2e-4)
    assert modelled[:, 1] == pytest.approx(X_MA, abs=2e-4)
    assert np.all(np.diff(modelled[:, 0]) > 0.)
    assert modelled[3:, 0].min() > 0.98 and modelled[3:, 1].min() > 0.49
    assert max(FIG3A_MA_YIELD[3:]) < 1.3


def test_model_at_the_temperatures_of_the_screening_experiments():
    # Fig. 3c (3-3_O2_experiments.csv of the deposit): 8 h under 3 bar O2 at
    # 25, 50, 75 and 100 C gave conversions of 8.5, 72.8, 95.5 and 84.6 %
    # and MA yields of 1.4, 41.1, 35.7 and 2.1 %. The model of the heat
    # maps, as loaded by default, at 480 min of reaction (NOTES.md):
    modelled = []
    for T_C in (25., 50., 75., 100.):
        sim = simulate_at_times(load_o2_simplified_reversible_model(
            T=T_C + 273.15), [480.], ['X_HP3', 'X_MA'])
        modelled.append((1. - sim['X_HP3'][0], sim['X_MA'][0]))
    assert np.array(modelled) == pytest.approx(np.array(
        [(0.5348, 0.3630), (0.9887, 0.4907), (0.9994, 0.0598),
         (0.9998, 0.0017)]), abs=2e-4)


def test_equilibrium_constants_of_steps_1_and_3_do_not_shape_the_maps():
    # K_6 only enters with x_H2O, which is 0; the reverse of R3 is of third
    # order in the by-products.
    T_C = [50., 95.]
    reference = simulate_heat_maps(load_o2_simplified_reversible_model(),
                                   T_C=T_C, time_offset=10.)
    for name, factor, limit in (('K_6_ref', 10., 1e-8), ('K_8_ref', 1.1, 1e-6)):
        km = load_o2_simplified_reversible_model()
        km.set_value(name, factor*km.get_value(name))
        changed = simulate_heat_maps(km, T_C=T_C, time_offset=10.)
        for i in (2, 3):
            assert np.abs(changed[i] - reference[i]).max() < limit, name


# --- 4. heat-map helper -------------------------------------------------------

def test_heat_map_grid_is_the_grid_of_the_script(maps):
    T_C, times, selectivity, MA_yield = maps['fig11b']
    assert T_C is FIG11_T_C and times is FIG11_TIMES
    assert np.array_equal(T_C, np.arange(30., 100.5, 1.))
    assert np.allclose(times, np.linspace(10., 1680., 300), rtol=0, atol=0)
    assert selectivity.shape == MA_yield.shape == (71, 300)
    # The feed sits at the first plotted time, 10 min.
    assert np.abs(selectivity[:, 0]).max() == 0.
    assert np.abs(MA_yield[:, 0]).max() == 0.
    assert ((selectivity >= 0.) & (selectivity < 1.)).all()
    assert (MA_yield <= selectivity + 1e-15).all()


def test_simulate_heat_maps_follows_the_model_and_restores_it():
    km = load_o2_simplified_reversible_model(T=333.15)
    T_C, times, selectivity, MA_yield = simulate_heat_maps(
        km, T_C=[50., 84.], times=[10., 49.1, 400.], time_offset=10.)
    assert km.get_value('T') == 333.15
    assert km.get_value('[HP3]') == pytest.approx(C0)
    direct = simulate_at_times(
        load_o2_simplified_reversible_model(T=84. + 273.15), [39.1, 390.],
        ['X_HP3', 'X_MA', 'fig11a_selectivity', 'fig11b_yield'])
    assert np.allclose(selectivity[1, 1:], direct['X_MA'], rtol=1e-9)
    assert np.allclose(selectivity[1, 1:], direct['fig11a_selectivity'],
                       rtol=1e-9)
    assert np.allclose(MA_yield[1, 1:], direct['fig11b_yield'], rtol=1e-9)
    assert np.allclose(MA_yield[1, 1:], (1. - direct['X_HP3'])*direct['X_MA'],
                       rtol=1e-9)
    # Without the offset the same columns are 10 min further along.
    later = simulate_heat_maps(km, T_C=[84.], times=[39.1, 390.])
    assert np.allclose(later[2][0], selectivity[1, 1:], rtol=1e-9)


@pytest.mark.parametrize('T_C', [[323.15], [50., 348.15], [-60.], [201.],
                                 [float('nan')]])
def test_simulate_heat_maps_rejects_temperatures_that_are_not_in_celsius(T_C):
    # 323.15 would otherwise be taken as 323.15 C
    km = load_o2_simplified_reversible_model()
    with pytest.raises(ValueError, match='degrees Celsius') as error:
        simulate_heat_maps(km, T_C=T_C, times=[30.])
    assert 'for 348.15 K pass 75' in str(error.value)
    assert km.get_value('T') == T_REF       # the model is left untouched


@pytest.mark.parametrize('load', [load_h2o2_model,
                                  load_o2_reversible_model])
def test_simulate_heat_maps_refuses_a_model_without_a_temperature(load):
    # these two models have no parameter T; without the check roadrunner
    # raises a RuntimeError that names neither the function nor the models
    with pytest.raises(ValueError, match='temperature parameter T') as error:
        simulate_heat_maps(load(), T_C=[50.], times=[30.])
    assert 'load_o2_model or load_o2_simplified_reversible_model' in str(
        error.value)


def test_simulate_heat_maps_accepts_minus_50_to_200_C():
    km = load_o2_model()
    T_C = simulate_heat_maps(km, T_C=[-50., 200.], times=[30.])[0]
    assert list(T_C) == [-50., 200.]


def test_fig11_heat_maps_uses_one_constant_set_per_panel(maps):
    T_C, times, selectivity, MA_yield = pkg.fig11_heat_maps()
    assert np.array_equal(selectivity, maps['fig11a'][2])
    assert np.array_equal(MA_yield, maps['fig11b'][3])
    assert not np.allclose(selectivity, maps['fig11b'][2], atol=1e-3)
    for panel, matrix in (('a', selectivity), ('b', MA_yield)):
        panel_T_C, panel_times, panel_matrix = pkg.fig11_panel_heat_map(panel)
        assert panel_T_C is FIG11_T_C and panel_times is FIG11_TIMES
        assert np.array_equal(panel_matrix, matrix)
    with pytest.raises(ValueError, match="'a' or 'b'"):
        pkg.fig11_panel_heat_map('c')


def test_contour_measures_on_a_plane():
    # value = t/100 + (T - 30)/10: the 3.0 contour is the line t = 300 - 10*(T - 30)
    T_C, times = np.arange(30., 61.), np.linspace(0., 600., 61)
    matrix = times[None, :]/100. + (T_C[:, None] - 30.)/10.
    t, T = np.array([210., 100., 700.]), np.array([40.5, 45.2, 70.])
    value = t/100. + (T - 30.)/10.
    per_px = np.hypot(0.02, 0.05)       # with 2 min and 0.5 C per pixel
    offsets = contour_offsets(T_C, times, matrix, 3., t, T, min_per_px=2.,
                              C_per_px=0.5)
    assert offsets == pytest.approx((value[:2] - 3.)/per_px)   # 3rd: outside
    shifts = contour_shifts(T_C, times, matrix, 3., 'early', t, T)
    assert shifts[:2] == pytest.approx(300. - 10.*(T[:2] - 30.) - t[:2])
    assert np.isnan(shifts[2])
    shifts = contour_shifts(T_C, times, matrix, 3., 'early', t, T,
                            along='T_C')
    assert shifts[:2] == pytest.approx(60. - t[:2]/10. - T[:2])
    # the plane never falls through the level: no 'late' crossing
    assert np.isnan(contour_shifts(T_C, times, matrix, 3., 'late', t, T)).all()
    with pytest.raises(ValueError, match="'time' or 'T_C'"):
        contour_shifts(T_C, times, matrix, 3., 'early', t, T, along='T')
    # offsets can be clipped, so that far-off samples do not dominate an rms
    clipped = contour_offsets(T_C, times, matrix, 3., t, T, min_per_px=2.,
                              C_per_px=0.5, max_offset=5.)
    assert clipped == pytest.approx(np.clip(offsets, -5., 5.))
    assert np.abs(offsets).max() > 5.


def test_contour_shifts_uses_the_crossing_nearest_to_each_sample():
    # Every row rises through 0.5 twice, at t = 100*pi/6 = 52.36 and at
    # t = 100*(2*pi + pi/6) = 680.68 min, and falls through it in between,
    # at t = 100*5*pi/6 = 261.80 min, and once more at 890.12 min.
    T_C, times = np.arange(30., 41.), np.linspace(0., 900., 9001)
    matrix = np.tile(np.sin(times/100.), (len(T_C), 1))
    rising = 100.*np.array([np.pi/6., 2.*np.pi + np.pi/6.])
    falling = 100.*np.array([5.*np.pi/6., 2.*np.pi + 5.*np.pi/6.])
    t, T = np.array([60., 300., 400., 670., 880.]), np.full(5, 35.3)
    shifts = contour_shifts(T_C, times, matrix, 0.5, 'early', t, T)
    assert shifts == pytest.approx(
        [rising[0] - 60., rising[0] - 300., rising[1] - 400.,
         rising[1] - 670., rising[1] - 880.], abs=1e-3)
    # the sample at 670 min is 10.7 min from the second crossing, not
    # 618 min from the first
    assert abs(shifts[3]) < 11.
    shifts = contour_shifts(T_C, times, matrix, 0.5, 'late', t, T)
    assert shifts == pytest.approx(
        [falling[0] - 60., falling[0] - 300., falling[0] - 400.,
         falling[1] - 670., falling[1] - 880.], abs=1e-3)


def test_raster_scale_is_the_calibration_of_panel_a():
    # The pixel scale behind every offset in px is the axis calibration of
    # the published raster recorded in reference_data/README.md.
    with open(os.path.join(pkg.reference_data_dir, 'README.md'),
              encoding='utf-8') as f:
        readme = f.read()
    (row,) = re.findall(r'^\| a \| t = ([\d.]+) col - [\d.]+ min '
                        r'\| T = -([\d.]+) row \+ [\d.]+ C \|$', readme,
                        flags=re.MULTILINE)
    assert pkg.FIG11_MIN_PER_PX == float(row[0]) == 3.0723
    assert pkg.FIG11_C_PER_PX == float(row[1]) == 0.1476
    # the main panel (half a cell beyond the grid on each side) is then
    # 545 x 481 px of the 1604 x 596 px raster
    assert (1682.8 - 7.2)/pkg.FIG11_MIN_PER_PX == pytest.approx(545.4,
                                                                abs=0.1)
    assert 71./pkg.FIG11_C_PER_PX == pytest.approx(481.0, abs=0.1)


# --- 5. the published Fig. 11 -------------------------------------------------

@pytest.mark.parametrize('panel, parameter_set', [('a', 'fig11a'),
                                                  ('b', 'fig11b')])
def test_heat_maps_reproduce_the_published_contours(maps, digitized, panel,
                                                    parameter_set):
    rms, per_branch = contour_scores(maps[parameter_set], digitized[panel],
                                     panel)
    assert rms <= 0.45, rms                     # 0.17 (a) and 0.27 (b)
    assert len(per_branch) == 6
    for branch, (_, mean) in per_branch.items():
        assert abs(mean) <= 0.4, (branch, mean)


@pytest.mark.parametrize('panel, parameter_set, at_least', [
    ('a', 'fig11b', 1.0),           # the other panel's constants: 1.07 px
    ('b', 'fig11a', 1.0),           # 1.45 px
    ('a', 'k0_printed', 4.0),       # the vector printed in the script: 4.6
    ('b', 'k0_printed', 6.0)])      # 7.1
def test_other_constant_sets_miss_the_published_contours(
        maps, digitized, panel, parameter_set, at_least):
    rms, _ = contour_scores(maps[parameter_set], digitized[panel], panel)
    assert rms >= at_least, rms


def test_maxima_quoted_in_the_text(maps):
    T_C, times, selectivity, _ = maps['fig11a']
    row, column = list(T_C).index(95.), 3
    assert times[column] == pytest.approx(26.76, abs=0.005)
    # "59.1 % (95 C, 27 min)"
    assert selectivity[row, column] == pytest.approx(0.5906, abs=0.003)
    assert selectivity[row, column] > selectivity[T_C <= 80.].max()

    _, _, _, MA_yield = maps['fig11b']
    row, column = list(T_C).index(84.), 7
    assert times[column] == pytest.approx(49.10, abs=0.005)
    # "55.2 % (84 C, 49 min)": the cell is reproduced, its value is 0.548.
    # The crest is flat, so the cell is not required to be the largest one.
    assert MA_yield[row, column] == pytest.approx(0.5483, abs=0.003)
    assert MA_yield.max() - MA_yield[row, column] < 1e-3
    assert MA_yield.max() < 0.552


def test_observed_maxima_quoted_in_the_text(maps):
    # "experimentally observed maximum was 50.5 %": the quantity of panel b,
    # (1 - X_HP3)*X_MA, of the deposited 50 C experiment is largest at
    # 360 min and equals that number. The paper does not say how it was
    # computed; X_MA alone is largest at 240 min (0.548, pinned in
    # test_pd_c_3hp_oxidation.py) and never 0.505.
    experiment = pkg.read_fig9b_o2()['experiment']
    time = experiment['time'].to_numpy()
    X_HP3 = experiment['X_HP3'].to_numpy()
    X_MA = experiment['X_MA'].to_numpy()
    product = (1. - X_HP3)*X_MA
    assert time[product.argmax()] == 360.
    assert product.max() == pytest.approx(0.5051, abs=5e-5)
    assert np.sort(product)[-2] < 0.49
    assert np.abs(X_MA - 0.505).min() > 0.01
    # The companion number, an "observed maximum of 56.9 % at 50 C (purple
    # star on plot)" for the selectivity of panel a, is not identified: it
    # is not X_MA of that experiment (at most 0.548), and MA per 3-HP
    # converted, X_MA/(1 - X_HP3), is 0.621, 0.577 and 0.513 at 240, 360 and
    # 450 min. The 50 C row of the reproduced map peaks at 0.5705.
    assert X_MA.max() < 0.569 - 0.02
    late = time >= 240.
    per_converted = X_MA[late]/(1. - X_HP3[late])
    assert per_converted[:3] == pytest.approx((0.6207, 0.5774, 0.5128),
                                              abs=5e-5)
    assert np.all(np.diff(per_converted) < 0.)
    assert np.abs(per_converted - 0.569).min() > 0.008
    T_C, _, selectivity, _ = maps['fig11a']
    assert selectivity[list(T_C).index(50.)].max() == pytest.approx(
        0.5705, abs=2e-4)


def test_50C_row_maxima(maps):
    T_C, times, selectivity, _ = maps['fig11a']
    row = list(T_C).index(50.)
    assert selectivity[row].max() == pytest.approx(0.5705, abs=2e-4)
    assert heat_map_ridge(times, selectivity)[row] == pytest.approx(
        255.8, abs=6.)
    _, _, _, MA_yield = maps['fig11b']
    assert MA_yield[row].max() == pytest.approx(0.5331, abs=2e-4)
    assert heat_map_ridge(times, MA_yield)[row] == pytest.approx(322.8,
                                                                 abs=6.)


def test_heat_map_maximum_and_ridge(maps):
    T_C, times, selectivity, _ = maps['fig11a']
    value, T_max, t_max = heat_map_maximum(T_C, times, selectivity)
    assert value == selectivity.max()
    assert selectivity[list(T_C).index(T_max),
                       list(times).index(t_max)] == value
    ridge = heat_map_ridge(times, selectivity)
    assert ridge.shape == (71,)
    assert (np.diff(ridge) <= 0.).all()         # hotter rows peak earlier


def ridge_offsets(T_C, ridge, digitized):
    """Model ridge minus the digitized red dashed line, along the time axis
    at the temperature of each digitized sample, in px of the raster."""
    samples = digitized['ridge']
    return (np.interp(samples['T_C'], T_C, ridge)
            - samples['time'])/pkg.FIG11_MIN_PER_PX


@pytest.mark.parametrize('panel, parameter_set, mean, rms', [
    ('a', 'fig11a', 0.60, 0.93), ('b', 'fig11b', 0.15, 1.48)])
def test_ridge_is_the_published_red_dashed_line(maps, digitized, panel,
                                                parameter_set, mean, rms):
    # The red dashed line joins the row maxima, so it is defined to one
    # grid column: 5.585 min, 1.82 px of the published raster. It is
    # therefore held to the digitized line to half a column in the mean
    # (a line that sits one column off is 1.82 px off in the mean) and to
    # one column in the rms. The largest single offsets (2.1 and 6.5 px)
    # are at 31 to 32 C, where the row maxima are flat and the line
    # advances several columns per row.
    T_C, times = maps[parameter_set][:2]
    matrix = maps[parameter_set][2 if panel == 'a' else 3]
    column = (times[1] - times[0])/pkg.FIG11_MIN_PER_PX
    assert column == pytest.approx(1.818, abs=1e-3)
    assert len(digitized[panel]['ridge']) == {'a': 118, 'b': 226}[panel]
    offsets = ridge_offsets(T_C, heat_map_ridge(times, matrix),
                            digitized[panel])
    assert offsets.mean() == pytest.approx(mean, abs=0.1)
    assert abs(offsets.mean()) <= 0.5*column
    assert np.sqrt(np.mean(offsets**2)) == pytest.approx(rms, abs=0.1)
    assert np.sqrt(np.mean(offsets**2)) <= column
    # negative control: the same line one column later or earlier fails
    # the limit on the mean (2.42 and -1.22 px in panel a, 1.97 and
    # -1.67 px in panel b)
    for shift in (1, -1):
        shifted = times[np.argmax(matrix, axis=1) + shift]
        assert abs(ridge_offsets(T_C, shifted, digitized[panel]).mean()
                   ) > 0.5*column + 0.2


@pytest.mark.parametrize('panel, parameter_set, largest', [
    ('a', FIG11A, 1.12e-4), ('b', FIG11B, 1.36e-4)])
def test_solver_of_the_script_does_not_matter_for_the_maps(maps, panel,
                                                           parameter_set,
                                                           largest):
    # The script integrates each of the 299 intervals of the time grid with
    # ode23s at default tolerances; the maps here are converged solutions.
    # The 100 C row, repeated as the script computes it, differs by
    # 1.1e-4 (a) and 1.4e-4 (b); in the other rows the difference is
    # smaller (measured for all rows outside the tests, NOTES.md).
    T_C, times, selectivity, MA_yield = maps['fig11' + panel]
    k = script_k_adj(parameter_set, 373.15)
    y = np.array([1., 0., 0., 0., 0., 0., 0., 1., 0.])     # unit 3-HP, O2
    rows = [y]
    for interval in np.diff(times):
        y = ode23s(lambda y: NU @ script_rates(k, y), interval, y)[1][-1]
        rows.append(y)
    liquid = np.array(rows)[:, :5]
    liquid = liquid/liquid.sum(axis=1, keepdims=True)
    script, converged = (
        (liquid[:, 1], selectivity[-1]) if panel == 'a'
        else ((1. - liquid[:, 0])*liquid[:, 1], MA_yield[-1]))
    assert T_C[-1] == 100.
    assert np.abs(script - converged).max() == pytest.approx(largest,
                                                             rel=0.05)


# --- 6. Table 3 O2 model against the published maps (numbers of NOTES.md) -----

def test_table3_o2_model_misses_the_published_contours(table3_maps,
                                                       digitized):
    for panel in 'ab':
        rms, per_branch = contour_scores(table3_maps, digitized[panel], panel)
        assert rms == pytest.approx(7.6, abs=0.2), (panel, rms)
    # Offsets this large are linearized estimates and depend on the 30 px
    # clip: without it panel a is unchanged and panel b reads 8.8 px.
    assert contour_scores(table3_maps, digitized['a'], 'a',
                          max_offset=None)[0] == pytest.approx(7.6, abs=0.2)
    assert contour_scores(table3_maps, digitized['b'], 'b',
                          max_offset=None)[0] == pytest.approx(8.75, abs=0.1)
    # With the 10 min offset of the published figure it is hardly closer.
    shifted = pkg.table3_o2_heat_maps(time_offset=10.)
    assert contour_scores(shifted, digitized['a'], 'a')[0] == pytest.approx(
        6.8, abs=0.6)
    assert contour_scores(shifted, digitized['b'], 'b')[0] == pytest.approx(
        6.1, abs=0.6)


def test_table3_o2_model_late_contours_are_early(table3_maps, digitized):
    T_C, times = table3_maps[:2]
    for panel, matrix in (('a', table3_maps[2]), ('b', table3_maps[3])):
        contours = digitized[panel]['contours']
        late = contours[contours['branch'] == 'late']
        for level, samples in late.groupby('level'):
            shift = contour_shifts(T_C, times, matrix, level, 'late',
                                   samples['time'], samples['T_C'])
            relative = shift/samples['time']
            # 8 to 11 % earlier than published (median per line) ...
            assert -0.12 < np.nanmedian(relative) < -0.07, (panel, level)
            # ... which is 1.3 to 1.8 C lower at a given time
            shift = contour_shifts(T_C, times, matrix, level, 'late',
                                   samples['time'], samples['T_C'],
                                   along='T_C')
            assert -2.0 < np.nanmedian(shift) < -1.1, (panel, level)


def test_table3_o2_model_maxima(table3_maps):
    T_C, times, selectivity, MA_yield = table3_maps
    value, T_max, t_max = heat_map_maximum(T_C, times, selectivity)
    assert value == pytest.approx(0.587, abs=0.003)     # at 91 C, 21 min
    assert 86. <= T_max <= 96. and 15. <= t_max <= 30.
    value, T_max, t_max = heat_map_maximum(T_C, times, MA_yield)
    assert value == pytest.approx(0.542, abs=0.003)     # at 84 C, 38 min
    assert 79. <= T_max <= 89. and 30. <= t_max <= 46.
    row = list(T_C).index(50.)
    assert selectivity[row].max() == pytest.approx(0.566, abs=0.003)
    assert MA_yield[row].max() == pytest.approx(0.527, abs=0.003)


def test_table3_o2_heat_maps_is_the_o2_model_at_each_temperature(table3_maps):
    T_C, times, selectivity, MA_yield = table3_maps
    row = list(T_C).index(50.)
    sim = simulate_at_times(load_o2_model(), times, ['X_HP3', 'X_MA'])
    assert np.allclose(selectivity[row], sim['X_MA'], rtol=1e-9)
    assert np.allclose(MA_yield[row], (1. - sim['X_HP3'])*sim['X_MA'],
                       rtol=1e-9)


# --- 7. reference data --------------------------------------------------------

def test_digitized_fig11_readers(digitized):
    for panel, n_contour, red_star, purple_star in (
            ('a', 1254, (27.5, 94.7), (254.3, 50.1)),
            ('b', 1809, (55.9, 83.9), (335.6, 49.8))):
        tables = digitized[panel]
        contours = tables['contours']
        assert list(contours.columns) == ['level', 'branch', 'time', 'T_C']
        assert len(contours) == n_contour
        assert sorted(set(contours['level'])) == [0.1, 0.3, 0.5]
        assert set(contours['branch']) == {'early', 'late'}
        assert contours['time'].between(7., 1683.).all()
        assert contours['T_C'].between(29.5, 100.5).all()
        assert len(tables['ridge']) > 100
        stars = tables['stars']
        assert tuple(stars.loc['red_star']) == pytest.approx(red_star,
                                                             abs=0.05)
        assert tuple(stars.loc['purple_star']) == pytest.approx(purple_star,
                                                                abs=0.05)
    size = sum(os.path.getsize(os.path.join(pkg.reference_data_dir, name))
               for name in ('fig11a_selectivity_digitized.csv',
                            'fig11b_yield_digitized.csv'))
    assert size < 150e3


def test_published_stars_against_the_cells_named_in_the_text(digitized):
    # Red stars: panel a sits on the cell of the text (95 C, 26.76 min);
    # panel b is drawn 7 min to the right of its cell (84 C, 49.10 min).
    star = digitized['a']['stars'].loc['red_star']
    assert (star['time'] - FIG11_TIMES[3], star['T_C'] - 95.) == pytest.approx(
        (0.8, -0.3), abs=0.1)
    star = digitized['b']['stars'].loc['red_star']
    assert (star['time'] - FIG11_TIMES[7], star['T_C'] - 84.) == pytest.approx(
        (6.8, -0.1), abs=0.1)
