# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests of the simplified reversible O2 model, of the fit that defines two of
its constants, and of the temperature-time heat maps (Fig. 11 of Kim et al.,
Appl. Catal. B: Environ. Energy 386 (2026) 126403) in
``nskinetics.models.pd_c_3hp_oxidation``.

The authors deposited the script that draws Fig. 11 but neither the heat-map
matrices nor the fitted constants. The constants of the model are defined
from the deposit: eight are values the deposit gives, and ``k_3`` and
``k_4`` are a least-squares fit to the 50 C data table of the script. The
model is checked in three steps:

* its rate laws, Arrhenius factors and initial conditions against an
  independent evaluation of the equations of the authors' script, written
  out again in this module;
* its constants: the eight pinned ones against the deposit, and ``k_3`` and
  ``k_4`` as the result of ``fit_simplified_reversible_constants``, from
  three starting points; and
* its heat maps against the published figure: the contour lines of the two
  main panels, digitized and shipped with the package, and the two maxima
  quoted in the text. The constants are not fitted to the figure, so this
  is a comparison, and its result is pinned as measured: an rms offset of
  0.35 px in panel b and 1.11 px in panel a, where 1 px of the published
  raster is 3.07 min or 0.148 C and the digitization is good to about
  0.3 px. The vector ``k0`` printed in the script is far off both
  panels. The red dashed line of the figure (the row maxima) is compared
  with the digitized one as well: it has to sit in the published grid
  column on average.

A heat map takes about 0.2 s and the fit less than a second, so this module
is not marked slow. It imports pytest at top level and is therefore not
registered in ``nskinetics/tests/__init__.py``.
"""
import os
import re

import numpy as np
import pytest

import nskinetics as nsk
import nskinetics.models.pd_c_3hp_oxidation as pkg
from nskinetics.models.pd_c_3hp_oxidation import (
    FIG11_T_C, FIG11_TIMES, contour_offsets, contour_shifts,
    fit_simplified_reversible_constants, heat_map_maximum, heat_map_ridge,
    load_h2o2_model, load_o2_model, load_o2_simplified_reversible_model,
    ode23s, simulate_at_times, simulate_heat_maps)

_SPECIES = ['HP3', 'MA', 'AAc', 'OXA', 'D']
_CONCENTRATIONS = ['[%s]' % s for s in _SPECIES]
_NAMES = (['k_%d_ref' % i for i in range(1, 6)]
          + ['K_%d_ref' % i for i in range(6, 11)])
C0 = 0.06                       # M
X_3HP_FITTING = 0.000339        # 3-HP mole fraction of the script's fits

# The documented constant sets at 50 C (NOTES.md), in the script's order.
DEPOSIT = (2.139999e-3, 5.199984e-5, 1.4200e-5, 2.3204e-3, 1.062501e-3,
           4.7840, 0.5918, 3.5116, 81.6870, 52.0472)
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
    """Heat maps of the two constant sets, run as the script runs them."""
    return {name: simulate_heat_maps(
                load_o2_simplified_reversible_model(parameter_set=name),
                time_offset=pkg.FIG11_TIME_OFFSET)
            for name in ('deposit', 'k0_printed')}


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
    assert list(sets) == ['deposit', 'k0_printed']
    for name, values in (('deposit', DEPOSIT), ('k0_printed', K0_PRINTED)):
        assert sets[name] == dict(zip(_NAMES, values)), name
        km = load_o2_simplified_reversible_model(parameter_set=name)
        assert [km.get_value(n) for n in _NAMES] == list(values), name
    # the default of the loader is the 'deposit' set
    km = load_o2_simplified_reversible_model()
    assert [km.get_value(n) for n in _NAMES] == list(DEPOSIT)
    # K_6 to K_10 are those of the vector k0 printed in the script
    assert DEPOSIT[5:] == K0_PRINTED[5:]


@pytest.mark.parametrize('parameter_set', ['table3', 'fitted', 'Deposit',
                                           None, ['deposit']])
def test_loader_rejects_an_unknown_parameter_set(parameter_set):
    with pytest.raises(ValueError, match='unknown parameter set') as error:
        load_o2_simplified_reversible_model(parameter_set=parameter_set)
    assert "'deposit', 'k0_printed'" in str(error.value)


@pytest.mark.parametrize('T', [75., 50, 0., -5., 249.9, 450.1,
                               float('nan')])
def test_loader_rejects_a_temperature_that_is_not_in_kelvin(T):
    with pytest.raises(ValueError, match='T is in kelvin') as error:
        load_o2_simplified_reversible_model(T=T)
    assert '348.15' in str(error.value)


def test_antimony_file_holds_the_deposit_set_and_the_heat_map_convention():
    as_shipped = pkg.pd_c_3hp_oxidation._load(
        pkg.o2_simplified_reversible_antimony_filename)
    assert [as_shipped.get_value(n) for n in _NAMES] == list(DEPOSIT)
    assert as_shipped.get_value('x_per_M') == pytest.approx(1./C0, rel=1e-14)
    activation = [as_shipped.get_value('Ea_%d' % i) for i in range(1, 6)]
    assert activation == list(EA[:5])
    assert as_shipped.get_value('R_gas') == R_GAS


def test_three_rate_constants_come_from_the_deposited_arrhenius_table():
    points = pkg.read_fig10d_arrhenius()['points']
    row_50C = points[points['inv_T'] == 0.00309].iloc[0]
    deposit = pkg.SIMPLIFIED_REVERSIBLE_PARAMETER_SETS['deposit']
    for name in ('k_1', 'k_2', 'k_5'):
        assert deposit[name + '_ref'] == pytest.approx(
            np.exp(row_50C['ln_' + name]), rel=1e-6), name


def test_equilibrium_constants_are_those_printed_in_the_script():
    # k0 = [0.0021 0.0001 0.0000 0.0023 0.0010 4.7840 0.5918 3.5116 81.6870
    # 52.0472] in APCATB_D_25_08770_heatmap.m
    assert K0_PRINTED[5:] == (4.7840, 0.5918, 3.5116, 81.6870, 52.0472)
    deposit = pkg.SIMPLIFIED_REVERSIBLE_PARAMETER_SETS['deposit']
    assert [deposit['K_%d_ref' % i] for i in range(6, 11)] == list(
        K0_PRINTED[5:])


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
    ('heatmap', 'deposit', 357.15), ('heatmap', 'deposit', 303.15),
    ('fitting', 'deposit', 323.15), ('fitting', 'deposit', 373.15)])
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
    adjusted = script_k_adj(DEPOSIT, T)
    names = (['k_%d' % i for i in range(1, 6)]
             + ['K_%d' % i for i in range(6, 11)])
    for name, value in zip(names, adjusted):
        assert km.get_value(name) == pytest.approx(value, rel=1e-12), name
    # Step i and its equilibrium constant share one factor, so the ratio
    # k_i/K_(i+5) does not depend on temperature.
    for i in range(1, 6):
        ratio = km.get_value('k_%d' % i)/km.get_value('K_%d' % (i + 5))
        assert ratio == pytest.approx(DEPOSIT[i - 1]/DEPOSIT[i + 4], rel=1e-12)
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
    # peak MA fraction 0.570 (heat-map convention) against 0.549 (fitting)
    assert X_MA[0].max() == pytest.approx(0.5704, abs=0.001)
    assert X_MA[1].max() == pytest.approx(0.5485, abs=0.001)
    assert np.abs(X_MA[0] - X_MA[1]).max() == pytest.approx(0.0295,
                                                            abs=0.001)


def test_pressure_dependence_is_that_of_the_rate_form():
    # The time courses behind the models (Figs. 9 and 10) are all at 3 bar
    # of O2 and the script uses 3.5 bar, so the dependence on P is that of
    # the script's rate form, which was fitted to no other pressure. It is
    # strong (NOTES.md): X_HP3 at 60 min and X_MA at 240 min, 50 C, for 1,
    # 3.5 (the default of the script) and 10 bar.
    for P, X_HP3, X_MA in ((1., 0.767, 0.268), (3.5, 0.559, 0.570),
                           (10., 0.246, 0.702)):
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
    ('heatmap', (0.8780, 0.9823, 0.9968, 0.9995, 0.9997),
     (0.2990, 0.4750, 0.5135, 0.5270, 0.5120)),
    ('fitting', (0.8561, 0.9589, 0.9787, 0.9892, 0.9946),
     (0.2718, 0.4443, 0.4887, 0.5109, 0.4977))])
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
    # maps at 480 min of reaction, as loaded by default and with the
    # script's fitting start (NOTES.md):
    for initial_condition, expected in (
            ('heatmap', [(0.5355, 0.3628), (0.9888, 0.4900),
                         (0.9994, 0.0597), (0.9998, 0.0017)]),
            ('fitting', [(0.5038, 0.3297), (0.9668, 0.4607),
                         (0.9988, 0.0580), (0.9999, 0.0014)])):
        modelled = []
        for T_C in (25., 50., 75., 100.):
            sim = simulate_at_times(load_o2_simplified_reversible_model(
                T=T_C + 273.15, initial_condition=initial_condition),
                [480.], ['X_HP3', 'X_MA'])
            modelled.append((1. - sim['X_HP3'][0], sim['X_MA'][0]))
        assert np.array(modelled) == pytest.approx(np.array(expected),
                                                   abs=2e-4)


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
    T_C, times, selectivity, MA_yield = maps['deposit']
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


@pytest.mark.parametrize('load', [load_h2o2_model, load_o2_model])
def test_simulate_heat_maps_refuses_a_model_without_a_temperature(load):
    # the two Table 3 models have no parameter T; without the check
    # roadrunner raises a RuntimeError that names neither the function nor
    # the models
    with pytest.raises(ValueError, match='temperature parameter T') as error:
        simulate_heat_maps(load(), T_C=[50.], times=[30.])
    assert 'load_o2_simplified_reversible_model' in str(error.value)


def test_simulate_heat_maps_accepts_minus_50_to_200_C():
    km = load_o2_simplified_reversible_model()
    T_C = simulate_heat_maps(km, T_C=[-50., 200.], times=[30.])[0]
    assert list(T_C) == [-50., 200.]


def test_fig11_heat_maps_uses_the_default_constant_set(maps):
    T_C, times, selectivity, MA_yield = pkg.fig11_heat_maps()
    assert T_C is FIG11_T_C and times is FIG11_TIMES
    assert np.array_equal(selectivity, maps['deposit'][2])
    assert np.array_equal(MA_yield, maps['deposit'][3])
    assert not np.allclose(selectivity, maps['k0_printed'][2], atol=1e-3)


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

@pytest.mark.parametrize('panel, rms, largest_mean', [('a', 1.109, 0.54),
                                                      ('b', 0.345, 0.29)])
def test_heat_maps_against_the_published_contours(maps, digitized, panel,
                                                  rms, largest_mean):
    # The constants are not fitted to the figure. Panel b is matched to
    # about the reading uncertainty of the digitized lines (0.3 px); panel a
    # is not: its late 0.1 and 0.5 contours are 1.3 and 1.5 px rms off.
    measured, per_branch = contour_scores(maps['deposit'], digitized[panel],
                                          panel)
    assert measured == pytest.approx(rms, abs=0.01)
    assert len(per_branch) == 6
    assert max(abs(mean) for _, mean in per_branch.values()) == pytest.approx(
        largest_mean, abs=0.01)
    # no sample is far enough for the 30 px clip to matter
    assert contour_scores(maps['deposit'], digitized[panel], panel,
                          max_offset=None)[0] == pytest.approx(measured,
                                                               rel=1e-12)
    if panel == 'a':
        assert per_branch[0.1, 'late'][0] == pytest.approx(1.25, abs=0.01)
        assert per_branch[0.5, 'late'][0] == pytest.approx(1.48, abs=0.01)
        assert max(per_branch[level, 'early'][0]
                   for level in (0.1, 0.3, 0.5)) == pytest.approx(0.72,
                                                                  abs=0.01)


@pytest.mark.parametrize('panel, at_least', [
    ('a', 4.0),       # the vector printed in the script: 4.6 px
    ('b', 6.0)])      # 7.1 px
def test_printed_vector_k0_misses_the_published_contours(
        maps, digitized, panel, at_least):
    rms, _ = contour_scores(maps['k0_printed'], digitized[panel], panel)
    assert rms >= at_least, rms


@pytest.mark.parametrize('change, rms_a, rms_b', [
    (dict(time_offset=0.), 2.27, 3.83),             # plotted time 0 = feed
    (dict(initial_condition='fitting'), 12.8, 14.5)])
def test_properties_of_the_script_that_are_in_the_published_figure(
        digitized, change, rms_a, rms_b):
    # NOTES.md: without the 10 min offset of the time axis, or with the
    # script's fitting convention in place of its heat-map convention, the
    # same constants miss the published contours by far more than the
    # 1.11 and 0.35 px of the model as the script runs it. (Offsets this
    # large are linearized estimates, clipped at 30 px per sample.)
    km = load_o2_simplified_reversible_model(
        initial_condition=change.get('initial_condition', 'heatmap'))
    changed = simulate_heat_maps(km, time_offset=change.get('time_offset',
                                                            10.))
    assert contour_scores(changed, digitized['a'], 'a')[0] == pytest.approx(
        rms_a, rel=0.02)
    assert contour_scores(changed, digitized['b'], 'b')[0] == pytest.approx(
        rms_b, rel=0.02)


def test_late_yield_contours_lie_earlier_than_selectivity_contours(
        maps, digitized):
    # Yield <= selectivity everywhere, so in one run the late side of a
    # yield contour lies earlier than that of the selectivity contour of the
    # same level: by 0.005 to 0.09 % for the levels 0.1 and 0.3 above 62 C.
    T_C, times, selectivity, MA_yield = maps['deposit']
    rows = T_C > 62.
    row_T_C = T_C[rows]
    for level in (0.1, 0.3):
        late = []
        for matrix in (selectivity, MA_yield):
            shift = contour_shifts(T_C, times, matrix, level, 'late',
                                   np.full(rows.sum(), 1680.), row_T_C)
            late.append(1680. + shift)
        relative = late[1]/late[0] - 1.
        assert np.isfinite(relative).all() and len(relative) == 38
        assert -9e-4 < relative.min() and relative.max() < -5e-5
    # In the published figure the late 0.1 and 0.3 contours of panel b lie
    # LATER than those of panel a above 62 C: by 2.2 % in the median and up
    # to 4.4 % (panel a interpolated to the temperatures of the panel b
    # samples). No single run can produce both panels.
    for level, median, largest, median_px, largest_px in (
            (0.1, 0.0224, 0.0410, 2.35, 4.56),
            (0.3, 0.0221, 0.0436, 1.19, 3.63)):
        a, b = (table['contours'] for table in (digitized['a'],
                                                digitized['b']))
        a = a[(a['level'] == level) & (a['branch'] == 'late')].sort_values(
            'T_C')
        b = b[(b['level'] == level) & (b['branch'] == 'late')
              & (b['T_C'] > 62.) & b['T_C'].between(a['T_C'].min(),
                                                    a['T_C'].max())]
        on_a = np.interp(b['T_C'], a['T_C'], a['time'])
        relative = (b['time'] - on_a)/on_a
        assert (np.median(relative), relative.max()) == pytest.approx(
            (median, largest), abs=5e-4)
        pixels = (b['time'] - on_a)/pkg.FIG11_MIN_PER_PX
        assert (np.median(pixels), pixels.max()) == pytest.approx(
            (median_px, largest_px), abs=0.02)


def test_maxima_quoted_in_the_text(maps):
    T_C, times, selectivity, MA_yield = maps['deposit']
    row, column = list(T_C).index(95.), 3
    assert times[column] == pytest.approx(26.76, abs=0.005)
    # "59.1 % (95 C, 27 min)": 0.592 in that cell, the largest of the map
    assert selectivity[row, column] == pytest.approx(0.5924, abs=5e-4)
    assert selectivity[row, column] == selectivity.max()

    row, column = list(T_C).index(84.), 7
    assert times[column] == pytest.approx(49.10, abs=0.005)
    # "55.2 % (84 C, 49 min)": the cell is the largest of the map, and its
    # value is 0.548. The crest is flat (13 cells between 76 and 90 C are
    # within 0.001 of the maximum), so the cell is not required to be the
    # largest one.
    assert MA_yield[row, column] == pytest.approx(0.5478, abs=5e-4)
    assert MA_yield.max() - MA_yield[row, column] < 1e-3
    assert MA_yield.max() < 0.552
    crest = np.argwhere(MA_yield > MA_yield.max() - 1e-3)
    assert len(crest) == 13
    assert (T_C[crest[:, 0]].min(), T_C[crest[:, 0]].max()) == (76., 90.)


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
    # 450 min. The 50 C row of the computed map peaks at 0.5704.
    assert X_MA.max() < 0.569 - 0.02
    late = time >= 240.
    per_converted = X_MA[late]/(1. - X_HP3[late])
    assert per_converted[:3] == pytest.approx((0.6207, 0.5774, 0.5128),
                                              abs=5e-5)
    assert np.all(np.diff(per_converted) < 0.)
    assert np.abs(per_converted - 0.569).min() > 0.008
    T_C, _, selectivity, _ = maps['deposit']
    assert selectivity[list(T_C).index(50.)].max() == pytest.approx(
        0.5704, abs=2e-4)


def test_50C_row_maxima(maps):
    T_C, times, selectivity, MA_yield = maps['deposit']
    row = list(T_C).index(50.)
    assert selectivity[row].max() == pytest.approx(0.5704, abs=2e-4)
    assert heat_map_ridge(times, selectivity)[row] == pytest.approx(
        261.3, abs=0.1)
    assert MA_yield[row].max() == pytest.approx(0.5325, abs=2e-4)
    assert heat_map_ridge(times, MA_yield)[row] == pytest.approx(322.8,
                                                                 abs=0.1)


def test_heat_map_maximum_and_ridge(maps):
    T_C, times, selectivity, _ = maps['deposit']
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


@pytest.mark.parametrize('panel, mean, rms, largest', [
    ('a', 0.05, 0.83, 1.83), ('b', 0.05, 1.52, 6.48)])
def test_ridge_is_the_published_red_dashed_line(maps, digitized, panel, mean,
                                                rms, largest):
    # The red dashed line joins the row maxima, so it is defined to one
    # grid column: 5.585 min, 1.82 px of the published raster. It is
    # therefore held to the digitized line to half a column in the mean
    # (a line that sits one column off is 1.82 px off in the mean) and to
    # one column in the rms. The largest single offsets (1.8 and 6.5 px)
    # are at 32 to 34 C, where the row maxima are flat and the line
    # advances several columns per row.
    T_C, times = maps['deposit'][:2]
    matrix = maps['deposit'][2 if panel == 'a' else 3]
    column = (times[1] - times[0])/pkg.FIG11_MIN_PER_PX
    assert column == pytest.approx(1.818, abs=1e-3)
    assert len(digitized[panel]['ridge']) == {'a': 118, 'b': 226}[panel]
    offsets = ridge_offsets(T_C, heat_map_ridge(times, matrix),
                            digitized[panel])
    assert offsets.mean() == pytest.approx(mean, abs=0.02)
    assert abs(offsets.mean()) <= 0.5*column
    assert np.sqrt(np.mean(offsets**2)) == pytest.approx(rms, abs=0.02)
    assert np.sqrt(np.mean(offsets**2)) <= column
    assert np.abs(offsets).max() == pytest.approx(largest, abs=0.02)
    assert digitized[panel]['ridge']['T_C'].to_numpy()[
        np.abs(offsets).argmax()] < 34.
    # negative control: the same line one column later or earlier fails
    # the limit on the mean (1.86 and -1.77 px in panel a, 1.87 and
    # -1.76 px in panel b)
    for shift in (1, -1):
        shifted = times[np.argmax(matrix, axis=1) + shift]
        assert abs(ridge_offsets(T_C, shifted, digitized[panel]).mean()
                   ) > 0.5*column + 0.2


@pytest.mark.parametrize('panel, largest', [('a', 1.12e-4), ('b', 1.36e-4)])
def test_solver_of_the_script_does_not_matter_for_the_maps(maps, panel,
                                                           largest):
    # The script integrates each of the 299 intervals of the time grid with
    # ode23s at default tolerances; the maps here are converged solutions.
    # The 100 C row, repeated as the script computes it, differs by
    # 1.1e-4 (a) and 1.4e-4 (b).
    T_C, times, selectivity, MA_yield = maps['deposit']
    k = script_k_adj(DEPOSIT, 373.15)
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


# --- 6. the fit that defines k_3 and k_4 ---------------------------------------

def test_data_table_of_the_script():
    # The 16 samples of the script's function data(), with its digits.
    data = pkg.read_heatmap_script_data()
    assert list(data.columns) == ['time', 'X_HP3', 'X_MA', 'X_AAc', 'X_OXA',
                                  'X_D']
    assert list(data['time']) == [0, 20, 40, 60, 90, 120, 180, 240, 360, 450,
                                  540, 1110, 1280, 1360, 1470, 1680]
    assert list(data.iloc[0]) == [0., 1., 0., 0., 0., 0.]
    assert list(data.iloc[1]) == [20., 0.91918, 0.06510, 0.00566, 0.00229,
                                  0.00777]
    assert list(data.iloc[-2]) == [1470., 0.0245, 0.1242, 0.0158, 0.0323,
                                   0.8033]
    assert list(data.iloc[-1]) == [1680., 0.0275, 0.0953, 0.0381, 0.0200,
                                   0.8191]
    # It is the experimental block of the deposited Fig. 9b table, except
    # that the script gives the last two rows to four decimals: they are the
    # rounded values of the file (largest difference 5e-5).
    deposited = pkg.read_fig9b_o2()['experiment']
    assert data.iloc[:14].equals(deposited.iloc[:14])
    last = deposited.iloc[14:].to_numpy()
    assert not np.array_equal(data.iloc[14:].to_numpy(), last)
    assert np.allclose(data.iloc[14:].to_numpy()[:, 1:],
                       np.round(last[:, 1:], 4), rtol=0., atol=1e-12)
    assert np.abs(data.to_numpy() - deposited.to_numpy()).max() == (
        pytest.approx(5e-5, abs=1e-9))
    # the CSV holds the digits of the script
    with open(os.path.join(pkg.reference_data_dir,
                           'heatmap_script_50C_data.csv')) as f:
        lines = f.read().splitlines()
    assert lines[0] == 't_min,3HP,MA,AA,OXA,D' and len(lines) == 17
    assert lines[1] == '0,1.00000,0.00000,0.00000,0.00000,0.00000'
    assert lines[-1] == '1680,0.0275,0.0953,0.0381,0.0200,0.8191'


@pytest.fixture(scope='module')
def default_fit():
    return fit_simplified_reversible_constants()


def test_default_k3_and_k4_are_the_result_of_the_fit(default_fit):
    # The shipped k_3 and k_4 ARE the fit, to the five significant digits
    # shipped: with every other constant at its deposited value, least
    # squares on the 80 residuals of the script's data table.
    deposit = pkg.SIMPLIFIED_REVERSIBLE_PARAMETER_SETS['deposit']
    assert list(default_fit) == ['fitted', 'standard_errors',
                                 'sum_of_squares', 'n_residuals', 'values']
    assert list(default_fit['fitted']) == ['k_3_ref', 'k_4_ref']
    for name in ('k_3_ref', 'k_4_ref'):
        assert float('%.5g' % default_fit['fitted'][name]) == deposit[name]
    assert (deposit['k_3_ref'], deposit['k_4_ref']) == (1.4200e-5, 2.3204e-3)
    # neither value is near a rounding boundary of its fifth digit
    assert default_fit['fitted']['k_3_ref'] == pytest.approx(1.42000e-5,
                                                             rel=1e-5)
    assert default_fit['fitted']['k_4_ref'] == pytest.approx(2.32041e-3,
                                                             rel=5e-6)
    assert default_fit['n_residuals'] == 80
    assert default_fit['sum_of_squares'] == pytest.approx(0.0337728,
                                                          abs=2e-7)
    # the other eight constants are returned as they were held
    assert {name: value for name, value in default_fit['values'].items()
            if name not in default_fit['fitted']} == {
                name: value for name, value in deposit.items()
                if name not in default_fit['fitted']}
    # standard errors: 31 % of k_3 and 3.7 % of k_4
    errors = default_fit['standard_errors']
    assert errors['k_3_ref'] == pytest.approx(4.42e-6, rel=0.01)
    assert errors['k_4_ref'] == pytest.approx(8.62e-5, rel=0.01)
    assert errors['k_3_ref']/deposit['k_3_ref'] == pytest.approx(0.311,
                                                                 abs=0.005)
    assert errors['k_4_ref']/deposit['k_4_ref'] == pytest.approx(0.037,
                                                                 abs=0.001)


@pytest.mark.parametrize('start', [(4.26e-5, 1.15e-3), (4.3e-6, 3.45e-3),
                                   (1e-4, 1e-4)])
def test_fit_is_reached_from_other_starting_points(default_fit, start):
    fit = fit_simplified_reversible_constants(start=start)
    for name in ('k_3_ref', 'k_4_ref'):
        # to 1e-5: k_3 has a standard error of 31 %, so the sum of squares
        # is flat about the optimum (the fifth digit shipped needs 3.5e-5)
        assert fit['fitted'][name] == pytest.approx(
            default_fit['fitted'][name], rel=1e-5)
    assert fit['sum_of_squares'] == pytest.approx(
        default_fit['sum_of_squares'], rel=1e-9)


def test_fit_is_deterministic_and_is_a_minimum(default_fit):
    again = fit_simplified_reversible_constants()
    assert again == default_fit
    # the sum of squares of the shipped set, and of sets with k_3 or k_4
    # moved by 1 %: the fit is the minimum
    data = pkg.read_heatmap_script_data()

    def sum_of_squares(**changes):
        km = load_o2_simplified_reversible_model(initial_condition='fitting')
        for name, factor in changes.items():
            km.set_value(name, factor*km.get_value(name))
        modeled = simulate_at_times(km, data['time'])
        return ((modeled[list(pkg.NORMALIZED)].to_numpy()
                 - data[list(pkg.NORMALIZED)].to_numpy())**2).sum()

    shipped = sum_of_squares()
    assert shipped == pytest.approx(default_fit['sum_of_squares'], abs=1e-8)
    for name in ('k_3_ref', 'k_4_ref'):
        for factor in (0.99, 1.01):
            assert sum_of_squares(**{name: factor}) > shipped + 1e-8


def test_fit_of_all_five_rate_constants_returns_the_arrhenius_row():
    # NOTES.md: with K_6 to K_10 as printed in the script and all five rate
    # constants fitted to the same data, k_1, k_2 and k_5 come out 0.09,
    # 1.5 and 0.37 % above the 50 C row of the deposited Arrhenius table,
    # well inside their standard errors (2.5, 31 and 3.9 %). This is the
    # support for taking that row as this model's constants.
    deposit = pkg.SIMPLIFIED_REVERSIBLE_PARAMETER_SETS['deposit']
    fit = fit_simplified_reversible_constants(
        free=('k_1', 'k_2', 'k_3', 'k_4', 'k_5'))
    fitted = fit['fitted']
    assert [fitted['k_%d_ref' % i] for i in range(1, 6)] == pytest.approx(
        [2.14192e-3, 5.27746e-5, 1.41648e-5, 2.31221e-3, 1.06641e-3],
        rel=2e-5)
    assert [fitted[name]/deposit[name] - 1.
            for name in ('k_1_ref', 'k_2_ref', 'k_5_ref')] == pytest.approx(
                [0.0009, 0.0149, 0.0037], abs=1e-4)
    assert [fit['standard_errors'][name]/fitted[name]
            for name in ('k_1_ref', 'k_2_ref', 'k_5_ref')] == pytest.approx(
                [0.025, 0.311, 0.039], abs=0.002)
    assert fit['sum_of_squares'] == pytest.approx(0.033768, abs=2e-6)
    # k_3 and k_4 move by 0.25 and 0.35 %
    assert fitted['k_3_ref']/deposit['k_3_ref'] == pytest.approx(0.9975,
                                                                 abs=2e-4)
    assert fitted['k_4_ref']/deposit['k_4_ref'] == pytest.approx(0.9965,
                                                                 abs=2e-4)


def test_fit_arguments():
    deposit = pkg.SIMPLIFIED_REVERSIBLE_PARAMETER_SETS['deposit']
    # model-parameter names are accepted, and held values can be given
    fit = fit_simplified_reversible_constants(
        free=('k_4_ref',), values={'k_3_ref': 2e-5})
    assert list(fit['fitted']) == ['k_4_ref']
    assert fit['values']['k_3_ref'] == 2e-5
    assert fit['values']['k_1_ref'] == deposit['k_1_ref']
    assert fit['sum_of_squares'] > 0.0337728
    # other data: the experimental block of the deposited Fig. 9b table,
    # which differs from the script's table in the rounding of two rows
    other = fit_simplified_reversible_constants(
        data=pkg.read_fig9b_o2()['experiment'])
    assert other['fitted']['k_3_ref'] == pytest.approx(1.41979e-5, rel=1e-5)
    assert other['fitted']['k_4_ref'] == pytest.approx(2.32047e-3, rel=1e-5)
    for arguments in (dict(free=('k_6',)), dict(free=('k_3', 'k_3_ref')),
                      dict(free=()), dict(values={'K_11_ref': 1.}),
                      dict(start=(1e-5,)), dict(start=(1e-5, 0.)),
                      dict(start=(1e-5, -1e-3))):
        with pytest.raises(ValueError):
            fit_simplified_reversible_constants(**arguments)


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
