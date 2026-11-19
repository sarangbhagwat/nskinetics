# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests of the reversible O2 variant and of the Fig. 10 data (Kim et al.,
Appl. Catal. B: Environ. Energy 386 (2026) 126403) in
``nskinetics.models.pd_c_3hp_oxidation``.

The reversible O2 variant is not in the paper or the authors' scripts: its
structure and constants were inferred from the model curves deposited for
the experiments at 25, 75 and 100 C (Fig. 10a-c). It is checked in three
ways:

* its rate laws against an independent evaluation of the documented
  equations (``NOTES.md``), written out again in this module on their own
  basis (mol per mol of initial 3-HP);
* against the deposited model curves, at their node times. Those curves are
  MATLAB ``ode23s`` output at default tolerances, so an accurate solution
  is expected to differ from them by up to 1.8e-4 (25 C), 2.5e-3 (75 C) and
  2.7e-3 (100 C). The limits are 1.5 to 2 times that solver error.
  Negative controls at 25 C, where the solver error is smallest, show what
  they exclude: water formed by reaction 1 only, no reverse of reaction 1,
  and a 10 % error in any of eight of the ten constants;
* against an accurate solution of the same equations with the same
  constants, computed in this module with scipy (``solve_ivp``, DOP853,
  relative tolerance 1e-13), independently of the Antimony file and of
  roadrunner, to 1e-7. The model at the loader's tolerances is within
  2.5e-8 of it, and a 1 % error in any single constant moves it by 1.5e-7
  or more. This guards the implementation, not the inference.

The straight lines of the Arrhenius plot (Fig. 10d) are recomputed from the
deposited points and compared with the deposited lines. The last section
pins the numbers ``NOTES.md`` reports for the per-temperature fits and for
the Table 3 O2 model against the experiments at 25, 75 and 100 C.
``test_pd_c_3hp_oxidation.py`` checks that the selected constant set
survives ``compile_events``.

The models are small, so this module is not marked slow. It imports pytest
at top level and is therefore not registered in
``nskinetics/tests/__init__.py``.
"""
import hashlib
import itertools
import os

import numpy as np
import pytest
import tellurium as te
from scipy.integrate import solve_ivp

import nskinetics as nsk
import nskinetics.models.pd_c_3hp_oxidation as pkg
from nskinetics.models.pd_c_3hp_oxidation import (
    NORMALIZED, SPECIES, load_o2_model, load_o2_reversible_model,
    pooled_r2, rms_difference, simulate_at_times)

_PKG_DIR = os.path.dirname(os.path.abspath(pkg.__file__))
_NORMALIZED = list(NORMALIZED)
_STATES = SPECIES + ('W',)
_CONCENTRATIONS = ['[%s]' % s for s in _STATES]
_NAMES = ('k_1', 'k_1r', 'k_2', 'k_2r', 'k_3', 'k_4', 'k_5', 'k_5r', 'k_6',
          'k_6r')
C0 = 0.06                       # M
TEMPERATURES_C = (25, 75, 100)

# The documented constant sets (NOTES.md), in the order of _NAMES, 1/min.
CONSTANTS = {
    25: (0.0465036, 51.0644, 0.0088933, 0.00136441, 2.43424e-4, 6.17632e-5,
         1.01563e-6, 1.47217e-7, 0.00284329, 3.12354e-4),
    75: (0.00918893, 1.88278e-5, 0.0490036, 0.00193586, 4.9301e-4,
         7.97777e-6, 0.00160383, 2.86632e-5, 0.00613222, 2.56106e-5),
    100: (0.0414885, 0.00474984, 0.122878, 0.0103261, 0.00169028,
          7.86978e-6, 9.32644e-7, 0., 0.0157177, 1.05066e-4),
    }

# Limits for an accurate solution against the deposited model curve:
# rms and maximum over all values, and maximum per species
# (3-HP, MA, AAc, OXA, D).
LIMITS = {
    25: (8e-5, 3.5e-4, (2e-4, 3.5e-4, 2e-5, 1.5e-5, 2e-4)),
    75: (1.2e-3, 4e-3, (3e-3, 4e-3, 3e-4, 1.5e-5, 3e-3)),
    100: (1.2e-3, 4.5e-3, (3e-3, 4.5e-3, 4e-4, 1e-5, 3e-3)),
    }

#              R1  R2  R3  R4  R5  R6
NU = np.array([[-1, 0, 0, 0, -1, 0],     # 3-HP
               [1, -1, 0, 0, 0, 0],      # OPA
               [0, 1, -1, -1, 0, -1],    # MA
               [0, 0, 1, 0, 0, 0],       # AAc
               [0, 0, 0, 1, 0, 0],       # OXA
               [0, 0, 0, 0, 1, 1],       # D
               [1, 0, 0, 1, 0, 0]])      # W, water formed by R1 and R4


def documented_rates(p, c):
    """Rates of the reversible variant as documented in ``NOTES.md``, in
    1/min, for ``c`` in mol per mol of initial 3-HP: (3-HP, OPA, MA, AAc,
    OXA, D, W)."""
    hp3, opa, ma, _, _, d, w = c
    k_1, k_1r, k_2, k_2r, k_3, k_4, k_5, k_5r, k_6, k_6r = p
    return np.array([k_1*hp3 - k_1r*opa*w,
                     k_2*opa - k_2r*ma,
                     k_3*ma,
                     k_4*ma,
                     k_5*hp3 - k_5r*d,
                     k_6*ma - k_6r*d])


def deviations(sim, nodes):
    """(rms over all, max over all, max per species) of sim - nodes."""
    diff = sim[_NORMALIZED].to_numpy() - nodes[_NORMALIZED].to_numpy()
    return (np.sqrt((diff**2).mean()), np.abs(diff).max(),
            np.abs(diff).max(axis=0))


def within_limits(sim, nodes, temperature_C):
    rms, max_abs, max_per_species = deviations(sim, nodes)
    rms_limit, max_limit, species_limits = LIMITS[temperature_C]
    return bool(rms <= rms_limit and max_abs <= max_limit
                and (max_per_species <= np.array(species_limits)).all())


def label(temperature_C):
    """Name of the constant set of the experiment at ``temperature_C``."""
    return '%dC' % temperature_C


@pytest.fixture(scope='module', params=TEMPERATURES_C)
def case(request):
    """(temperature in C, model, deposited tables) of one Fig. 10 panel."""
    temperature_C = request.param
    return (temperature_C, load_o2_reversible_model(label(temperature_C)),
            pkg.read_fig10_o2(temperature_C))


# --- 1. loader ----------------------------------------------------------------

def test_loader_returns_a_configured_kinetic_model():
    km = load_o2_reversible_model()
    assert isinstance(km, nsk.KineticModel)
    assert km._units == {'time': 'min', 'conc': 'M'}
    km.validate_units()
    assert tuple(km._te.getFloatingSpeciesIds()) == _STATES
    assert km._te.integrator.absolute_tolerance == 1e-10
    assert km._te.integrator.relative_tolerance == 1e-9
    assert load_o2_reversible_model() is not km
    # the feed of Fig. 10: 60 mM 3-HP, no products, no water of reaction
    assert km.get_value('[HP3]') == pytest.approx(C0)
    for product in _STATES[1:]:
        assert km.get_value('[%s]' % product) == 0.
    assert km.get_value('X_HP3') == pytest.approx(1.)
    assert km.get_value('x_OPA') == 0.


def test_constant_sets_are_the_documented_ones():
    assert list(pkg.O2_REVERSIBLE_PARAMETER_SETS) == ['25C', '75C', '100C']
    for temperature_C, values in CONSTANTS.items():
        expected = dict(zip(_NAMES, values))
        assert pkg.O2_REVERSIBLE_PARAMETER_SETS[label(temperature_C)] == (
            expected)
        km = load_o2_reversible_model(label(temperature_C))
        for name, value in expected.items():
            assert km.get_value(name) == value, (temperature_C, name)


def test_antimony_file_holds_the_75C_set():
    path = os.path.join(_PKG_DIR, pkg.o2_reversible_antimony_filename)
    r = te.loadAntimonyModel(path)
    default = load_o2_reversible_model()
    for name, value in zip(_NAMES, CONSTANTS[75]):
        assert r[name] == value, name
        assert default.get_value(name) == value, name


@pytest.mark.parametrize('parameter_set', ['50C', '75', '75 C', 75, 348.15,
                                           25, None, ['75C']])
def test_loader_has_no_other_parameter_set(parameter_set):
    # The sets are named, not addressed by a temperature: a number is
    # refused whatever its unit.
    with pytest.raises(ValueError, match='unknown parameter set') as error:
        load_o2_reversible_model(parameter_set)
    message = str(error.value)
    assert "'25C', '75C', '100C'" in message
    # points to the model of the 50 C curve (Fig. 9b)
    assert 'load_o2_model' in message


# --- 2. rate laws and invariants ------------------------------------------------

@pytest.mark.parametrize('temperature_C', TEMPERATURES_C)
def test_rate_laws_are_the_documented_ones(temperature_C):
    km = load_o2_reversible_model(label(temperature_C))
    # a state with every species present, in mol per mol of initial 3-HP
    c = np.array([0.55, 0.07, 0.20, 0.03, 0.01, 0.14, 0.31])
    for species, value in zip(_STATES, c):
        km.set_value('[%s]' % species, C0*value)
    rates = documented_rates(CONSTANTS[temperature_C], c)
    # every term is active in this state (except the zero k_5r at 100 C)
    p = np.array(CONSTANTS[temperature_C])
    assert abs(p[1]*c[1]*c[6]) > 1e-9 and abs(p[3]*c[2]) > 1e-9
    for j, rate in enumerate(rates, start=1):
        assert km.get_value('R%d' % j) == pytest.approx(C0*rate, rel=1e-12), j
    for species, derivative in zip(_STATES, NU @ rates):
        assert km.get_value("[%s]'" % species) == pytest.approx(
            C0*derivative, rel=1e-11, abs=1e-18), species
    total = sum(km.get_value('[%s]' % s) for s in ('HP3', 'MA', 'AAc', 'OXA',
                                                   'D'))
    for name, species in zip(_NORMALIZED, ('HP3', 'MA', 'AAc', 'OXA', 'D')):
        assert km.get_value(name) == pytest.approx(
            km.get_value('[%s]' % species)/total, rel=1e-12)
    assert km.get_value('x_OPA') == pytest.approx(0.000339*c[1], rel=1e-12)


def test_reverse_of_reaction_1_needs_the_water_of_reaction():
    km = load_o2_reversible_model('25C')
    km.set_value('[HP3]', 0.05)
    km.set_value('[OPA]', 0.002)
    forward = CONSTANTS[25][0]*0.05
    assert km.get_value('[W]') == 0.
    assert km.get_value('R1') == pytest.approx(forward, rel=1e-12)
    km.set_value('[W]', 0.002)
    assert km.get_value('R1') == pytest.approx(
        forward - CONSTANTS[25][1]*0.002*0.002/C0, rel=1e-12)


def test_carbon_species_are_conserved_and_water_counts_R1_and_R4(case):
    temperature_C, km, tables = case
    t_end = tables['model']['time'].iloc[-1]
    times = np.arange(0., t_end + 0.125, 0.25)
    sim = simulate_at_times(km, times, _CONCENTRATIONS + ['R1', 'R4'])
    carbon = sim[_CONCENTRATIONS[:-1]].sum(axis=1)
    assert np.allclose(carbon, C0, rtol=1e-6, atol=0.)
    assert (sim[_CONCENTRATIONS].to_numpy() > -1e-9).all()
    # W is the time integral of R1 + R4 (trapezoidal rule on a fine grid)
    rate = (sim['R1'] + sim['R4']).to_numpy()
    integral = (0.5*(rate[1:] + rate[:-1])*np.diff(times)).sum()
    assert sim['[W]'].iloc[-1] == pytest.approx(integral, rel=1e-4)
    assert sim['[W]'].iloc[-1] > 0.1*C0


# --- 3. against the deposited model curves (Fig. 10a-c) -------------------------

def test_model_matches_the_deposited_curve(case):
    temperature_C, km, tables = case
    nodes = tables['model']
    sim = simulate_at_times(km, nodes['time'])
    rms, max_abs, max_per_species = deviations(sim, nodes)
    rms_limit, max_limit, species_limits = LIMITS[temperature_C]
    assert rms <= rms_limit, rms
    assert max_abs <= max_limit, max_abs
    assert (max_per_species <= np.array(species_limits)).all(), max_per_species
    assert within_limits(sim, nodes, temperature_C)
    # the difference is the solver error of the deposited curve (NOTES.md)
    expected = {25: (3.8e-5, 1.8e-4), 75: (7.1e-4, 2.5e-3),
                100: (6.9e-4, 2.7e-3)}[temperature_C]
    assert (rms, max_abs) == pytest.approx(expected, rel=0.03)


def test_straight_segments_between_the_nodes_of_the_deposited_curve(case):
    # The published lines join the solver's steps; between the steps they
    # are further from the solution than at the steps (NOTES.md).
    temperature_C, km, tables = case
    nodes = tables['model']
    times = np.arange(0., nodes['time'].iloc[-1] + 0.05, 0.1)
    sim = simulate_at_times(km, times)
    largest = max(np.abs(sim[name] - np.interp(times, nodes['time'],
                                               nodes[name])).max()
                  for name in _NORMALIZED)
    expected = {25: 6.6e-4, 75: 9.7e-3, 100: 8.7e-3}[temperature_C]
    assert largest == pytest.approx(expected, rel=0.05)


def _variant_25C(antimony_edit=None, **changes):
    """The 25 C model with an edited Antimony text or changed constants."""
    path = os.path.join(_PKG_DIR, pkg.o2_reversible_antimony_filename)
    with open(path, encoding='utf-8') as f:
        text = f.read()
    if antimony_edit is not None:
        old, new = antimony_edit
        assert text.count(old) == 1
        text = text.replace(old, new)
    km = nsk.KineticModel(te.loada(text), units={'time': 'min', 'conc': 'M'})
    km._te.integrator.absolute_tolerance = 1e-10
    km._te.integrator.relative_tolerance = 1e-9
    for name, value in dict(zip(_NAMES, CONSTANTS[25]), **changes).items():
        km.set_value(name, value)
    return km


def test_unedited_variant_is_the_25C_model():
    nodes = pkg.read_fig10_o2(25)['model']
    assert simulate_at_times(_variant_25C(), nodes['time']).equals(
        simulate_at_times(load_o2_reversible_model('25C'), nodes['time']))


def test_water_from_reaction_1_only_misses_the_25C_curve():
    nodes = pkg.read_fig10_o2(25)['model']
    km = _variant_25C(('R4: MA => OXA + W;', 'R4: MA => OXA;'))
    sim = simulate_at_times(km, nodes['time'])
    assert not within_limits(sim, nodes, 25)
    rms, max_abs, _ = deviations(sim, nodes)
    assert rms > 2e-4 and max_abs > 8e-4


def test_no_reverse_of_reaction_1_misses_the_25C_curve():
    nodes = pkg.read_fig10_o2(25)['model']
    sim = simulate_at_times(_variant_25C(k_1r=0.), nodes['time'])
    assert not within_limits(sim, nodes, 25)
    assert deviations(sim, nodes)[1] > 0.9


@pytest.mark.parametrize('factor', [0.9, 1.1])
@pytest.mark.parametrize('name', ['k_1', 'k_1r', 'k_2', 'k_2r', 'k_3', 'k_4',
                                  'k_6', 'k_6r'])
def test_a_10_percent_error_in_a_constant_misses_the_25C_curve(name, factor):
    # k_5 and k_5r are left out: they are too small at 25 C to be seen.
    nodes = pkg.read_fig10_o2(25)['model']
    value = factor*dict(zip(_NAMES, CONSTANTS[25]))[name]
    sim = simulate_at_times(_variant_25C(**{name: value}), nodes['time'])
    assert not within_limits(sim, nodes, 25)


# --- 4. against an accurate, independent solution of the same equations ----------

# Limit for the model against the accurate solution. CVODE at the loader's
# tolerances is within 2.5e-8 of it.
EXACT_LIMIT = 1e-7


@pytest.fixture(scope='module')
def exact():
    """Accurate solutions of the documented equations, by temperature in C.

    The equations of ``documented_rates`` with the constants of
    ``CONSTANTS``, on their own basis (mol per mol of initial 3-HP),
    integrated with scipy (DOP853, relative tolerance 1e-13): columns
    ``time``, the five plotted fractions, and OPA and W per initial 3-HP.
    Times: 601 equally spaced values from 0 to the end of the deposited
    curve, with its node times and the sampling times of the experiment.
    """
    solutions = {}
    for temperature_C in TEMPERATURES_C:
        tables = pkg.read_fig10_o2(temperature_C)
        t_end = tables['model']['time'].iloc[-1]
        times = np.unique(np.concatenate((
            np.linspace(0., t_end, 601), tables['model']['time'],
            tables['experiment']['time'])))
        p = CONSTANTS[temperature_C]
        solution = solve_ivp(lambda t, c: NU @ documented_rates(p, c),
                             (0., t_end), [1., 0., 0., 0., 0., 0., 0.],
                             method='DOP853', t_eval=times, rtol=1e-13,
                             atol=1e-15)
        assert solution.success
        c = solution.y
        measured = c[[0, 2, 3, 4, 5]]
        solutions[temperature_C] = (
            times, (measured/measured.sum(axis=0)).T, c[1], c[6])
    return solutions


def test_accurate_solution_is_converged(exact):
    # The reference does not depend on its solver: an implicit multistep
    # method (LSODA) at a relative tolerance of 1e-12 gives the same
    # fractions to 1e-10 (they agree to 1e-12).
    for temperature_C in TEMPERATURES_C:
        times, fractions, _, _ = exact[temperature_C]
        p = CONSTANTS[temperature_C]
        solution = solve_ivp(lambda t, c: NU @ documented_rates(p, c),
                             (0., times[-1]), [1., 0., 0., 0., 0., 0., 0.],
                             method='LSODA', t_eval=times, rtol=1e-12,
                             atol=1e-15)
        measured = solution.y[[0, 2, 3, 4, 5]]
        assert np.abs((measured/measured.sum(axis=0)).T
                      - fractions).max() < 1e-10


def test_model_matches_the_accurate_solution(case, exact):
    temperature_C, km, _ = case
    times, fractions, opa, water = exact[temperature_C]
    assert len(times) > 600 and (np.diff(times) > 0.).all()
    sim = simulate_at_times(km, times, _NORMALIZED + ['[OPA]', '[W]'])
    assert np.abs(sim[_NORMALIZED].to_numpy() - fractions).max() <= EXACT_LIMIT
    assert np.abs(sim['[OPA]']/C0 - opa).max() <= EXACT_LIMIT
    assert np.abs(sim['[W]']/C0 - water).max() <= EXACT_LIMIT
    # NOTES.md: the model is within 2.5e-8 of the accurate solution
    assert np.abs(sim[_NORMALIZED].to_numpy() - fractions).max() <= 4e-8


@pytest.mark.parametrize('error, smallest_effect', [(0.1, 1.57e-6),
                                                    (0.01, 1.57e-7)])
def test_accurate_solution_detects_an_error_in_any_constant(
        exact, error, smallest_effect):
    # The constant with the smallest effect is k_5r at 25 C (NOTES.md); the
    # effect of a 1 % error in it is still above the limit.
    smallest = np.inf
    for temperature_C in TEMPERATURES_C:
        times, fractions, _, _ = exact[temperature_C]
        km = load_o2_reversible_model(label(temperature_C))
        for name, value in zip(_NAMES, CONSTANTS[temperature_C]):
            if value == 0.:
                continue
            for factor in (1. - error, 1. + error):
                km.set_value(name, factor*value)
                sim = simulate_at_times(km, times)
                smallest = min(smallest, np.abs(
                    sim[_NORMALIZED].to_numpy() - fractions).max())
            km.set_value(name, value)
    assert smallest > EXACT_LIMIT
    assert smallest == pytest.approx(smallest_effect, rel=0.2)


def test_unmeasured_opa_pool_and_water_of_reaction():
    # peak of OPA and final W, in mol per mol of initial 3-HP (NOTES.md)
    for temperature_C, opa_peak, t_peak, w_end in (
            (25, 0.0294, 2.3, 0.133), (75, 0.1255, 42.2, 0.881),
            (100, 0.2010, 14.6, 1.000)):
        t_end = pkg.read_fig10_o2(temperature_C)['model']['time'].iloc[-1]
        times = np.concatenate((np.arange(0., 60., 0.1), [t_end]))
        sim = simulate_at_times(
            load_o2_reversible_model(label(temperature_C)), times,
            ['[OPA]', '[W]'])
        peak = sim['[OPA]'].idxmax()
        assert sim['[OPA]'][peak]/C0 == pytest.approx(opa_peak, abs=2e-4)
        assert sim['time'][peak] == pytest.approx(t_peak, abs=0.15)
        assert sim['[W]'].iloc[-1]/C0 == pytest.approx(w_end, abs=1e-3)


def test_ratio_of_the_constants_of_reaction_1_is_not_monotonic():
    ratio = [CONSTANTS[T][0]/CONSTANTS[T][1] for T in TEMPERATURES_C]
    assert ratio == pytest.approx([9.107e-4, 488.1, 8.735], rel=1e-3)


def test_forward_constants_in_the_units_of_table_3():
    # k/(m_cat*x_O2_ref) for R1, R2 and R4, k/m_cat for R3, R5 and R6
    # (NOTES.md; for comparison with Table 3 only)
    m_cat, x_O2_ref = 0.2, 15*0.000339
    expected = {25: (45.7, 8.74, 1.22e-3, 6.07e-2, 5.08e-6, 1.42e-2),
                75: (9.04, 48.2, 2.47e-3, 7.84e-3, 8.02e-3, 3.07e-2),
                100: (40.8, 121., 8.45e-3, 7.74e-3, 4.66e-6, 7.86e-2)}
    for temperature_C, values in expected.items():
        p = dict(zip(_NAMES, CONSTANTS[temperature_C]))
        converted = [p['k_1']/(m_cat*x_O2_ref), p['k_2']/(m_cat*x_O2_ref),
                     p['k_3']/m_cat, p['k_4']/(m_cat*x_O2_ref),
                     p['k_5']/m_cat, p['k_6']/m_cat]
        assert converted == pytest.approx(values, rel=3e-3)


# --- 5. reference data of Fig. 10a-c ---------------------------------------------

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
    # the argument is a temperature in degrees Celsius; 348.15 is 75 C in K
    for temperature_C in (50, 20, '25', '75C', 348.15):
        with pytest.raises(ValueError, match='no Fig. 10 data'):
            pkg.read_fig10_o2(temperature_C)


# --- 6. Arrhenius plot (Fig. 10d) ------------------------------------------------

# The plotted rate constants in 1/min, taken to be constants of the
# simplified reversible model (NOTES.md): the shortest decimals whose
# logarithms are the deposited points.
FIG10D_K = {'k_1': (7.75e-5, 2.14e-3, 2.43685e-3, 1.00763e-2),
            'k_2': (6.64e-4, 5.2e-5, 4.44096e-4, 1.59645e-3),
            'k_5': (1.676e-3, 1.0625e-3, 5.4797e-3, 1.41038e-2)}
# Slope (K), value at 1/373.15 K (the first row of the deposited line) and
# activation energy (kJ/mol) of the three lines.
FIG10D_LINES = {'k_1': (-7208.69, -4.61102, 59.93),
                'k_2': (-8293.97, -6.32951, 68.96),
                'k_5': (-6265.49, -4.17048, 52.09)}


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
        assert np.abs(np.log(k) - points['ln_' + name]).max() <= 1e-5


def test_fig10d_lines_are_least_squares_through_three_points():
    assert pkg.FIG10D_REGRESSION_TEMPERATURES_C == {
        'k_1': (25, 75, 100), 'k_2': (50, 75, 100), 'k_5': (50, 75, 100)}
    lines = pkg.fig10d_regression_lines()
    deposited = pkg.read_fig10d_arrhenius()
    for name, (slope, at_100C, Ea) in FIG10D_LINES.items():
        line = lines[name]
        assert line['slope'] == pytest.approx(slope, rel=1e-4)
        assert line['Ea'] == -8.314*line['slope']
        assert line['Ea']/1000. == pytest.approx(Ea, abs=0.005)
        assert line['temperatures_C'] == (
            pkg.FIG10D_REGRESSION_TEMPERATURES_C[name])
        drawn = deposited['line_' + name]['ln_' + name]
        assert drawn.iloc[0] == at_100C
        assert line['slope']/373.15 + line['intercept'] == pytest.approx(
            at_100C, abs=2e-5)
        # the published line ends where the deposited line ends
        T_end = pkg.FIG10D_LINE_ENDS_C[name][1]
        assert (line['slope']/(T_end + 273.15) + line['intercept']
                == pytest.approx(drawn.iloc[-1], abs=2e-4))
    k_1 = lines['k_1']
    assert k_1['slope']/298.15 + k_1['intercept'] == pytest.approx(
        -9.47060, abs=2e-5)


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


# --- 7. models against the experiments at 25, 75 and 100 C (NOTES.md) -----------

def test_per_temperature_sets_against_the_experiments():
    for temperature_C, r2, rms in ((25, 0.9996, 0.0079), (75, 0.9899, 0.0293),
                                   (100, 0.9926, 0.0257)):
        experiment = pkg.read_fig10_o2(temperature_C)['experiment']
        sim = simulate_at_times(
            load_o2_reversible_model(label(temperature_C)),
            experiment['time'])
        assert pooled_r2(experiment, sim) == pytest.approx(r2, abs=2e-4)
        assert rms_difference(experiment, sim) == pytest.approx(rms, abs=2e-4)


def test_table3_o2_model_against_the_experiments():
    for temperature_C, r2, rms in ((25, 0.579, 0.243), (75, 0.644, 0.174),
                                   (100, 0.746, 0.150)):
        experiment = pkg.read_fig10_o2(temperature_C)['experiment']
        sim = simulate_at_times(load_o2_model(T=temperature_C + 273.15),
                                experiment['time'])
        assert pooled_r2(experiment, sim) == pytest.approx(r2, abs=2e-3)
        assert rms_difference(experiment, sim) == pytest.approx(rms, abs=2e-3)


def time_of_half_conversion(times, X_HP3):
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
    # Time at which half of the 3-HP is converted and maximum of the MA
    # fraction: Table 3 O2 model against the per-temperature constant set.
    # Both are interpolated on a 0.1 min grid, so the documented values
    # (NOTES.md) are pinned to their last digit.
    half, peak = {}, {}
    times = np.arange(0., 1890.01, 0.1)
    for temperature_C in TEMPERATURES_C:
        for name, km in (
                ('table3', load_o2_model(T=temperature_C + 273.15)),
                ('set', load_o2_reversible_model(label(temperature_C)))):
            sim = simulate_at_times(km, times, ['X_HP3', 'X_MA'])
            half[name, temperature_C] = time_of_half_conversion(
                times, sim['X_HP3'].to_numpy())
            peak[name, temperature_C] = maximum(times, sim['X_MA'].to_numpy())
    assert [half['table3', T] for T in TEMPERATURES_C] == pytest.approx(
        [458.1, 16.04, 4.16], abs=0.02)
    assert half['set', 25] == np.inf
    assert [half['set', T] for T in (75, 100)] == pytest.approx(
        [75.0, 22.0], abs=0.02)
    for T, table3, per_set in ((25, (0.5218, 1479.1), (0.0278, 374.3)),
                               (75, (0.5842, 51.7), (0.4034, 134.6)),
                               (100, (0.5866, 12.8), (0.5488, 42.2))):
        for name, (value, time) in (('table3', table3), ('set', per_set)):
            assert peak[name, T][0] == pytest.approx(value, abs=1e-4)
            assert peak[name, T][1] == pytest.approx(time, abs=0.06)
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
