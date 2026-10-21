# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Event-aware local sensitivity and the private FD engine (F2 items 2-3;
spec tests a, b, c, e, j, k). Expected values come from closed forms."""
import math
import warnings

import matplotlib
import numpy as np
import pandas as pd
import pytest
import tellurium as te

import nskinetics as nsk
from nskinetics.analysis import (At, CallableMetric, EventCount, Final,
                                 LocalSensitivityResult, Productivity,
                                 TimeAtMax, Yield, local_sensitivity)
from nskinetics.analysis import _fd
from nskinetics.exceptions import KineticSimulationError

matplotlib.use('Agg')

UNITS = {'time': 'h', 'conc': 'g/L'}
TIGHT = (1e-12, 1e-10)          # integrator (atol, rtol) well below rel_step
CVODE = (1e-12, 1e-6)           # CVODE's own default tolerances

# Every time S falls to 5 it is topped up by 5, so between resets
# S = 10 e^(-k (t - t_n)) with t_n = n ln2 / k. With 4 resets before t = 6
# at k = 0.5: S(6) = 10 e^(-6k) 2^4 = 160 e^-3, hence
# d ln S(6) / d ln k = -k t = -3 and dS/dk = -t S = -960 e^-3.
EVENT_TOY = """
model ev()
  S = 10; n = 0; k = 0.5;
  J: S => ; k*S;
  E1: at S < 5: S = S + 5, n = n + 1;
end
"""

# S => P without events: S = S0 e^(-kt), P = P0 + S0 (1 - e^(-kt)).
DECAY_TOY = """
model decay()
  S = 10; P = 0; k = 0.5;
  J: S => P; k*S;
end
"""

FEEDSPIKE_TOY = """
model fb()
  compartment env = 1; species S in env, P in env;
  S = 50; P = 0; k = 0.5; threshold = 5; target = 50; feed_conc = 500;
  n_max = 3; n_spk = 0; last_vol = 0; tot_vol = 0;
  J: S => P; k*S*env;
end
"""


def _feedspike_model():
    km = nsk.KineticModel(te.loadAntimonyModel(FEEDSPIKE_TOY), units=UNITS)
    spike = nsk.FeedSpike('S', 'S <= threshold', 'target', 'feed_conc',
                          volume_var='env', max_count='n_max',
                          count_var='n_spk', last_vol_var='last_vol',
                          tot_vol_var='tot_vol')
    for event in spike.expand():
        km.add_event(event)
    km.compile_events()
    return km


def _model_state(km):
    """Everything A7 says an analysis must leave unchanged."""
    integrator = km._te.getIntegrator()
    return (km._user_parameter_values(),
            {v: km._te[v] for v in km.baseline_variables},
            list(km.events), km._events_compiled,
            {k: integrator.getValue(k) for k in integrator.getSettings()},
            km._te.getSBML())


@pytest.fixture
def close_figures():
    import matplotlib.pyplot as plt
    try:
        yield
    finally:
        plt.close('all')


@pytest.fixture(scope='module')
def feedspike_ls():
    km = _feedspike_model()
    metrics = {'P_titer': Final('[P]'),
               'P_yield': Yield('[P]', '[S]', volume='env',
                                fed=('tot_vol', 'feed_conc'), basis='added'),
               'P_prod': Productivity('[P]', turnaround=12.0),
               'n_spikes': EventCount('n_spk')}
    cond = nsk.Condition(initial={'[S]': 50}, t_end=48, n_steps=97)
    ls = local_sensitivity(km, cond, ['k', 'threshold', 'target', '[S]'],
                           outputs=['[S]', '[P]'], metrics=metrics,
                           tolerances=TIGHT)
    return km, cond, ls


# --- (a) analytic event toy -----------------------------------------------------
def test_event_toy_matches_closed_form_without_flags():
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    k, t = 0.5, 6.0
    ls = local_sensitivity(km, nsk.Condition(t_end=t, n_steps=61), ['k'],
                           outputs=['[S]'], metrics={'S6': Final('[S]')},
                           tolerances=TIGHT)
    assert ls.counters == ('n',)
    assert np.isclose(ls.nominal[-1, 0], 160 * math.exp(-3), rtol=1e-6)
    assert np.isclose(ls.scaled[-1, 0, 0], -k * t, rtol=1e-4)
    assert np.isclose(ls.raw[-1, 0, 0], -960 * math.exp(-3), rtol=1e-4)
    flags = ls.flags.loc['k']
    assert not flags['event_count_changed'] and flags['reliable']
    assert flags['n_refinements'] == 0 and flags['scheme'] == 'central'
    assert flags['error'] is None and flags['step'] == 1e-3
    ms = ls.metric_sensitivities.set_index('metric').loc['S6']
    assert np.isclose(ms['elasticity'], -3.0, rtol=1e-4)
    assert np.isclose(ms['derivative'], -960 * math.exp(-3), rtol=1e-4)


def test_event_toy_default_settings_are_accurate():
    # No tolerances given: the private copy is tightened to (1e-12, 1e-10);
    # the caller's integrator keeps CVODE's defaults.
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=6, n_steps=61), ['k'],
                           outputs=['[S]'])
    assert ls.tolerances == TIGHT and ls.max_refinements == 2
    assert np.isclose(ls.scaled[-1, 0, 0], -3.0, rtol=1e-5)
    assert ls.flags.loc['k', 'reliable']
    assert km._te.getIntegrator().getValue('relative_tolerance') == 1e-6


def test_default_settings_match_the_closed_form_on_a_sum_rate():
    # S => P at rate (k1 + k2) S: dS/dk1 = dS/dk2 = -t S. At CVODE's own
    # tolerances the k1 column is ~1% off; the default tightening fixes it.
    km = nsk.KineticModel.from_antimony("""
    model s()
      S = 10; P = 0; k1 = 0.2; k2 = 0.3;
      J: S => P; (k1 + k2)*S;
    end
    """, units=UNITS)
    cond = nsk.Condition(t_end=6, n_steps=61)
    ls = local_sensitivity(km, cond, ['k1', 'k2'], outputs=['[S]'])
    t = ls.time
    expected = -t * 10 * np.exp(-0.5 * t)
    for j in range(2):
        assert np.allclose(ls.raw[:, 0, j], expected, rtol=1e-5, atol=1e-12)
    # a model whose integrator is already tighter keeps its tolerances
    km._te.getIntegrator().setValue('relative_tolerance', 1e-12)
    km._te.getIntegrator().setValue('absolute_tolerance', 1e-14)
    tighter = local_sensitivity(km, cond, ['k1'], outputs=['[S]'])
    assert tighter.tolerances == (1e-14, 1e-12)


# --- (b) an event-count change is refined and flagged ---------------------------
def test_event_count_change_is_refined_twice_and_flagged():
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    k_edge = 4 * math.log(2) / 6           # the 4th reset lands exactly at t=6
    cond = nsk.Condition(parameters={'k': k_edge}, t_end=6, n_steps=61)
    ls = local_sensitivity(km, cond, ['k'], outputs=['[S]'],
                           metrics={'S6': Final('[S]')}, tolerances=TIGHT)
    flags = ls.flags.loc['k']
    # Whichever way the numerical tie at t=6 resolves, one side changes the
    # count at every step: flagged after two refinements.
    assert flags['event_count_changed'] and not flags['reliable']
    assert flags['n_refinements'] == 2
    # rank() and metric_sensitivities show the flags inline
    row = ls.rank('[S]').iloc[0]
    assert not row['reliable'] and row['event_count_changed']
    assert row['n_refinements'] == 2
    ms = ls.metric_sensitivities.iloc[0]
    assert not ms['reliable'] and ms['event_count_changed']
    assert ms['n_refinements'] == 2


def test_event_count_change_keeps_the_agreeing_side_at_its_largest_step():
    # Just below the edge the 4th reset lands 6e-6 h after t = 6 (and the
    # earlier ones 1.5e-6 to 4.5e-6 h after the output times 1.5, 3, 4.5),
    # so the nominal run has 3 resets by a clear margin: k- keeps them, k+
    # moves them across those times at every step down to 1e-5. The value is
    # the backward difference at the largest agreeing step (1e-3),
    # -6 S(6) = -30 from below up to its O(h) truncation error (~1e-3
    # relative).
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    k = 4 * math.log(2) / 6 * (1 - 1e-6)
    cond = nsk.Condition(parameters={'k': k}, t_end=6, n_steps=61)
    ls = local_sensitivity(km, cond, ['k'], outputs=['[S]'],
                           tolerances=TIGHT)
    flags = ls.flags.loc['k']
    assert flags['event_count_changed'] and not flags['reliable']
    assert flags['n_refinements'] == 2 and flags['scheme'] == 'backward'
    assert flags['step'] == pytest.approx(1e-3, rel=1e-12)
    assert np.isclose(ls.raw[-1, 0, 0], -6 * ls.nominal[-1, 0], rtol=2e-3)


def test_event_counters_override():
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    k_edge = 4 * math.log(2) / 6
    cond = nsk.Condition(parameters={'k': k_edge}, t_end=6, n_steps=61)
    off = local_sensitivity(km, cond, ['k'], outputs=['[S]'],
                            event_counters=())
    assert off.counters == () and off.flags.loc['k', 'reliable']
    assert off.flags.loc['k', 'n_refinements'] == 0
    explicit = local_sensitivity(km, cond, ['k'], outputs=['[S]'],
                                 event_counters=['n'], max_refinements=1)
    assert explicit.counters == ('n',)
    assert explicit.flags.loc['k', 'n_refinements'] == 1
    assert explicit.flags.loc['k', 'event_count_changed']


def test_forward_method_on_the_edge_uses_the_agreeing_side():
    # Just below the edge (see above) k+ adds the 4th reset at every step
    # while the nominal run has 3, so the forward differences are refined
    # twice and fall back to backward runs, made for the whole ladder of
    # steps; the largest agreeing one is used.
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    k = 4 * math.log(2) / 6 * (1 - 1e-6)
    ls = local_sensitivity(km, nsk.Condition(parameters={'k': k},
                                             t_end=6, n_steps=61),
                           ['k'], outputs=['[S]'], method='forward',
                           tolerances=TIGHT)
    flags = ls.flags.loc['k']
    assert flags['event_count_changed'] and not flags['reliable']
    assert flags['n_refinements'] == 2 and flags['scheme'] == 'backward'
    assert flags['step'] == pytest.approx(1e-3, rel=1e-12)
    assert np.isclose(ls.raw[-1, 0, 0], -6 * ls.nominal[-1, 0], rtol=2e-3)


@pytest.mark.parametrize('method', ['central', 'forward'])
def test_edge_fallback_at_cvode_tolerances_is_within_one_percent(method):
    # At CVODE's own tolerances the nominal count at the edge is numerically
    # marginal, so either side may be the agreeing one; whichever is used,
    # the one-sided value is taken at the largest agreeing step and stays
    # within 1% of -6 S(6) (the refined step 1e-5 gave 34% before).
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    cond = nsk.Condition(parameters={'k': 4 * math.log(2) / 6}, t_end=6,
                         n_steps=61)
    ls = local_sensitivity(km, cond, ['k'], outputs=['[S]'], method=method,
                           tolerances=CVODE)
    assert abs(ls.raw[-1, 0, 0] / (-6 * ls.nominal[-1, 0]) - 1) < 0.01


def test_a_firing_moved_across_an_output_time_is_refined_or_flagged():
    # The first reset fires at t_e = ln2/k; a step h in ln k moves it by
    # ~h t_e = 1.4e-3 h. An output time 1e-3 h after t_e is crossed by k- at
    # h = 1e-3 (the final count is unchanged: only the counter column at
    # that row differs), so the step is refined once and the row then
    # matches d ln S/d ln k = -k t.
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    k = 0.5
    t_e = math.log(2) / k
    ls = local_sensitivity(km, nsk.Condition(times=[0, t_e + 1e-3, 6]), ['k'],
                           outputs=['[S]'], tolerances=TIGHT)
    flags = ls.flags.loc['k']
    assert flags['reliable'] and not flags['event_count_changed']
    assert flags['n_refinements'] == 1 and flags['scheme'] == 'central'
    assert np.allclose(ls.scaled[1:, 0, 0], -k * ls.time[1:], rtol=1e-3)
    # 1e-6 h after t_e even the step 1e-5 (a 1.4e-5 h shift) crosses it: the
    # parameter is flagged and the row comes from the side that does not
    # cross (k+ fires earlier), still -k t, not the jump over 2h.
    ls = local_sensitivity(km, nsk.Condition(times=[0, t_e + 1e-6, 6]), ['k'],
                           outputs=['[S]'], tolerances=TIGHT)
    flags = ls.flags.loc['k']
    assert flags['event_count_changed'] and not flags['reliable']
    assert flags['n_refinements'] == 2 and flags['scheme'] == 'forward'
    assert flags['step'] == pytest.approx(1e-3, rel=1e-12)
    assert np.allclose(ls.scaled[1:, 0, 0], -k * ls.time[1:], rtol=2e-3)


def test_fd_counters_are_compared_at_every_output_time():
    # Without refinement the k- run of the grid above has the same final
    # count as the nominal run, but its counter column differs at the row
    # just after t_e: that alone is a count change.
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    private = _fd.private_model(km, tolerances=TIGHT)
    t_e = math.log(2) / 0.5
    cond = nsk.Condition(times=[0, t_e + 1e-3, 6])
    jac = _fd.fd_jacobian(private, [cond], [_fd.FDParameter('k', 'k', 0.5)],
                          lambda rs, values: rs[0]['[S]'], max_refinements=0)
    flags = jac.flags.loc['k']
    assert flags['event_count_changed'] and flags['scheme'] == 'forward'
    nominal = _fd._counter_trace(jac.results[0], ('n',))
    assert nominal[:, 0].tolist() == [0, 1, 4]
    late = private.run(cond.with_parameters(k=0.5 * math.exp(-1e-3)),
                       extra_selections=['n'])
    trace = _fd._counter_trace(late, ('n',))
    assert trace[-1, 0] == nominal[-1, 0] and trace[1, 0] == 0
    assert not _fd._same_trace(trace, nominal)
    assert _fd._same_trace(nominal.copy(), nominal)
    assert _fd._same_trace(None, None)


def test_time_at_max_has_a_meaningful_local_derivative():
    # A -> B -> C: B peaks at t* = ln(k1/k2)/(k1 - k2). TimeAtMax refines
    # the grid maximum by a parabola, so its finite difference is the
    # derivative of t*, not 0 or a grid jump; the elasticities sum to -1
    # (scaling both rates by c scales t* by 1/c).
    km = nsk.KineticModel.from_antimony("""
    model chain()
      A = 1; B = 0; C = 0; k1 = 1; k2 = 0.25;
      J1: A => B; k1*A;
      J2: B => C; k2*B;
    end
    """, units=UNITS)
    k1, k2 = 1.0, 0.25
    t_star = math.log(k1 / k2) / (k1 - k2)
    dt_dk1 = ((k1 - k2) / k1 - math.log(k1 / k2)) / (k1 - k2) ** 2
    e1 = k1 / t_star * dt_dk1                       # -0.6120
    ls = local_sensitivity(km, nsk.Condition(t_end=10, n_steps=1001),
                           ['k1', 'k2'], outputs=[],
                           metrics={'tpk': TimeAtMax('[B]')})
    assert np.isclose(ls.metric_values[0], t_star, rtol=1e-4)
    ms = ls.metric_sensitivities.set_index('parameter')
    assert np.isclose(ms.loc['k1', 'elasticity'], e1, rtol=1e-2)
    assert np.isclose(ms.loc['k2', 'elasticity'], -1 - e1, rtol=1e-2)
    assert np.isclose(ms['elasticity'].sum(), -1, atol=1e-2)
    assert ms['reliable'].all()


def test_tightened_tolerances_raise_the_step_budget_and_name_failures():
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    integrator = km._te.getIntegrator()
    private = _fd.private_model(km)               # rtol 1e-6 -> 1e-10
    budget = private._te.getIntegrator().getValue('maximum_num_steps')
    assert budget == 10 * integrator.getValue('maximum_num_steps')
    same = _fd.private_model(km, tolerances=(1e-12, 1e-6))
    assert same._te.getIntegrator().getValue('maximum_num_steps') == \
        integrator.getValue('maximum_num_steps')
    # A budget too small even after the raise: the start-value run (k is not
    # set by the condition) and the nominal run both fail, naming the
    # tolerances in use and the remedies; the caller's model is unchanged.
    integrator.setValue('maximum_num_steps', 2)
    before = _model_state(km)
    with pytest.raises(KineticSimulationError,
                       match=r'start-value run.*tolerances \(atol, rtol\) = '
                             r'\(1e-12, 1e-10\).*tolerances='):
        local_sensitivity(km, nsk.Condition(t_end=100, n_steps=3), ['k'])
    with pytest.raises(KineticSimulationError,
                       match=r'nominal run.*tolerances='):
        local_sensitivity(km, nsk.Condition(t_end=100, n_steps=3,
                                            parameters={'k': 0.5}), ['k'])
    assert _model_state(km) == before


# --- (c) no-event decay ---------------------------------------------------------
def test_decay_log_sensitivity_is_minus_kt_at_every_time():
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=6, n_steps=61),
                           ['k', '[S]'], tolerances=TIGHT)
    assert ls.outputs == ['[S]', '[P]'] and ls.counters == ()
    k, t = 0.5, ls.time
    assert np.allclose(ls.scaled[:, 0, 0], -k * t, rtol=1e-5, atol=0)
    # initial value: d ln S(t) / d ln S0 = 1, dP/dS0 = 1 - e^(-kt)
    assert np.allclose(ls.scaled[:, 0, 1], 1.0, rtol=1e-5)
    assert np.allclose(ls.raw[:, 1, 1], 1 - np.exp(-k * t), rtol=1e-5,
                       atol=1e-9)
    # P(0) = 0: the log scaling is undefined there
    assert np.isnan(ls.scaled[0, 1, 0]) and not np.isnan(ls.raw[0, 1, 0])
    assert ls.flags['reliable'].all()


def test_scalings_semi_and_none():
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    cond = nsk.Condition(t_end=6, n_steps=61)
    semi = local_sensitivity(km, cond, ['k'], outputs=['[S]'], scaled='semi',
                             tolerances=TIGHT)
    none = local_sensitivity(km, cond, ['k'], outputs=['[S]'], scaled='none',
                             tolerances=TIGHT)
    k, t = 0.5, semi.time
    s = 10 * np.exp(-k * t)
    assert np.allclose(semi.scaled[:, 0, 0], -k * t * s, rtol=1e-5, atol=1e-12)
    assert np.allclose(none.scaled[:, 0, 0], -t * s, rtol=1e-5, atol=1e-12)
    assert np.array_equal(none.raw, none.scaled)
    assert semi.rank('[S]').columns[3] == 'semi_elasticity'


def test_log_sensitivities_at_integrator_noise_level_are_nan(tmp_path):
    # S = 10 e^(-t) reaches ~1e-25: there its difference is integration
    # noise (~3 atol / h), so d ln S / d ln k2 (exactly 0) came out O(1).
    km = nsk.KineticModel.from_antimony("""
    model d2()
      S = 10; P = 0; k = 1; k2 = 0.3;
      J: S => P; k*S;
      J2: P => ; k2*P;
    end
    """, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=60, n_steps=61),
                           ['k', 'k2'], outputs=['[S]', '[P]'])
    s, t = ls.nominal[:, 0], ls.time
    # RoadRunner scales atol by the start value S0 = 10: 30 atol S0 / h
    cut = _fd.FD_NOISE * ls.tolerances[0] * 10 / (0.1 * 1e-3)     # 3e-7
    noisy, clean = s < cut, s > 1e-3
    assert noisy.sum() > 30 and clean.sum() > 5
    assert np.isnan(ls.scaled[noisy, 0, :]).all()
    assert np.isfinite(ls.scaled[s > 1.01 * cut, 0, :]).all()
    with np.errstate(all='ignore'):
        direct = ls.raw[:, 0, :] * ls.values / s[:, None]
    assert np.array_equal(ls.scaled[clean, 0, :], direct[clean])
    assert np.allclose(ls.scaled[clean, 0, 0], -t[clean], rtol=1e-5)
    assert np.allclose(ls.scaled[clean, 0, 1], 0.0, atol=1e-5)
    # raw derivatives, other scalings and other outputs are untouched
    assert np.isfinite(ls.raw).all() and not np.isnan(ls.scaled[1:, 1]).any()
    assert ls.time_course('[S]').iloc[-1].isna().all()
    assert ls.rank('[S]')['elasticity'].isna().all()
    assert ls.rank('[S]', at=2.0)['elasticity'].notna().all()
    semi = local_sensitivity(km, nsk.Condition(t_end=60, n_steps=61),
                             ['k2'], outputs=['[S]'], scaled='semi')
    assert np.isfinite(semi.scaled).all()
    # the mask comes from the recorded tolerances and steps: exact round trip
    path = tmp_path / 'ls.csv'
    ls.to_csv(path)
    back = LocalSensitivityResult.from_csv(path)
    assert np.array_equal(back.scaled, ls.scaled, equal_nan=True)


def _scaled_decay(S0, units=None):
    return nsk.KineticModel.from_antimony(f"""
    model d()
      S = {S0!r}; P = 0; k = 0.5;
      J: S => P; k*S;
    end
    """, units=units or {'time': 'h', 'conc': 'mol/L'})


def test_log_sensitivity_noise_level_is_relative_to_the_start_value():
    # RoadRunner multiplies a scalar absolute tolerance by each state's start
    # value, so a 1 uM substrate in mol/L is as accurate as 10 g/L: its log
    # sensitivities (exactly -k t) were masked as NaN by an absolute cut-off.
    assert np.array_equal(
        _fd.effective_atol(1e-12, [10.0, 0.0, -2.0, np.nan]),
        [1e-11, 1e-12, 2e-12, 1e-12])
    micro = local_sensitivity(_scaled_decay(1e-6),
                              nsk.Condition(t_end=10, n_steps=11),
                              ['k', '[S]'], outputs=['[S]', '[P]'])
    t = micro.time
    assert micro.nominal[-1, 0] < 1e-8                  # below 30 atol / h
    assert np.isfinite(micro.scaled[:, 0, :]).all()
    assert np.allclose(micro.scaled[:, 0, 0], -0.5 * t, rtol=1e-4, atol=1e-9)
    assert np.allclose(micro.scaled[:, 0, 1], 1.0, rtol=1e-4)
    rank = micro.rank('[S]')
    assert rank['elasticity'].notna().all()
    assert rank['parameter'].iloc[0] == 'k'
    assert rank['elasticity'].iloc[0] == pytest.approx(-5.0, rel=1e-4)
    # [P] starts at 0 and peaks at ~1e-6: the default tolerances scale atol
    # to that peak, and the floors use the recorded (scaled) atol.
    atol = micro.tolerances[0]
    assert atol == pytest.approx(1e-12 * micro.nominal[:, 1].max(),
                                 rel=1e-6)
    floor = micro._log_floor()                 # per (output, parameter)
    assert floor.shape == (2, 2)
    assert np.allclose(floor[0], 30 * atol * 1e-6 / 1e-3)      # atol S0
    assert np.allclose(floor[1], 30 * atol / 1e-3)     # [P] starts at 0
    # The other way round: from S0 = 1e4, S reaches ~3e-8 at t = 54, where
    # its difference is integration noise (once +20 instead of -27).
    large = local_sensitivity(_scaled_decay(1e4, {'time': 'h',
                                                  'conc': 'mg/L'}),
                              nsk.Condition(t_end=60, n_steps=31), ['k'],
                              outputs=['[S]'])
    t, s = large.time, large.nominal[:, 0]
    assert np.isnan(large.scaled[t >= 54, 0, 0]).all()
    early = s / 1e4 > 1e-6
    assert early.sum() > 10
    assert np.allclose(large.scaled[early, 0, 0], -0.5 * t[early],
                       rtol=1e-4, atol=1e-9)
    assert large.rank('[S]')['elasticity'].isna().all()


DECOUPLED_TOY = """
model m()
  B = 0; v = 1e-9; k = 1;
  J1: => B; v;
  J2: B => ; k*B;
end
"""


def test_default_tolerances_scale_atol_to_outputs_starting_at_zero():
    # B starts at 0 and nothing else sets the step size: RoadRunner keeps a
    # fixed absolute tolerance atol for it, so at atol = 1e-12 a nanomolar B
    # is 10% off in d ln B / d ln k (and was masked as NaN). The default
    # tolerances scale atol to B's peak; explicit ones are used as given.
    km = nsk.KineticModel.from_antimony(DECOUPLED_TOY,
                                        units={'time': 'h', 'conc': 'mol/L'})
    cond = nsk.Condition(t_end=10, n_steps=21)
    ls = local_sensitivity(km, cond, ['k'], outputs=['[B]'])
    assert ls.tolerances[0] == pytest.approx(1e-12 * ls.nominal[:, 0].max(),
                                             rel=1e-6)
    assert ls.tolerances[1] == 1e-10
    t = ls.time[1:]
    exact = -1 + t * np.exp(-t) / (1 - np.exp(-t))     # k = 1
    assert np.isfinite(ls.scaled[1:, 0, 0]).all()
    assert np.allclose(ls.scaled[1:, 0, 0], exact, rtol=0, atol=1e-4)
    tight = local_sensitivity(km, cond, ['k'], outputs=['[B]'],
                              tolerances=TIGHT)
    assert tight.tolerances == TIGHT
    assert np.isnan(tight.scaled[1:, 0, 0]).all()
    # The columns a metric reads are scaled too: a metric-only call (the
    # kinetic leg E_m,k) once kept atol = 1e-12 and was 1.2e-2 off.
    m = local_sensitivity(km, cond, ['k'], outputs=[],
                          metrics={'Bend': Final('[B]')})
    assert m.tolerances == pytest.approx(ls.tolerances, rel=1e-6)
    assert m.tolerances[0] < 1e-12
    elasticity = m.metric_sensitivities.set_index('metric').loc[
        'Bend', 'elasticity']
    assert elasticity == pytest.approx(
        -1 + 10 * math.exp(-10) / (1 - math.exp(-10)), abs=1e-4)
    both = local_sensitivity(km, cond, ['k'], outputs=['[B]'],
                             metrics={'Bend': Final('[B]')})
    assert both.metric_sensitivities['elasticity'].iloc[0] == pytest.approx(
        elasticity, abs=1e-5)
    # The helper itself (for calibration's several conditions): one run per
    # condition; outputs that start nonzero, stay 0 or peak above 1 change
    # nothing, and the caller's model is not touched.
    private = _fd.private_model(km)
    assert _fd.scale_atol_to_outputs(private, [cond], []) == TIGHT
    assert _fd.scale_atol_to_outputs(
        private, [cond.with_initial({'[B]': 1e-9})], ['[B]']) == TIGHT
    assert _fd.scale_atol_to_outputs(
        private, [cond.with_parameters(v=0.0)], ['[B]']) == TIGHT
    assert _fd.scale_atol_to_outputs(
        private, [cond.with_parameters(v=5.0)], ['[B]']) == TIGHT
    atol, rtol = _fd.scale_atol_to_outputs(
        private, [cond.with_parameters(v=1e-3), cond], ['[B]', '[B]'])
    assert rtol == 1e-10 and atol == pytest.approx(
        1e-12 * ls.nominal[:, 0].max(), rel=1e-6)
    assert _fd.tolerances_of(private) == (atol, rtol)
    assert _fd.tolerances_of(km) == (1e-12, 1e-6)


ZERO_STEP_WARNING = r"\['\[P\]'\] are at 0.*rel_step=0.001 in their own units"


def test_zero_parameter_is_differenced_in_linear_space():
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    cond = nsk.Condition(t_end=6, n_steps=61)
    # The default step at p = 0 (rel_step in the id's own units) is
    # unrelated to the model's scale: a warning says so.
    with pytest.warns(UserWarning, match=ZERO_STEP_WARNING) as record:
        ls = local_sensitivity(km, cond, ['[P]'], outputs=['[P]'],
                               tolerances=TIGHT)
    assert len(record) == 1 and 'abs_step={id: step}' in str(
        record[0].message)
    flags = ls.flags.loc['[P]']
    assert flags['value'] == 0.0 and not flags['log']
    assert flags['step'] == 1e-3                 # rel_step, absolute at p=0
    assert np.allclose(ls.raw[:, 0, 0], 1.0, rtol=1e-5)   # dP(t)/dP0 = 1
    assert np.all(ls.scaled[1:, 0, 0] == 0.0)             # p = 0
    with warnings.catch_warnings():
        warnings.simplefilter('error')          # an abs_step entry: silent
        custom = local_sensitivity(km, cond, ['[P]', 'k'], outputs=['[P]'],
                                   abs_step={'[P]': 0.01, 'k': 0.005},
                                   tolerances=TIGHT)
        local_sensitivity(km, cond, ['[P]'], outputs=['[P]'],
                          abs_step=0.01, tolerances=TIGHT)    # a scalar too
    assert custom.flags.loc['[P]', 'step'] == 0.01
    assert not custom.flags.loc['k', 'log']        # listed: linear space
    assert custom.flags.loc['k', 'step'] == 0.005
    t = custom.time
    # dP/dk = 10 t e^(-kt); central truncation error h^2/6 |d3P/dk3| with
    # the coarse absolute step h = 0.005 is below 2e-4 relative
    assert np.allclose(custom.raw[:, 0, 1], 10 * t * np.exp(-0.5 * t),
                       rtol=1e-3, atol=1e-9)


def test_default_step_at_zero_is_unrelated_to_the_model_scale():
    # Product inhibition at micromolar levels in mol/L: the default step of
    # [P] = 0 is 1e-3 mol/L, 10 x S0 and 100 x Ki (the minus side runs the
    # denominator at -99), and dS/dP0 comes out reliable but with the wrong
    # sign; abs_step on the quantity's scale gives the derivative.
    km = nsk.KineticModel.from_antimony("""
    model pinh()
      S = 1e-4; P = 0; k = 0.5; Ki = 1e-5;
      J: S => P; k*S/(1 + P/Ki);
    end
    """, units={'time': 'h', 'conc': 'mol/L'})
    cond = nsk.Condition(t_end=10, n_steps=11)
    with pytest.warns(UserWarning, match=ZERO_STEP_WARNING):
        coarse = local_sensitivity(km, cond, ['[P]'], outputs=['[S]'])
    assert coarse.flags.loc['[P]', 'reliable']
    assert (coarse.raw[1:, 0, 0] < 0).all()                 # wrong sign
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        fine = local_sensitivity(km, cond, ['[P]'], outputs=['[S]'],
                                 abs_step={'[P]': 1e-8})
    assert (fine.raw[1:, 0, 0] > 0).all()
    assert fine.raw[1, 0, 0] == pytest.approx(0.60, abs=0.01)   # t = 1


# --- (e) the caller's model is never mutated ------------------------------------
def test_user_model_unchanged_and_base_condition_kept():
    km = _feedspike_model()
    km._te['k'] = 0.4                          # a user edit must persist
    before = _model_state(km)
    cond = nsk.Condition(initial={'[S]': 40}, parameters={'threshold': 10},
                         t_end=48, n_steps=97)
    ls = local_sensitivity(km, cond, ['k', 'threshold', '[S]', 'env'],
                           outputs=['[S]', '[P]'],
                           metrics=[Final('[P]'), EventCount('n_spk')],
                           tolerances=TIGHT)
    assert _model_state(km) == before
    assert km._te['k'] == 0.4
    # nominal values: the condition's own (threshold, [S]) or the model's
    assert dict(zip(ls.parameters, ls.values)) == {
        'k': 0.4, 'threshold': 10.0, '[S]': 40.0, 'env': 1.0}
    # the nominal run is the base condition's run (up to the tolerances)
    ref = km.run(cond, extra_selections=['n_spk'])
    assert np.allclose(ls.nominal, ref[['[S]', '[P]']], rtol=1e-4, atol=1e-6)
    assert ls.counters == ('n_spk',)
    assert ls.metric_values[1] == ref.final['n_spk']


def test_invalid_arguments_raise_before_any_run(monkeypatch):
    km = _feedspike_model()
    before = _model_state(km)
    cond = nsk.Condition(t_end=48, n_steps=97)
    runs = []
    original = nsk.KineticModel.run

    def counting_run(self, *args, **kwargs):
        runs.append(1)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(nsk.KineticModel, 'run', counting_run)
    bad = [
        (dict(params=['nope']), ValueError, '^params: '),
        (dict(params=[]), ValueError, None),
        (dict(params=['k', 'k']), ValueError, None),
        (dict(params='k', outputs=['time']), ValueError, None),
        (dict(params=['k'], outputs=['[nope]']), ValueError, '^outputs: '),
        (dict(params=['k'], outputs=[]), ValueError, None),
        (dict(params=['k'], metrics={'m': Final('[nope]')}), ValueError,
         '^metrics: '),
        (dict(params=['k'], method='backward'), ValueError, None),
        (dict(params=['k'], scaled='ln'), ValueError, None),
        (dict(params=['k'], rel_step=0), ValueError, None),
        (dict(params=['k'], abs_step={'other': 1}), ValueError, None),
        (dict(params=['k'], tolerances=(1e-12,)), TypeError, None),
        (dict(params=['k'], event_counters=['nope']), ValueError,
         '^event_counters: '),
        (dict(params=['k'], workers=0), ValueError, None),
        (dict(params=['k'], max_refinements='2'), TypeError,
         'max_refinements'),
        (dict(params=['k'], max_refinements=-1), ValueError,
         'max_refinements'),
        (dict(params=['[S]', 'S']), ValueError, '^params: .*once'),
        (dict(params=['k'], metrics={'late': At('[S]', 100)}), ValueError,
         r'^metrics: .*outside the simulated range \[0.0, 48.0\]'),
    ]
    for kwargs, error, match in bad:
        params = kwargs.pop('params')
        with pytest.raises(error, match=match):
            local_sensitivity(km, cond, params, **kwargs)
    # a species listed in the other form than the condition's own value
    with pytest.raises(ValueError, match='^params: .*same form'):
        local_sensitivity(km, cond.with_initial({'[S]': 5}), ['S'])
    with pytest.raises(ValueError, match='^metrics: .*outside'):
        local_sensitivity(km, nsk.Condition(times=[1, 2, 3]), ['k'],
                          metrics={'early': At('[S]', 0.5)})
    with pytest.raises(ValueError, match='reset=False'):
        local_sensitivity(km, cond.replace(reset=False), ['k'])
    with pytest.raises(TypeError):
        local_sensitivity(km, {'t_end': 1}, ['k'])
    assert runs == []
    assert _model_state(km) == before


def test_tolerance_warning():
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    cond = nsk.Condition(t_end=1, n_steps=11)
    with pytest.warns(UserWarning, match='relative tolerance'):
        local_sensitivity(km, cond, ['k'], rel_step=1e-5, tolerances=CVODE)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        # the default tightening (rtol 1e-10) suits rel_step=1e-5
        ls = local_sensitivity(km, cond, ['k'], rel_step=1e-5)
        local_sensitivity(km, cond, ['k'], rel_step=1e-5,
                          tolerances=(1e-14, 1e-8))
    assert ls.tolerances == TIGHT
    assert km._te.getIntegrator().getValue('relative_tolerance') == 1e-6


def test_refined_steps_warn_when_the_tolerance_is_too_loose():
    # Just above the edge (nominal: 4 resets) the k- run loses the 4th reset
    # at step 1e-3 but keeps it at 1e-4: resolved by one refinement, at a
    # step where rtol 2e-6 exceeds step/100 (the initial 1e-3 passed).
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    cond = nsk.Condition(parameters={'k': 4 * math.log(2) / 6 * (1 + 5e-4)},
                         t_end=6, n_steps=61)
    with pytest.warns(UserWarning, match=r"\['k'\] use a step \(refined"):
        loose = local_sensitivity(km, cond, ['k'], outputs=['[S]'],
                                  tolerances=(1e-12, 2e-6))
    assert loose.flags.loc['k', 'n_refinements'] == 1
    assert loose.flags.loc['k', 'reliable']
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        ls = local_sensitivity(km, cond, ['k'], outputs=['[S]'])
    assert ls.flags.loc['k', 'step'] == pytest.approx(1e-4, rel=1e-12)
    assert np.isclose(ls.raw[-1, 0, 0], -6 * ls.nominal[-1, 0], rtol=1e-5)
    # the helper behind it
    assert _fd.noisy_parameters(ls.flags, (1e-12, 2e-6)) == ['k']
    assert _fd.noisy_parameters(ls.flags, None) == []


def test_failed_nominal_run_raises():
    km = nsk.KineticModel.from_antimony("""
    model blow()
      species S; S = 1; k = 1;
      J: => S; k*S^2;
    end
    """, units=UNITS)
    with pytest.raises(KineticSimulationError, match='nominal run'):
        local_sensitivity(km, nsk.Condition(t_end=2, n_steps=21,
                                            parameters={'k': 1.0}),
                          ['k'])


# --- rank, time_course, workers -------------------------------------------------
def test_rank_and_time_course(feedspike_ls):
    _, _, ls = feedspike_ls
    table = ls.rank('P_titer')
    assert list(table.columns) == [
        'parameter', 'value', 'derivative', 'semi_elasticity', 'elasticity',
        'reliable', 'event_count_changed', 'n_refinements']
    mag = table['elasticity'].abs().to_numpy()
    assert np.all(np.diff(mag) <= 0)
    assert set(table['parameter']) == {'k', 'threshold', 'target', '[S]'}
    # elasticity = derivative * p / m
    m = ls.metric_values[ls.metric_names.index('P_titer')]
    assert np.allclose(table['elasticity'],
                       table['derivative'] * table['value'] / m)
    # an output at a grid time equals the stored row
    t = ls.time[40]
    at = ls.rank('[P]', at=t).set_index('parameter')
    j = ls.parameters.index('k')
    assert at.loc['k', 'derivative'] == ls.raw[40, 1, j]
    assert at.loc['k', 'elasticity'] == pytest.approx(ls.scaled[40, 1, j],
                                                      rel=1e-12)
    with pytest.raises(KeyError):
        ls.rank('[nope]')
    with pytest.raises(ValueError):
        ls.rank('P_titer', at=1.0)
    with pytest.raises(ValueError):
        ls.rank('[P]', at=100.0)
    course = ls.time_course('[P]')
    assert list(course.columns) == ls.parameters
    assert np.array_equal(course.index.to_numpy(), ls.time)
    assert np.array_equal(course.to_numpy(), ls.scaled[:, 1, :],
                          equal_nan=True)


def test_feedspike_sensitivities_have_the_expected_signs(feedspike_ls):
    _, _, ls = feedspike_ls
    assert ls.counters == ('n_spk',) and ls.flags['reliable'].all()
    titer = ls.metric_sensitivities.set_index(['metric', 'parameter'])
    # more feed per spike (target) and more starting substrate raise the
    # titer; a higher threshold leaves more substrate unconverted
    assert titer.loc[('P_titer', 'target'), 'derivative'] > 0
    assert titer.loc[('P_titer', '[S]'), 'derivative'] > 0
    assert titer.loc[('P_titer', 'threshold'), 'derivative'] < 0
    # the spike count does not move under small perturbations
    assert np.all(titer.loc['n_spikes', 'derivative'] == 0.0)


def test_two_workers_match_serial(feedspike_ls):
    km, cond, serial = feedspike_ls
    pooled = local_sensitivity(km, cond, serial.parameters,
                               outputs=serial.outputs,
                               metrics={'P_titer': Final('[P]')},
                               tolerances=TIGHT, workers=2)
    assert np.allclose(pooled.raw, serial.raw, rtol=1e-9, atol=1e-12,
                       equal_nan=True)
    assert pooled.flags.equals(serial.flags)


# --- (j) CSV round trip ---------------------------------------------------------
def test_csv_round_trip_is_exact(feedspike_ls, tmp_path):
    _, cond, ls = feedspike_ls
    path = tmp_path / 'ls.csv'
    ls.to_csv(path)
    back = LocalSensitivityResult.from_csv(path)
    for name in ('time', 'values', 'nominal', 'raw', 'scaled',
                 'metric_values', 'metric_derivatives'):
        a, b = getattr(ls, name), getattr(back, name)
        assert a.shape == b.shape and np.array_equal(a, b, equal_nan=True), name
    assert np.isnan(back.scaled[0, 1]).all()          # P(0) = 0
    pd.testing.assert_frame_equal(back.flags, ls.flags, check_exact=True)
    pd.testing.assert_frame_equal(back.metric_sensitivities,
                                  ls.metric_sensitivities, check_exact=True)
    for name in ('outputs', 'parameters', 'metric_names', 'scaling',
                 'method', 'rel_step', 'counters', 'units', 'condition',
                 'tolerances', 'max_refinements'):
        assert getattr(back, name) == getattr(ls, name), name
    assert back.condition == cond and back.tolerances == TIGHT
    ms = back.metric_sensitivities
    assert list(ms.columns) == [
        'metric', 'parameter', 'metric_value', 'derivative',
        'semi_elasticity', 'elasticity', 'reliable', 'event_count_changed',
        'n_refinements']
    # semi_elasticity = p dm/dp, as in rank()
    row = ms.iloc[1]
    value = ls.values[ls.parameters.index(row['parameter'])]
    assert row['semi_elasticity'] == row['derivative'] * value
    assert np.allclose(ms['semi_elasticity'],
                       ms['elasticity'] * ms['metric_value'], rtol=1e-12,
                       equal_nan=True)
    back.to_csv(tmp_path / 'again.csv')
    assert (tmp_path / 'again.csv').read_text() == path.read_text()


def test_csv_round_trip_with_flags_and_errors(tmp_path):
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    cond = nsk.Condition(parameters={'k': 4 * math.log(2) / 6}, t_end=6,
                         n_steps=11)
    ls = local_sensitivity(km, cond, ['k', '[S]'], outputs=['[S]'],
                           scaled='semi')
    ls.flags.loc['[S]', 'error'] = 'condition 0: something failed'
    ls.to_csv(tmp_path / 'f.csv')
    back = LocalSensitivityResult.from_csv(tmp_path / 'f.csv')
    pd.testing.assert_frame_equal(back.flags, ls.flags, check_exact=True)
    assert back.scaling == 'semi' and back.metric_names == []
    assert np.array_equal(back.scaled, ls.scaled, equal_nan=True)


def test_from_csv_rejects_other_files(tmp_path):
    path = tmp_path / 'x.csv'
    path.write_text('a,b\n1,2\n')
    with pytest.raises(ValueError, match='not written'):
        LocalSensitivityResult.from_csv(path)


def test_from_csv_rejects_incomplete_files(feedspike_ls, tmp_path):
    # A truncated or edited file is refused, never filled with
    # uninitialized memory (or a KeyError for a missing meta key).
    _, _, ls = feedspike_ls
    path = tmp_path / 'ls.csv'
    ls.to_csv(path)
    lines = path.read_text().splitlines(keepends=True)
    for prefix, match in (('raw,', 'raw: expected'),
                          ('nominal,', 'nominal: expected'),
                          ('time,', 'time: expected'),
                          ('flag,', 'flag: expected'),
                          ('metric_value,', 'metric_value: expected'),
                          ('meta,tolerances,', r"meta keys \['tolerances'\]")):
        first = next(line for line in lines if line.startswith(prefix))
        cut = tmp_path / 'cut.csv'
        cut.write_text(''.join(line for line in lines if line != first))
        with pytest.raises(ValueError, match=f'incomplete.*{match}'):
            LocalSensitivityResult.from_csv(cut)


# --- (k) plots ------------------------------------------------------------------
def test_plot_smoke(feedspike_ls, close_figures):
    _, _, ls = feedspike_ls
    fig, ax = ls.plot_time_course('[P]', top=2)
    assert len(ax.get_lines()) == 3                 # 2 curves + the zero line
    assert ax.get_ylabel() == 'd ln [P] / d ln p'
    fig, ax = ls.plot_heatmap()
    image = ax.get_images()[0]
    assert image.get_array().shape == (2 + 4, 4)    # outputs+metrics x params
    fig, ax = ls.plot_heatmap(t=24.0)
    # Drawing on given axes returns them with their figure.
    fig0, ax0 = matplotlib.pyplot.subplots()
    assert ls.plot_time_course('[P]', ax=ax0) == (fig0, ax0)
    fig0, ax0 = matplotlib.pyplot.subplots()
    assert ls.plot_heatmap(ax=ax0) == (fig0, ax0)
    with pytest.raises(ValueError):
        ls.plot_time_course('[P]', top=0)
    with pytest.raises(TypeError, match='top'):
        ls.plot_time_course('[P]', top='3')
    with pytest.raises(TypeError, match='^t must be a time'):
        ls.plot_heatmap(t='x')
    # The messages name the argument of the method called.
    with pytest.raises(ValueError, match=r'^t=100'):
        ls.plot_heatmap(t=100)
    with pytest.raises(ValueError, match=r'^at=100'):
        ls.rank('[P]', at=100)
    # NaN is not a time at all, rather than outside the range.
    with pytest.raises(ValueError, match='^at must be a finite time'):
        ls.rank('[P]', at=float('nan'))
    with pytest.raises(ValueError, match='^t must be a finite time'):
        ls.plot_heatmap(t=float('nan'))


def test_unreliable_parameters_are_marked_in_plots(close_figures):
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    cond = nsk.Condition(parameters={'k': 4 * math.log(2) / 6}, t_end=6,
                         n_steps=61)
    ls = local_sensitivity(km, cond, ['k'], outputs=['[S]'])
    fig, ax = ls.plot_time_course('[S]')
    assert ax.get_lines()[0].get_label() == 'k (unreliable)'
    assert ax.get_lines()[0].get_linestyle() == '--'
    fig, ax = ls.plot_heatmap()
    assert [t.get_text() for t in ax.get_xticklabels()] == ['k*']


# --- the private FD engine (F3's entry point) -----------------------------------
def test_detect_event_counters():
    assert _fd.detect_event_counters(_feedspike_model()) == ('n_spk',)
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    assert _fd.detect_event_counters(km) == ('n',)     # not 'S = S + 5'
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    assert _fd.detect_event_counters(km) == ()
    km = nsk.KineticModel.from_antimony("""
    model c()
      x = 0; m = 0; S = 10;
      J: S => ; S;
      E1: at time > 1: m = 1 + m;
      E2: at time > 2: x = x + 2;
    end
    """, units=UNITS)
    assert _fd.detect_event_counters(km) == ('m',)


GROWING_TOY = """
model g()
  compartment V = 1; species S in V, n in V;
  S = 10; n = 0; k = 0.5; q = 0.1;
  V' = q;
  J: S => ; k*S;
  E1: at time > 1: n = n + 1;
end
"""


def test_a_species_counter_in_a_growing_compartment_is_not_watched():
    # n = n + 1 adds 1 to the concentration of n, so its amount (the plain
    # id a counter is read as) is V(t_fire) = 1.1, off by round-off between
    # runs: every parameter came out refined and unreliable.
    km = nsk.KineticModel.from_antimony(GROWING_TOY, units=UNITS)
    with pytest.warns(UserWarning, match=r"\['n'\] are incremented.*not "
                                         'watched'):
        assert _fd.detect_event_counters(km) == ()
    with pytest.warns(UserWarning, match=r"\['n'\] are incremented"):
        ls = local_sensitivity(km, nsk.Condition(t_end=4, n_steps=41),
                               ['q', 'k'], outputs=['[S]'])
    assert ls.counters == ()
    assert ls.flags['reliable'].all()
    assert (ls.flags['n_refinements'] == 0).all()
    # In a constant compartment the same species is a counter, read as its
    # concentration (which counts firings whatever the compartment's size).
    km = nsk.KineticModel.from_antimony(
        GROWING_TOY.replace("V' = q;", ''), units=UNITS)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        assert _fd.detect_event_counters(km) == ('[n]',)


def test_species_counter_in_a_perturbed_constant_compartment():
    # n = n + 1 adds 1 to [n]; its amount is n x V. Perturbing V once moved
    # the amount trace by count x dV although no firing changed: V came out
    # refined twice, flagged and unreliable. dS/dV is the same as with n a
    # parameter.
    toy = """
    model c()
      compartment V = 2; species S in V; species n in V;
      S = 10; n = 0; k = 0.5;
      J: S => ; k*S*V;
      E1: at S < 5: S = S + 5, n = n + 1;
    end
    """
    km = nsk.KineticModel.from_antimony(toy, units=UNITS)
    cond = nsk.Condition(t_end=6, n_steps=13)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        ls = local_sensitivity(km, cond, ['V', 'k'], outputs=['[S]'])
    assert ls.counters == ('[n]',)
    assert ls.flags['reliable'].all()
    assert not ls.flags['event_count_changed'].any()
    assert (ls.flags['n_refinements'] == 0).all()
    ref_km = nsk.KineticModel.from_antimony(
        toy.replace('species n in V;', ''), units=UNITS)
    ref = local_sensitivity(ref_km, cond, ['V', 'k'], outputs=['[S]'])
    assert ref.counters == ('n',)
    assert np.allclose(ls.raw, ref.raw, rtol=1e-6, atol=1e-9)
    # A firing at t_end ([S] decays at rate k, so the fourth top-up lands at
    # t = 6 for k = 4 ln2 / 6) is still caught through '[n]'.
    edge = local_sensitivity(km, cond.with_parameters(k=4 * math.log(2) / 6),
                             ['k'], outputs=['[S]'])
    assert edge.flags.loc['k', 'event_count_changed']
    assert not edge.flags.loc['k', 'reliable']


def test_perturbing_a_counter_start_value_is_no_count_change():
    # The firings are compared, not the counter values: n0 +- step shifts
    # the whole trace but fires nothing more.
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    cond = nsk.Condition(t_end=6, n_steps=61)
    with pytest.warns(UserWarning, match=r"\['n'\] are at 0"):
        ls = local_sensitivity(km, cond, ['n', 'k'], outputs=['[S]'])
    assert ls.counters == ('n',)
    flags = ls.flags.loc['n']
    assert flags['reliable'] and not flags['event_count_changed']
    assert flags['n_refinements'] == 0 and flags['scheme'] == 'central'
    assert np.allclose(ls.raw[:, 0, 0], 0.0, atol=1e-9)
    assert ls.flags.loc['k', 'reliable']
    assert ls.scaled[-1, 0, 1] == pytest.approx(-3.0, rel=1e-4)
    # A moved firing is still caught (case b): only the start offset goes.
    trace = np.array([[0.0], [0.0], [1.0]])
    assert _fd._same_trace(trace + 0.001, trace)
    assert not _fd._same_trace(np.array([[0.0], [1.0], [1.0]]), trace)


def test_private_model_caches_the_event_counters(monkeypatch):
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    private = _fd.private_model(km)
    assert private._fd_counters == _fd.detect_event_counters(km) == ('n',)
    assert private._fd_unwatched == ()

    def boom(*args, **kwargs):
        raise AssertionError('event counters re-detected')

    monkeypatch.setattr(_fd, 'detect_event_counters', boom)
    cond = nsk.Condition(t_end=6, n_steps=7)
    params = [_fd.FDParameter('k', 'k', 0.5)]
    fn = lambda rs, values: rs[0]['[S]']              # noqa: E731
    with monkeypatch.context() as m:
        m.setattr(_fd, '_scan_counters', boom)
        point = _fd.evaluate(private, [cond], params, fn)
        assert point.ok and point.counters == ('n',)
        assert _fd.fd_jacobian(private, [cond], params, fn).counters == ('n',)
    assert local_sensitivity(km, cond, ['k'], outputs=['[S]']).counters == (
        'n',)
    # an explicit list (even empty) still wins
    assert _fd.evaluate(private, [cond], params, fn,
                        event_counters=()).counters == ()


PULSE_TOY = """
model pulse()
  S = 10; P = 0; k = 0.5;
  J: S => P; k*S;
  E1: at time > 2: S = S + 1;
  E2: at time > 4: S = S + 1;
end
"""


def test_a_pulsed_reacting_species_is_not_an_event_counter():
    # S = S + 1 on a species that also reacts: its final value differs
    # between runs whatever the event counts, so it must not be watched.
    km = nsk.KineticModel.from_antimony(PULSE_TOY, units=UNITS)
    assert _fd.detect_event_counters(km) == ()
    ls = local_sensitivity(km, nsk.Condition(t_end=6, n_steps=61), ['k'],
                           outputs=['[S]'], tolerances=TIGHT)
    assert ls.counters == ()
    flags = ls.flags.loc['k']
    assert flags['reliable'] and flags['n_refinements'] == 0
    # S(6) = 10 e^(-6k) + e^(-4k) + e^(-2k) at k = 0.5
    expected = -60 * math.exp(-3) - 4 * math.exp(-2) - 2 * math.exp(-1)
    assert np.isclose(ls.raw[-1, 0, 0], expected, rtol=1e-5)
    # rule variables are not counters either; a boundary species is (read
    # as its concentration, in the constant default compartment)
    km = nsk.KineticModel.from_antimony("""
    model r()
      S = 10; x = 0; c = 0; species $B; B = 0;
      J: S => ; S;
      x' = 0.1;
      z := c + 1;
      E1: at time > 1: x = x + 1, c = c + 1, B = B + 1;
    end
    """, units=UNITS)
    assert _fd.detect_event_counters(km) == ('c', '[B]')


IA_TOY = """
model ia()
  S_init = 10; S = S_init; P = 0; k = 0.5;
  J: S => P; k*S;
end
"""


def test_initial_assignment_inputs_warn_and_init_is_the_remedy():
    # S = S_init is an SBML initial assignment: a perturbed S_init is set
    # after the reset computed S, so dP/dS_init comes out 0 (with a warning);
    # 'init(S_init)' re-evaluates the assignment: dP(6)/dS_init = 1 - e^-3.
    km = nsk.KineticModel.from_antimony(IA_TOY, units=UNITS)
    assert _fd.initial_assignment_inputs(km) == {'S_init': ('S',)}
    cond = nsk.Condition(t_end=6, n_steps=61)
    with pytest.warns(UserWarning, match=r"init\(S_init\)") as record:
        ls = local_sensitivity(km, cond, ['S_init', 'k'], outputs=['[P]'],
                               tolerances=TIGHT)
    assert record[0].filename == __file__        # points at the caller
    assert ls.raw[-1, 0, 0] == 0.0
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        fixed = local_sensitivity(km, cond, ['init(S_init)'], outputs=['[P]'],
                                  tolerances=TIGHT)
    assert np.isclose(fixed.raw[-1, 0, 0], 1 - math.exp(-3), rtol=1e-5)


def test_fd_jacobian_over_two_conditions_with_a_local_parameter():
    # Two experiments differing in S0; a global k (log space) and a local
    # initial value '[S]' of the second experiment only (linear space).
    # f = [S_0(T), S_1(T)] with S_i(T) = S0_i e^(-kT).
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    private = _fd.private_model(km, tolerances=TIGHT, rel_step=1e-3)
    T, k = 2.0, 0.3
    conds = [nsk.Condition(initial={'[S]': 10.0}, t_end=T, n_steps=3),
             nsk.Condition(initial={'[S]': 4.0}, t_end=T, n_steps=3,
                           label='second')]
    params = [_fd.FDParameter('k', 'k', k, log=True, step=1e-3),
              _fd.FDParameter('[S]@second', '[S]', 6.0, log=False, step=1e-3,
                              conditions=(1,))]
    fn = lambda results, values: np.array(  # noqa: E731
        [r['[S]'][-1] for r in results])
    jac = _fd.fd_jacobian(private, conds, params, fn)
    s0 = 10 * math.exp(-k * T)
    s1 = 6.0 * math.exp(-k * T)            # pinned at the parameter's value
    assert np.allclose(jac.f0, [s0, s1], rtol=1e-8)
    expected = np.array([[-k * T * s0, 0.0],
                         [-k * T * s1, math.exp(-k * T)]])
    assert np.allclose(jac.jacobian, expected, rtol=1e-5, atol=1e-12)
    natural = jac.natural
    assert np.allclose(natural[:, 0], [-T * s0, -T * s1], rtol=1e-5)
    assert jac.jacobian[0, 1] == 0.0       # the local one leaves cond 0 alone
    assert jac.conditions[0].parameters == {'k': k}
    assert jac.conditions[1].initial == {'[S]': 6.0}
    assert jac.conditions[1].label == 'second'
    assert list(jac.flags.index) == ['k', '[S]@second']
    assert jac.flags['reliable'].all()


def test_fd_jacobian_virtual_parameter_needs_no_run():
    # f = scale * S(T): an observable scale (target None) is differenced by
    # re-evaluating fn on the nominal results; k needs runs.
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    private = _fd.private_model(km, tolerances=TIGHT)
    T, k, scale = 2.0, 0.5, 3.0
    cond = nsk.Condition(t_end=T, n_steps=3)
    calls = []

    def fn(results, values):
        calls.append(dict(values))
        return np.array([values['scale'] * results[0]['[S]'][-1]])

    params = [_fd.FDParameter('k', 'k', k),
              _fd.FDParameter('scale', None, scale, log=False, step=1e-4)]
    jac = _fd.fd_jacobian(private, [cond], params, fn)
    s = 10 * math.exp(-k * T)
    assert np.allclose(jac.jacobian, [[-k * T * scale * s, s]], rtol=1e-6)
    assert jac.conditions[0].parameters == {'k': k}      # scale not pinned
    assert calls[0] == {'k': k, 'scale': scale}          # nominal point
    assert {'k': k, 'scale': scale + 1e-4} in calls
    assert jac.flags.loc['scale', 'reliable']


def test_fd_jacobian_records_a_failed_perturbation():
    # dS/dt = k S^2, S(0) = 1 blows up at t = 1/k: at k = 0.5 and T = 1.99
    # the k+ run (k = 0.505, blow-up at 1.98) fails while the nominal and
    # k- runs succeed.
    km = nsk.KineticModel.from_antimony("""
    model blow()
      species S; S = 1; k = 0.5;
      J: => S; k*S^2;
    end
    """, units=UNITS)
    private = _fd.private_model(km)
    cond = nsk.Condition(t_end=1.99, n_steps=3)
    params = [_fd.FDParameter('k', 'k', 0.5, step=0.01)]
    fn = lambda rs, values: rs[0]['[S]'][-1:]          # noqa: E731
    # Forward differences have no other side: a NaN column with the error.
    jac = _fd.fd_jacobian(private, [cond], params, fn, method='forward')
    flags = jac.flags.loc['k']
    assert not flags['reliable'] and flags['error'].startswith('condition 0')
    assert np.isnan(jac.jacobian).all()
    # Central differences fall back to the side that ran (k-), reliably.
    jac = _fd.fd_jacobian(private, [cond], params, fn)
    flags = jac.flags.loc['k']
    assert flags['reliable'] and flags['error'] is None
    assert flags['scheme'] == 'backward' and flags['step'] == 0.01
    s0 = 1 / (1 - 0.5 * 1.99)
    s_minus = 1 / (1 - 0.5 * math.exp(-0.01) * 1.99)
    assert np.isclose(jac.jacobian[0, 0], (s0 - s_minus) / 0.01, rtol=1e-4)


def test_one_failed_side_of_a_central_difference_falls_back_to_the_other():
    # [E] = 0 is differenced linearly; its central minus run (E = -0.001)
    # makes sqrt(E + P) undefined and fails. The plus side is fine.
    km = nsk.KineticModel.from_antimony("""
    model sq()
      S = 4; P = 0; E = 0; k = 0.5;
      J: S => P; k*S;
      J2: => E; k*sqrt(E + P);
    end
    """, units=UNITS)
    cond = nsk.Condition(t_end=2, n_steps=5)
    with pytest.warns(UserWarning) as record:
        ls = local_sensitivity(km, cond, ['[E]', 'k'], outputs=['[E]'])
    messages = [str(w.message) for w in record]
    assert any("central difference of ['[E]'] failed" in m for m in messages)
    assert any("['[E]'] are at 0" in m for m in messages)
    flags = ls.flags.loc['[E]']
    assert flags['scheme'] == 'forward' and flags['reliable']
    assert flags['error'] is None and not flags['event_count_changed']
    assert ls.flags.loc['k', 'scheme'] == 'central'
    with pytest.warns(UserWarning, match=r"\['\[E\]'\] are at 0"):
        forward = local_sensitivity(km, cond, ['[E]', 'k'], outputs=['[E]'],
                                    method='forward')
    assert np.allclose(ls.raw[:, 0, 0], forward.raw[:, 0, 0], rtol=1e-12)
    assert np.isfinite(ls.raw).all()


def test_a_failed_side_with_a_count_change_on_the_other_is_refined():
    # q = 0 is differenced linearly: q - 0.001 makes sqrt(q) undefined and
    # fails, and q + 0.001 moves the firing to t = 5.5, before t_end; at
    # step 1e-4 it fires at 6.4, after t_end. The column was NaN, unrefined,
    # with event_count_changed=False.
    km = nsk.KineticModel.from_antimony("""
    model ev()
      S = 10; n = 0; q = 0; k = 0.5;
      J: S => ; k*S + sqrt(q)^2*S;
      E1: at time > 6.5 - 1000*q: n = n + 1;
    end
    """, units=UNITS)
    cond = nsk.Condition(t_end=6, n_steps=61)
    with pytest.warns(UserWarning) as record:
        ls = local_sensitivity(km, cond, ['q'], outputs=['[S]'])
    messages = [str(w.message) for w in record]
    assert any("central difference of ['q'] failed" in m for m in messages)
    assert any("['q'] are at 0" in m for m in messages)
    assert ls.counters == ('n',)
    flags = ls.flags.loc['q']
    assert flags['n_refinements'] == 1 and flags['scheme'] == 'forward'
    assert flags['step'] == pytest.approx(1e-4, rel=1e-12)
    assert flags['reliable'] and not flags['event_count_changed']
    assert flags['error'] is None
    # dS/dq = -t S (the rate is (k + q) S for q >= 0); forward error O(h t)
    assert ls.raw[-1, 0, 0] == pytest.approx(-6 * 10 * math.exp(-3),
                                             rel=1e-3)
    # Without refinement the running side never agrees: flagged, NaN, and
    # the failed side's error.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        flat = local_sensitivity(km, cond, ['q'], outputs=['[S]'],
                                 max_refinements=0)
    flags = flat.flags.loc['q']
    assert flags['event_count_changed'] and not flags['reliable']
    assert flags['n_refinements'] == 0 and 'condition 0' in flags['error']
    assert np.isnan(flat.raw).all()


def test_a_small_abs_step_warns_about_integration_noise():
    # abs_step sets a linear step that private_model's rel_step check never
    # saw: 1e-9 on k1 = 0.5 is a relative step of 2e-9 < 100 rtol (1e-8).
    km = nsk.KineticModel.from_antimony("""
    model s()
      S = 10; P = 0; k1 = 0.5; k2 = 1;
      J: S => P; (k1 + k2)*S;
    end
    """, units=UNITS)
    cond = nsk.Condition(t_end=6, n_steps=13)
    with pytest.warns(UserWarning, match=r"\['k1'\] use a step \(refined, "
                                         r'or set by abs_step\)'):
        local_sensitivity(km, cond, ['k1', 'k2'], outputs=['[S]'],
                          abs_step={'k1': 1e-9})
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        ls = local_sensitivity(km, cond, ['k1', 'k2'], outputs=['[S]'],
                               abs_step={'k1': 1e-6})
    t = ls.time        # FD noise ~ 3 rtol / 2e-6 = 1.5e-4 relative
    assert np.allclose(ls.raw[:, 0, 0], -t * 10 * np.exp(-1.5 * t),
                       rtol=1e-3, atol=1e-9)


def test_fd_jacobian_records_fn_failing_on_a_perturbed_run():
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    private = _fd.private_model(km)
    cond = nsk.Condition(t_end=2, n_steps=3)

    def fn(results, values):
        if values['k'] != 0.5:
            raise RuntimeError('boom')
        return results[0]['[S]'][-1:]

    jac = _fd.fd_jacobian(private, [cond], [_fd.FDParameter('k', 'k', 0.5)],
                          fn)
    flags = jac.flags.loc['k']
    assert np.isnan(jac.jacobian).all() and not flags['reliable']
    assert flags['error'] == 'fn: RuntimeError: boom'

    # fn failing on one side only: the other side's one-sided difference.
    def fn_minus(results, values):
        if values['k'] < 0.5:
            raise RuntimeError('boom')
        return results[0]['[S]'][-1:]

    jac = _fd.fd_jacobian(private, [cond], [_fd.FDParameter('k', 'k', 0.5)],
                          fn_minus)
    flags = jac.flags.loc['k']
    assert flags['reliable'] and flags['scheme'] == 'forward'
    assert np.isclose(jac.jacobian[0, 0], -0.5 * 2 * 10 * math.exp(-1),
                      rtol=1e-2)                  # dS/d ln k = -k t S


def _reset_doubling_g(km, double=False, **kwargs):
    km.reset(mode='baseline')
    if double:
        km._te['g'] = 2


def test_fd_jacobian_starts_from_the_settings_as_built():
    # A reset_func writes g = 2 for one experiment's reset kwargs only. A
    # run of that experiment on the private model (as a calibration residual
    # evaluation would make) must not leak g = 2 into the other experiment.
    km = nsk.KineticModel.from_antimony("""
    model gain()
      S = 10; P = 0; k = 0.5; g = 1;
      J: S => P; g*k*S;
    end
    """, units=UNITS)
    km.reset_func = _reset_doubling_g
    c0 = nsk.Condition(t_end=2, n_steps=3)
    c1 = c0.replace(reset_kwargs={'double': True})
    private = _fd.private_model(km, tolerances=TIGHT)
    assert private._fd_settings == km._user_parameter_values()
    private.run(c1)
    assert private._te['g'] == 2.0
    fn = lambda rs, values: np.array([r['[P]'][-1] for r in rs])  # noqa
    params = [_fd.FDParameter('k', 'k', 0.5)]
    jac = _fd.fd_jacobian(private, [c0, c1], params, fn)
    fresh = _fd.fd_jacobian(_fd.private_model(km, tolerances=TIGHT),
                            [c0, c1], params, fn)
    assert jac.f0[0] == fresh.f0[0]
    assert np.isclose(jac.f0[0], 10 * (1 - math.exp(-1)), rtol=1e-8)
    assert np.isclose(jac.f0[1], 10 * (1 - math.exp(-2)), rtol=1e-8)
    assert np.array_equal(jac.jacobian, fresh.jacobian)
    assert jac.tolerances == TIGHT
    # an explicit snapshot wins
    other = _fd.fd_jacobian(private, [c0], params, fn,
                            settings={**private._fd_settings, 'g': 3.0})
    assert np.isclose(other.f0[0], 10 * (1 - math.exp(-3)), rtol=1e-8)


def test_a_model_built_from_a_blueprint_records_its_fd_settings(monkeypatch):
    # copy() and blueprint builds (a process-pool worker's model) do not
    # carry private_model's records: the first call's snapshot must be kept,
    # or g = 2 written by c1's reset_func leaked into the next c0 run.
    km = nsk.KineticModel.from_antimony("""
    model gain()
      S = 10; P = 0; k = 0.5; g = 1;
      J: S => P; g*k*S;
    end
    """, units=UNITS)
    km.reset_func = _reset_doubling_g
    c0 = nsk.Condition(t_end=2, n_steps=21)
    c1 = c0.replace(reset_kwargs={'double': True})
    private = _fd.private_model(km)
    # built as _map_with_model builds a worker's model (settings as built)
    worker, other = (private._blueprint().build() for _ in range(2))
    params = [_fd.FDParameter('k', 'k', 0.5)]
    fn = lambda rs, values: np.array([r['[P]'][-1] for r in rs])  # noqa
    want = _fd.evaluate(private, [c0, c1], params, fn).f
    want_jac = _fd.fd_jacobian(private, [c0, c1], params, fn).jacobian
    assert np.allclose(want, [10 * (1 - math.exp(-1)),
                              10 * (1 - math.exp(-2))], rtol=1e-8)
    assert not hasattr(worker, '_fd_settings')
    assert not hasattr(worker, '_fd_counters')
    first = _fd.evaluate(worker, [c0, c1], params, fn).f
    assert worker._fd_settings == private._fd_settings
    assert worker._fd_counters == () and worker._fd_unwatched == ()

    def boom(*args, **kwargs):
        raise AssertionError('event counters re-detected')

    monkeypatch.setattr(_fd, '_scan_counters', boom)
    second = _fd.evaluate(worker, [c0, c1], params, fn).f
    jac = _fd.fd_jacobian(worker, [c0, c1], params, fn)
    assert np.array_equal(first, want) and np.array_equal(second, want)
    assert np.allclose(jac.jacobian, want_jac, rtol=1e-12)
    monkeypatch.undo()
    # An invalid call, or explicit settings/counters, records nothing.
    with pytest.raises(ValueError, match='Duplicate'):
        _fd.evaluate(other, [c0], params * 2, fn)
    _fd.evaluate(other, [c0], params, fn, settings=private._fd_settings,
                 event_counters=())
    assert not hasattr(other, '_fd_settings')
    assert not hasattr(other, '_fd_counters')


def test_fd_helpers_for_linear_steps_and_fixed_step_integrators():
    # a scalar abs_step applies to a parameter at 0 (linear space)
    p, = _fd.make_parameters({'z': 0.0}, rel_step=1e-3, abs_step=0.01)
    assert not p.log and p.step == 0.01
    assert _fd.relative_step(2.0, False, 0.01) == 0.005
    assert math.isnan(_fd.relative_step(0.0, False, 0.01))
    # tolerances= needs an integrator that has them
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    km._te.setIntegrator('rk4')
    with pytest.raises(ValueError, match='needs an integrator with '
                                         'tolerances'):
        local_sensitivity(km, nsk.Condition(t_end=1, n_steps=3), ['k'],
                          tolerances=TIGHT)
    assert _fd.tolerances_of(_fd.private_model(km)) is None


def test_fd_jacobian_validation():
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    cond = nsk.Condition(t_end=1, n_steps=3)
    fn = lambda rs, values: rs[0]['[S]']             # noqa: E731
    p = _fd.FDParameter('k', 'k', 0.5)
    with pytest.raises(ValueError, match='applies to condition'):
        _fd.fd_jacobian(km, [cond], [_fd.FDParameter('k', 'k', 0.5,
                                                     conditions=(1,))], fn)
    with pytest.raises(ValueError, match='several parameters'):
        _fd.fd_jacobian(km, [cond], [p, _fd.FDParameter('k2', 'k', 0.5)], fn)
    with pytest.raises(ValueError, match='method'):
        _fd.fd_jacobian(km, [cond], [p], fn, method='complex')
    with pytest.raises(ValueError, match='log space'):
        _fd.FDParameter('k', 'k', 0.0)
    with pytest.raises(ValueError):
        _fd.FDParameter('k', 'k', 1.0, step=0.0)
    with pytest.raises(TypeError, match='max_refinements'):
        _fd.fd_jacobian(km, [cond], [p], fn, max_refinements=1.5)
    with pytest.raises(TypeError, match='settings'):
        _fd.fd_jacobian(km, [cond], [p], fn, settings=[('k', 1.0)])
    # The shared argument checks F3 relies on, before any run.
    with pytest.raises(ValueError, match='At least one condition'):
        _fd.fd_jacobian(km, [], [p], fn)
    with pytest.raises(TypeError, match='FDParameter objects'):
        _fd.fd_jacobian(km, [cond], ['k'], fn)
    with pytest.raises(ValueError, match='Duplicate parameter names'):
        _fd.fd_jacobian(km, [cond], [p, _fd.FDParameter('k', '[S]', 10.0)],
                        fn)
    with pytest.raises(ValueError, match='at least one parameter'):
        _fd.fd_jacobian(km, [cond], [], fn)
    with pytest.raises(TypeError, match='callable'):
        _fd.fd_jacobian(km, [cond], [p], None)
    # fn must return a 1-D array (checked on the nominal evaluation).
    with pytest.raises(ValueError, match='1-D'):
        _fd.fd_jacobian(km, [cond], [p],
                        lambda rs, values: np.ones((2, 2)))
    # A bare Condition is one condition.
    private = _fd.private_model(km)
    jac = _fd.fd_jacobian(private, cond, [p], fn)
    assert jac.jacobian.shape == (3, 1) and len(jac.conditions) == 1
    # An output length that changes away from the nominal point only spoils
    # that column (both sides, so no one-sided fallback either).
    jac = _fd.fd_jacobian(private, [cond], [p],
                          lambda rs, v: np.ones(1 if v['k'] == 0.5 else 2))
    flags = jac.flags.loc['k']
    assert np.isnan(jac.jacobian).all() and not flags['reliable']
    assert 'returned 2 values' in flags['error']


def test_evaluate_is_the_nominal_point_and_fd_jacobian_reuses_it(
        monkeypatch):
    # The calibration pattern: evaluate (the residual function) at a point,
    # then the Jacobian there without repeating the nominal runs.
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    private = _fd.private_model(km, tolerances=TIGHT)
    conds = [nsk.Condition(initial={'[S]': 10.0}, t_end=2, n_steps=3),
             nsk.Condition(initial={'[S]': 4.0}, t_end=2, n_steps=3)]
    params = [_fd.FDParameter('k', 'k', 0.3),
              _fd.FDParameter('[S]@1', '[S]', 6.0, log=False, step=1e-3,
                              conditions=(1,))]
    fn = lambda rs, values: np.array(  # noqa: E731
        [r['[S]'][-1] for r in rs])
    runs = []
    original = nsk.KineticModel.run

    def counting_run(self, *args, **kwargs):
        runs.append(1)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(nsk.KineticModel, 'run', counting_run)
    point = _fd.evaluate(private, conds, params, fn)
    assert point.ok and point.error is None and len(runs) == 2
    assert point.status == ['ok', 'ok'] and point.counters == ()
    assert point.values == {'k': 0.3, '[S]@1': 6.0}
    assert point.conditions[1].initial == {'[S]': 6.0}
    del runs[:]
    full = _fd.fd_jacobian(private, conds, params, fn)
    n_full = len(runs)
    del runs[:]
    reused = _fd.fd_jacobian(private, conds, params, fn, nominal=point)
    assert len(runs) == n_full - len(conds)
    assert np.array_equal(point.f, full.f0)
    assert np.array_equal(reused.f0, full.f0)
    assert np.array_equal(reused.jacobian, full.jacobian)
    pd.testing.assert_frame_equal(reused.flags, full.flags)
    assert reused.conditions == full.conditions
    # a point of another problem, or a failed one, is refused before runs
    del runs[:]
    other = [_fd.FDParameter('k', 'k', 0.31), params[1]]
    with pytest.raises(ValueError, match="different \\['parameters', "
                                         "'conditions'\\]"):
        _fd.fd_jacobian(private, conds, other, fn, nominal=point)
    with pytest.raises(ValueError, match="different \\['conditions'\\]"):
        _fd.fd_jacobian(private, [conds[0], conds[1].replace(t_end=3)],
                        params, fn, nominal=point)
    with pytest.raises(TypeError, match='FDPoint'):
        _fd.fd_jacobian(private, conds, params, fn, nominal=point.f)
    assert runs == []


def test_evaluate_records_a_failed_run_without_raising():
    # dS/dt = k S^2 blows up at t = 1/k: condition 1 (k = 0.6, T = 1.99)
    # fails, condition 0 (T = 1) does not; the point records it.
    km = nsk.KineticModel.from_antimony("""
    model blow()
      species S; S = 1; k = 0.6;
      J: => S; k*S^2;
    end
    """, units=UNITS)
    private = _fd.private_model(km)
    conds = [nsk.Condition(t_end=1, n_steps=3),
             nsk.Condition(t_end=1.99, n_steps=3, label='late')]
    calls = []

    def fn(rs, values):
        calls.append(1)
        return np.array([r['[S]'][-1] for r in rs])

    params = [_fd.FDParameter('k', 'k', 0.6)]
    point = _fd.evaluate(private, conds, params, fn)
    assert not point.ok and point.f is None and calls == []
    assert point.status == ['ok', 'failed']
    assert point.results[0] is not None and point.results[1] is None
    assert point.error.startswith('condition 1: ')
    with pytest.raises(ValueError, match='failed evaluation'):
        _fd.fd_jacobian(private, conds, params, fn, nominal=point)
    with pytest.raises(KineticSimulationError,
                       match="condition 1 \\('late'\\) failed"):
        _fd.fd_jacobian(private, conds, params, fn)


# --- shipped model (slow) -------------------------------------------------------
@pytest.mark.slow
def test_shipped_model_sensitivities_equal_direct_differences():
    # The shipped fed-batch model: a custom reset_func with reset kwargs, a
    # FeedSpike counter, and k_13 = 0 in its baseline (linear differencing).
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import te_r
    km = te_r.copy()
    # te_r is shared by the test session and other tests leave parameter
    # edits on it; start this private copy from its as-loaded values so the
    # test does not depend on test order.
    for p in km._user_parameter_values():
        km._te[p] = km._te[f'init({p})']
    before = _model_state(km)
    cond = nsk.Condition(t_end=60, n_steps=121,
                         reset_kwargs={'reset_spike_cap': True})
    params = ['k_7ie', 'k_7', 'k_10', 'k_13']
    # k_13 = 0 has no abs_step: its default step (rel_step in its own
    # units) is reported.
    with pytest.warns(UserWarning, match=r"\['k_13'\] are at 0"):
        ls = local_sensitivity(km, cond, params,
                               outputs=['[s_glu]', '[s_EtOH]'],
                               metrics={'EtOH': Final('[s_EtOH]'),
                                        'n_spk': EventCount('n_glu_spikes')})
    assert _model_state(km) == before
    assert ls.counters == ('n_glu_spikes',)
    flags = ls.flags
    assert flags['reliable'].all() and not flags['event_count_changed'].any()
    assert flags['log'].tolist() == [True, True, True, False]
    assert flags.loc['k_13', 'value'] == 0.0
    assert flags.loc['k_13', 'step'] == 1e-3
    # 4 spikes at every perturbation: the count's derivative is exactly 0
    n_spk = ls.metric_sensitivities.query("metric == 'n_spk'")
    assert (n_spk['metric_value'] == 4).all()
    assert (n_spk['derivative'] == 0).all() and n_spk['reliable'].all()
    # Each column equals a central difference of two direct runs on a copy
    # of the model with the tolerances the analysis used (it runs on such a
    # copy, so they agree to round-off), and the ethanol metric equals the
    # trajectory's last row.
    assert ls.tolerances == (1e-12, 1e-10)
    copy = km.copy()
    copy._te.getIntegrator().setValue('absolute_tolerance', ls.tolerances[0])
    copy._te.getIntegrator().setValue('relative_tolerance', ls.tolerances[1])
    for j, p in enumerate(params):
        v = ls.values[j]
        if flags.loc[p, 'log']:     # d/dp = d/d(ln p) / p, h = 1e-3 in ln p
            lo, hi, width = v * math.exp(-1e-3), v * math.exp(1e-3), 2e-3 * v
        else:                       # p = 0: h = rel_step in model units
            lo, hi, width = v - 1e-3, v + 1e-3, 2e-3
        up, down = (copy.run(cond.with_parameters({p: x})).final['[s_EtOH]']
                    for x in (hi, lo))
        assert np.isclose(ls.raw[-1, 1, j], (up - down) / width,
                          rtol=1e-9, atol=0)
        assert ls.metric_derivatives[0, j] == ls.raw[-1, 1, j]
    # every parameter moves the final ethanol titer across the 4 spikes
    assert (np.abs(ls.raw[-1, 1]) > 0).all()
    top = ls.rank('EtOH').iloc[0]
    assert top['parameter'] == 'k_7ie' and top['elasticity'] < -0.5
