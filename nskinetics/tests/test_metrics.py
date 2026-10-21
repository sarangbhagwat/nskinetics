# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Kinetics-level metrics (F2 item 1; spec test d), each checked against a
closed form on an inline toy."""
import math
import pickle

import numpy as np
import pandas as pd
import pytest
import tellurium as te

import nskinetics as nsk
from nskinetics.analysis import (At, CallableMetric, EventCount, Final,
                                 Integral, Max, Metric, Min, Productivity,
                                 TimeAtMax, Yield)

UNITS = {'time': 'h', 'conc': 'g/L'}

FEEDSPIKE_TOY = """
model fb()
  compartment env = 1; species S in env, P in env;
  S = 50; P = 0; k = 0.5; threshold = 5; target = 50; feed_conc = 500;
  n_max = 3; n_spk = 0; last_vol = 0; tot_vol = 0;
  J: S => P; k*S*env;
end
"""

# S => P with S(0) = 10: S = 10 e^(-kt), P = 10 (1 - e^(-kt)).
DECAY_TOY = """
model decay()
  S = 10; P = 0; k = 0.5;
  J: S => P; k*S;
end
"""

# A -> B -> C: B peaks at t* = ln(k1/k2) / (k1 - k2).
CHAIN_TOY = """
model chain()
  A = 1; B = 0; C = 0; k1 = 1; k2 = 0.25;
  J1: A => B; k1*A;
  J2: B => C; k2*B;
end
"""

# Every time S falls to 5 it is topped up by 5: the n-th reset is at
# t = n ln2 / k, so by t = 6 at k = 0.5 there are floor(6k / ln2) = 4.
EVENT_TOY = """
model ev()
  S = 10; n = 0; k = 0.5;
  J: S => ; k*S;
  E1: at S < 5: S = S + 5, n = n + 1;
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


@pytest.fixture(scope='module')
def decay():
    """One decay run over [0, 6] with 601 points, k = 0.5."""
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    km._te.getIntegrator().setValue('relative_tolerance', 1e-10)
    km._te.getIntegrator().setValue('absolute_tolerance', 1e-12)
    return km.run(nsk.Condition(t_end=6, n_steps=601))


@pytest.fixture(scope='module')
def feedspike():
    km = _feedspike_model()
    cond = nsk.Condition(initial={'[S]': 50}, t_end=48, n_steps=481,
                         selections=['time', '[S]', '[P]', 'env', 'tot_vol',
                                     'feed_conc', 'n_spk'])
    return km, cond, km.run(cond)


# --- (d) metrics on the FeedSpike toy -------------------------------------------
def test_yield_added_matches_closed_form_from_result_columns(feedspike):
    _, _, res = feedspike
    f = res.final
    expected = (f['[P]'] * f['env']
                / (50.0 * 1.0 + f['tot_vol'] * f['feed_conc']))
    y = Yield('[P]', '[S]', volume='env', fed=('tot_vol', 'feed_conc'),
              basis='added')
    assert y(res) == pytest.approx(expected, rel=1e-12)
    # feed_conc as a number gives the same value
    y_num = Yield('[P]', '[S]', volume='env', fed=('tot_vol', 500.0),
                  basis='added')
    assert y_num(res) == pytest.approx(expected, rel=1e-12)
    assert f['n_spk'] == 3 and f['tot_vol'] > 0     # the feed was used


def test_yield_consumed_is_one_by_mass_balance(feedspike):
    # S => P conserves amount and each spike adds exactly
    # last_vol * feed_conc of S, so P formed / S consumed is exactly 1.
    _, _, res = feedspike
    y = Yield('[P]', '[S]', volume='env', fed=('tot_vol', 'feed_conc'))
    assert y.basis == 'consumed'
    assert y(res) == pytest.approx(1.0, rel=1e-8)


def test_productivity_with_turnaround(feedspike):
    _, _, res = feedspike
    p_end = res.final['[P]']
    assert Productivity('[P]')(res) == pytest.approx(p_end / 48.0, rel=1e-12)
    assert Productivity('[P]', turnaround=12.0)(res) == pytest.approx(
        p_end / 60.0, rel=1e-12)


def test_zero_denominator_is_nan_under_numpy_raise():
    # importing nskinetics (biosteam) sets divide/invalid to 'raise'; the
    # metric must still return NaN, not raise or give an infinity.
    df = pd.DataFrame({'time': [0.0, 1.0], '[P]': [0.0, 1.0],
                       '[S]': [0.0, 0.0], 'V': [1.0, 1.0]})
    res = nsk.SimulationResult(df)
    old = np.seterr(all='raise')
    try:
        assert math.isnan(Yield('[P]', '[S]', basis='added')(res))   # 1 / 0
        assert math.isnan(Yield('[P]', '[S]')(res))                  # 1 / 0
        flat = nsk.SimulationResult(df.assign(**{'[P]': [0.0, 0.0]}))
        assert math.isnan(Yield('[P]', '[S]', volume='V')(flat))     # 0 / 0
        # a zero-length run plus zero turnaround
        one_row = nsk.SimulationResult(df.iloc[:1].copy())
        assert math.isnan(Productivity('[P]')(one_row))
    finally:
        np.seterr(**old)


def test_metrics_plug_into_run_many_via_name_and_selections():
    km = _feedspike_model()
    # The condition records only time and [P]; the metrics' own selections
    # (env, tot_vol, feed_conc, [S], n_spk) are recorded by run_many.
    conds = [nsk.Condition(t_end=48, n_steps=97, selections=['time', '[P]'],
                           parameters={'n_max': n}) for n in (1, 3)]
    metrics = [Final('[P]'), EventCount('n_spk'),
               Yield('[P]', '[S]', volume='env', fed=('tot_vol', 'feed_conc'))]
    batch = nsk.run_many(km, conds, metrics=metrics)
    m = batch.metrics
    assert list(m.columns[3:]) == ['Final([P])', 'EventCount(n_spk)',
                                   'Yield([P]/[S], consumed)']
    assert list(m['status']) == ['ok', 'ok']
    assert list(m['EventCount(n_spk)']) == [1.0, 3.0]
    assert np.allclose(m['Yield([P]/[S], consumed)'], 1.0, rtol=1e-8)
    assert 'n_spk' in batch.results[0] and 'tot_vol' in batch.results[0]


# --- each metric against a closed form ------------------------------------------
def test_final_max_min_on_decay(decay):
    k, T = 0.5, 6.0
    s_end = 10 * math.exp(-k * T)
    assert Final('[S]')(decay) == pytest.approx(s_end, rel=1e-8)
    assert Max('[P]')(decay) == pytest.approx(10 - s_end, rel=1e-8)
    assert Min('[S]')(decay) == pytest.approx(s_end, rel=1e-8)
    assert Max('[S]')(decay) == 10.0
    assert Min('[P]')(decay) == 0.0


def test_at_on_and_between_grid_points(decay):
    s = lambda t: 10 * math.exp(-0.5 * t)                      # noqa: E731
    assert At('[S]', 2.0)(decay) == pytest.approx(s(2.0), rel=1e-8)
    # between grid points (dt = 0.01): linear interpolation of the exact
    # values at the neighbouring points
    t = 2.004
    expected = s(2.0) + (s(2.01) - s(2.0)) * 0.4
    assert At('[S]', t)(decay) == pytest.approx(expected, rel=1e-8)
    with pytest.raises(ValueError, match='outside'):
        At('[S]', 7.0)(decay)


def test_at_checks_the_condition_time_span_before_running():
    # the hook the analyses call before any run
    at = At('[S]', 7.0)
    at._check_condition(nsk.Condition(t_end=7, n_steps=3))      # on the edge
    at._check_condition(nsk.Condition(times=[6, 7, 8]))
    for cond in (nsk.Condition(t_end=6, n_steps=3),
                 nsk.Condition(times=[7.5, 8])):
        with pytest.raises(ValueError, match='outside the simulated range'):
            at._check_condition(cond)
    Final('[S]')._check_condition(nsk.Condition(t_end=1, n_steps=3))


def test_integral_is_trapezoid_of_exponential(decay):
    k, T = 0.5, 6.0
    exact = 10 / k * (1 - math.exp(-k * T))
    # trapezoid error on dt = 0.01 is ~ dt^2/12 * f'(t) integrated: < 1e-5
    assert Integral('[S]')(decay) == pytest.approx(exact, rel=1e-5)


def test_time_at_max_of_chain_intermediate():
    km = nsk.KineticModel.from_antimony(CHAIN_TOY, units=UNITS)
    res = km.run(nsk.Condition(t_end=10, n_steps=1001))
    t_star = math.log(1 / 0.25) / (1 - 0.25)
    assert abs(TimeAtMax('[B]')(res) - t_star) <= 0.01 / 2 + 1e-12
    # the parabolic refinement is far closer than the grid (O(dt^2))
    assert abs(TimeAtMax('[B]')(res) - t_star) < 1e-4


def test_time_at_max_refines_the_grid_maximum_by_a_parabola():
    # exact for a parabola, also on an uneven grid
    t = np.array([0.0, 0.5, 1.0, 1.7, 2.0, 3.0])
    df = pd.DataFrame({'time': t, 'y': 5 - (t - 1.234) ** 2})
    assert TimeAtMax('y')(df) == pytest.approx(1.234, rel=1e-12)
    # a maximum on the first or last row, or tied with the next row
    # (a plateau), is the grid time itself
    assert TimeAtMax('y')(pd.DataFrame({'time': t, 'y': -t})) == 0.0
    assert TimeAtMax('y')(pd.DataFrame({'time': t, 'y': t})) == 3.0
    flat = np.array([0.0, 1.0, 2.0, 2.0, 2.0, 1.0])
    assert TimeAtMax('y')(pd.DataFrame({'time': t, 'y': flat})) == 1.0
    nan = pd.DataFrame({'time': t, 'y': [0, 1, np.nan, 1, 0, 0]})
    assert math.isnan(TimeAtMax('y')(nan))


def test_productivity_and_yield_without_volume(decay):
    k, T = 0.5, 6.0
    p_end = 10 * (1 - math.exp(-k * T))
    assert Productivity('[P]')(decay) == pytest.approx(p_end / T, rel=1e-8)
    assert Productivity('[P]', turnaround=2.0)(decay) == pytest.approx(
        p_end / (T + 2.0), rel=1e-8)
    assert Yield('[P]', '[S]')(decay) == pytest.approx(1.0, rel=1e-8)
    assert Yield('[P]', '[S]', basis='added')(decay) == pytest.approx(
        1 - math.exp(-k * T), rel=1e-8)


def test_event_count_is_number_of_firings():
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    cond = nsk.Condition(t_end=6, n_steps=61, selections=['time', '[S]', 'n'])
    expected = math.floor(6 * 0.5 / math.log(2))
    assert expected == 4
    assert EventCount('n')(km.run(cond)) == expected
    # the basis is final minus initial, whatever the counter starts at
    assert EventCount('n')(km.run(cond.with_initial(n=2))) == expected


def test_callable_metric_and_from_function(decay):
    fn = lambda res: res.final['[P]'] / res.final['[S]']        # noqa: E731
    m = CallableMetric(fn, 'ratio', selections=['[P]', '[S]'])
    k, T = 0.5, 6.0
    assert m(decay) == pytest.approx(math.exp(k * T) - 1, rel=1e-8)
    assert m.selections == ('[P]', '[S]') and m.name == 'ratio'
    m2 = Metric.from_function(fn, 'ratio', ['[P]', '[S]'])
    assert isinstance(m2, CallableMetric) and m2 == m
    with pytest.raises(TypeError, match='real scalar'):
        CallableMetric(lambda res: res['[P]'], 'col')(decay)
    with pytest.raises(TypeError, match='callable'):
        CallableMetric(3, 'x')
    # The base class computes nothing; its error says how to wrap a function.
    with pytest.raises(NotImplementedError, match='CallableMetric'):
        Metric('x')(decay)


def test_metrics_accept_a_dataframe(decay):
    assert Final('[S]')(decay.df) == Final('[S]')(decay)
    with pytest.raises(TypeError, match='SimulationResult'):
        Final('[S]')({'time': [0]})


# --- protocol, names, validation ------------------------------------------------
def test_default_names_selections_and_equality():
    y = Yield('[P]', '[S]', volume='env', fed=('tot_vol', 'feed_conc'),
              basis='added')
    assert y.name == 'Yield([P]/[S], added)'
    assert y.selections == ('[P]', '[S]', 'env', 'tot_vol', 'feed_conc')
    assert Yield('[P]', '[S]', volume='env', fed=('tot_vol', 2.0)
                 ).selections == ('[P]', '[S]', 'env', 'tot_vol')
    assert At('[S]', 12).name == 'At([S], 12)'
    assert Final('[P]', name='titer').name == 'titer'
    assert Productivity('[P]', turnaround=12).selections == ('[P]',)
    assert EventCount('n_spk').selections == ('n_spk',)
    assert Final('[P]') == Final('[P]') and Final('[P]') != Max('[P]')
    assert hash(Final('[P]')) == hash(Final('[P]'))
    assert len({Final('[P]'), Final('[P]'), Final('[S]')}) == 2
    namespace = {'Final': Final, 'Yield': Yield, 'At': At,
                 'Productivity': Productivity, 'EventCount': EventCount}
    for m in (Final('[P]'), y, At('[S]', 1.5), Productivity('[P]',
                                                             turnaround=3),
              EventCount('n', name='spikes')):
        assert pickle.loads(pickle.dumps(m)) == m
        assert eval(repr(m), namespace) == m                # noqa: S307
    assert repr(Final('[P]')) == "Final(selection='[P]', name='Final([P])')"


@pytest.mark.parametrize('make, error', [
    (lambda: Final(3), TypeError),
    (lambda: Final(''), ValueError),
    (lambda: Final('time'), ValueError),
    (lambda: At('[S]', float('nan')), ValueError),
    (lambda: At('[S]', '1'), TypeError),
    (lambda: Productivity('[P]', turnaround=-1), ValueError),
    (lambda: Productivity('[P]', turnaround=True), TypeError),
    (lambda: Yield('[P]', '[S]', basis='theoretical'), ValueError),
    (lambda: Yield('[P]', '[S]', fed=('tot_vol', 500)), ValueError),
    (lambda: Yield('[P]', '[S]', volume='env', fed='tot_vol'), TypeError),
    (lambda: Yield('[P]', '[S]', volume='env', fed=('tot_vol',)), TypeError),
    (lambda: Yield('[P]', '[S]', volume='env', fed=('tot_vol', None)),
     TypeError),
])
def test_invalid_arguments_raise(make, error):
    with pytest.raises(error):
        make()


def test_missing_column_raises_key_error(decay):
    with pytest.raises(KeyError, match='not recorded'):
        Final('[X]')(decay)
