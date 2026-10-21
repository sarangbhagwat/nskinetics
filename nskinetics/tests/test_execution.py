# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""run_many / scan, RunBatch / ScanResult (F1 item 11; spec tests i, j, k)."""
import functools
import math
import multiprocessing
import os

import matplotlib
import numpy as np
import pandas as pd
import pytest
import tellurium as te

import nskinetics as nsk
from nskinetics.engine import execution
from nskinetics.engine.kinetic_model import _ModelBlueprint
from nskinetics.exceptions import KineticSimulationError

matplotlib.use('Agg')

FEEDSPIKE_TOY = """
model fb()
  compartment env = 1; species S in env, P in env;
  S = 50; P = 0; k = 0.5; threshold = 5; target = 50; feed_conc = 500;
  n_max = 3; n_spk = 0; last_vol = 0; tot_vol = 0;
  J: S => P; k*S*env;
end
"""

# dS/dt = k*S^2 with S(0) = 1: S = 1/(1 - k*t), which blows up at t = 1/k.
BLOW_UP_TOY = """
model blow()
  species S; S = 1; k = 1;
  J: => S; k*S^2;
end
"""

UNITS = {'time': 'h', 'conc': 'g/L'}


def _feedspike_model(reset=None):
    km = nsk.KineticModel(te.loadAntimonyModel(FEEDSPIKE_TOY), units=UNITS,
                          reset=reset)
    spike = nsk.FeedSpike('S', 'S <= threshold', 'target', 'feed_conc',
                          volume_var='env', max_count='n_max',
                          count_var='n_spk', last_vol_var='last_vol',
                          tot_vol_var='tot_vol')
    for event in spike.expand():
        km.add_event(event)
    km.compile_events()
    return km


def _p_end(res):
    return res.final['[P]']


def _model_state(km):
    """Everything A7 says a batch must leave unchanged."""
    return (km._user_parameter_values(),
            {v: km._te[v] for v in km.baseline_variables},
            list(km.events), km._events_compiled, km.results_array,
            km._te.getCurrentSBML())


class SpikeCount:
    """A metric object in the analysis Metric protocol (name, selections)."""
    name = 'n_spikes'
    selections = ('n_spk',)

    def __call__(self, res):
        return res['n_spk'][-1]


def _reset_with_lambda_attribute(model, **kwargs):
    model._te.reset()


@pytest.fixture
def close_figures():
    """Close every figure after the test, also when an assertion fails."""
    import matplotlib.pyplot as plt
    try:
        yield
    finally:
        plt.close('all')


# --- (i) failures are recorded --------------------------------------------------
def test_run_many_records_a_blow_up_as_failed():
    km = nsk.KineticModel.from_antimony(BLOW_UP_TOY, units=UNITS)
    conds = [nsk.Condition(t_end=0.5, n_steps=11, label='ok'),
             nsk.Condition(t_end=2, n_steps=21, label='blow-up'),
             nsk.Condition(parameters={'k': 0.25}, t_end=2, n_steps=21)]
    batch = nsk.run_many(km, conds, metrics={'S_end': lambda r: r.final['[S]']})
    m = batch.metrics
    assert list(m.columns) == ['label', 'status', 'error', 'S_end']
    assert list(m['label']) == ['ok', 'blow-up', None]
    assert list(m['status']) == ['ok', 'failed', 'ok']
    assert m['error'][0] is None and m['error'][2] is None
    assert m['error'][1].startswith('KineticSimulationError')
    assert math.isnan(m['S_end'][1])
    assert np.isclose(m['S_end'][0], 2.0, rtol=1e-4)      # 1/(1 - 0.5)
    assert np.isclose(m['S_end'][2], 2.0, rtol=1e-4)      # 1/(1 - 0.25*2)
    assert batch.results[1] is None and batch.results[0].label == 'ok'
    assert batch.n_failed == 1 and list(batch.ok) == [True, False, True]
    assert batch.conditions == conds and len(batch) == 3
    with pytest.raises(KineticSimulationError):
        nsk.run_many(km, conds, on_error='raise')


def test_metric_failures_are_recorded_per_metric():
    km = _feedspike_model()
    cond = nsk.Condition(t_end=10, n_steps=11)
    metrics = {'P_end': _p_end, 'missing': lambda r: r['[nope]'][-1],
               'unrecorded': lambda r: r.final['n_spk'],    # plain Series
               'vector': lambda r: r['[P]'],
               'ratio': lambda r: r['[S]'][-1] * 0 / np.float64(0.0),
               'low_S': lambda r: r.final['[S]'] < 1e3,     # numpy bool
               'flag': lambda r: False}
    batch = nsk.run_many(km, [cond], metrics=metrics)
    row = batch.metrics.iloc[0]
    assert row['status'] == 'failed'
    assert "metric 'missing': KeyError" in row['error']
    assert "metric 'vector': TypeError" in row['error']
    # A KeyError from a plain Series says how to record the column (once:
    # SimulationResult's own message for '[nope]' gets no second hint).
    unrecorded = row['error'].split("metric 'unrecorded': ", 1)[1]
    assert unrecorded.startswith("KeyError: 'n_spk' (column 'n_spk' was not "
                                 "recorded; recorded: ['time', '[S]', '[P]']")
    assert row['error'].count('Condition.selections') == 1
    assert math.isnan(row['unrecorded'])
    # A column returned by mistake is reported by type and shape, not repr.
    assert len(row['error']) < 800 and 'shape (11,)' in row['error']
    assert row['low_S'] == 1.0 and row['flag'] == 0.0   # bools count as 1/0
    assert np.isfinite(row['P_end'])                 # other metrics still count
    assert math.isnan(row['missing']) and math.isnan(row['vector'])
    assert math.isnan(row['ratio'])                  # 0/0 is NaN, not an error
    assert batch.results[0] is not None
    with pytest.raises(KeyError):
        nsk.run_many(km, [cond], metrics={'missing': metrics['missing']},
                     on_error='raise')


# --- the caller's model is never touched (A7) and nothing runs on bad input ------
def test_run_many_runs_on_a_private_copy():
    km = _feedspike_model()
    km._te.k = 0.45                                   # a user edit
    km.simulate(0, 5, 6, reset=True)                  # a dirty state
    before = _model_state(km)
    conds = [nsk.Condition(parameters={'k': k, 'n_max': 2},
                           initial={'env': 1.5}, t_end=48, n_steps=97)
             for k in (0.2, 0.9)]
    batch = nsk.run_many(km, conds, metrics=[SpikeCount()])
    after = _model_state(km)
    assert after[:5] == before[:5] and after[4] is before[4]
    assert after[5] == before[5]
    # Each result equals a direct run on a fresh copy; the metric's
    # selections were recorded on top of the default columns.
    for cond, res in zip(conds, batch.results):
        assert res.col_names == ['time', '[S]', '[P]', 'n_spk']
        assert np.array_equal(res.array, km.copy().run(
            cond, extra_selections=['n_spk']).array)
    assert list(batch.metrics['n_spikes']) == [2.0, 2.0]


@pytest.mark.parametrize('kwargs, error', [
    (dict(conditions=[nsk.Condition(parameters={'nope': 1}, t_end=1)]),
     ValueError),
    (dict(conditions=[nsk.Condition(t_end=1, selections=['[nope]'])]),
     ValueError),
    (dict(conditions=[nsk.Condition(t_end=1, reset=False)]), ValueError),
    (dict(conditions=[nsk.Condition(t_end=1,
                                    reset_kwargs={'mode': 'roadrunner'})]),
     ValueError),
    (dict(conditions=[nsk.Condition(t_end=1,
                                    reset_kwargs={'mode': 'custom'})]),
     ValueError),
    (dict(conditions=[nsk.Condition(t_end=1,
                                    reset_kwargs={'preserve': ['nope']})]),
     ValueError),
    # Preserving state carries the previous run's end state into the next.
    (dict(conditions=[nsk.Condition(t_end=1),
                      nsk.Condition(t_end=1,
                                    reset_kwargs={'preserve': ('n_spk',)})]),
     ValueError),
    (dict(conditions=[nsk.Condition(t_end=1,
                                    reset_kwargs={'preserve': 'env'})]),
     ValueError),
    (dict(conditions=[nsk.Condition(t_end=1,
                                    reset_kwargs={'preserve': ('[S]',)})]),
     ValueError),
    (dict(metrics={'status': _p_end}), ValueError),
    (dict(metrics={'P': 1.0}), TypeError),
    (dict(metrics=[_p_end]), TypeError),             # no `name`
    (dict(metrics=_p_end), TypeError),
    (dict(on_error='ignore'), ValueError),
    (dict(workers=0), ValueError),
    (dict(workers=True), TypeError),
    (dict(keep_results='no'), TypeError),
    (dict(conditions=[{'t_end': 1}]), TypeError),
    (dict(conditions=nsk.Condition(t_end=1)), TypeError),
])
def test_run_many_validates_before_running(monkeypatch, kwargs, error):
    km = _feedspike_model()
    calls = []
    real_run = nsk.KineticModel.run
    monkeypatch.setattr(nsk.KineticModel, 'run',
                        lambda self, *a, **k: calls.append(1)
                        or real_run(self, *a, **k))
    args = dict(conditions=[nsk.Condition(t_end=1)] * 2)
    args.update(kwargs)
    with pytest.raises(error):
        nsk.run_many(km, args.pop('conditions'), **args)
    assert calls == []


def test_batch_preserve_names_settings_only():
    km = _feedspike_model()
    with pytest.raises(ValueError, match=r"conditions \[1\] preserve state "
                                         r"variable\(s\) \['n_spk', 'env'\]"):
        nsk.run_many(km, [nsk.Condition(t_end=1),
                          nsk.Condition(t_end=1, reset_kwargs={
                              'preserve': ('n_spk', 'env', 'n_max')})])
    # A setting (the spike cap) may be preserved.
    batch = nsk.run_many(km, [nsk.Condition(
        t_end=48, n_steps=49, reset_kwargs={'preserve': ('n_max',)})],
        metrics={'P': _p_end})
    assert batch.n_failed == 0
    # A list of plain functions is rejected with a pointer to the dict form.
    with pytest.raises(TypeError, match=r"Pass a mapping instead, e\.g\. "
                                        r"metrics=\{'P_end': fn\}"):
        nsk.run_many(km, [nsk.Condition(t_end=1)], metrics=[_p_end])


def test_run_many_options_and_inputs():
    km = _feedspike_model()
    conds = [nsk.Condition(parameters={'k': 0.2}, t_end=10, n_steps=11),
             nsk.Condition(initial={'[S]': 30}, t_end=10, n_steps=11)]
    batch = nsk.run_many(km, conds, metrics=[SpikeCount()], keep_results=False)
    assert isinstance(batch, nsk.RunBatch)               # top-level exports
    assert all(name in nsk.__all__ for name in (
        'Condition', 'SimulationResult', 'run_many', 'scan', 'RunBatch',
        'ScanResult'))
    assert batch.results == [None, None]
    assert batch.metric_names == ['n_spikes'] and batch.n_failed == 0
    inputs = batch.inputs
    assert list(inputs.columns) == ['k', '[S]']
    assert inputs['k'][0] == 0.2 and math.isnan(inputs['k'][1])
    assert inputs['[S]'][1] == 30.0
    empty = nsk.run_many(km, [], metrics={'P': _p_end})
    assert len(empty) == 0 and list(empty.metrics.columns) == [
        'label', 'status', 'error', 'P']
    assert 'n_failed=0' in repr(batch)


# --- (j) scan -------------------------------------------------------------------
def test_scan_product_grid_matches_direct_runs(close_figures):
    km = _feedspike_model()
    base = nsk.Condition(parameters={'threshold': 8}, t_end=24, n_steps=49,
                         label='FB')
    grid = {'k': [0.2, 0.5], '[S]': [30, 50, 70]}
    before = _model_state(km)
    s = nsk.scan(km, base, grid, metrics={'P_end': _p_end})
    after = _model_state(km)
    assert after[:4] == before[:4] and after[4] is before[4]
    assert after[5] == before[5]                      # the model is untouched
    assert isinstance(s, nsk.ScanResult) and isinstance(s, nsk.RunBatch)
    m = s.metrics
    assert list(m.columns) == ['label', 'k', '[S]', 'status', 'error', 'P_end']
    assert len(m) == 6 and len(s.results) == 6
    assert list(m['k']) == [0.2] * 3 + [0.5] * 3          # first key slowest
    assert list(m['[S]']) == [30.0, 50.0, 70.0] * 2
    assert m['label'][4] == 'FB: k=0.5, [S]=50'
    assert s.placement == {'k': 'parameters', '[S]': 'initial'}
    assert s.metric_names == ['P_end'] and s.variables == ['k', '[S]']
    for (k, s0), cond, res, p_end in zip(
            [(k, s0) for k in grid['k'] for s0 in grid['[S]']],
            s.conditions, s.results, m['P_end']):
        # Grid values are merged over the base condition's own values.
        assert cond.parameters == {'threshold': 8.0, 'k': k}
        assert cond.initial == {'[S]': s0}
        direct = _feedspike_model().run(
            base.with_parameters(k=k).with_initial({'[S]': s0}))
        # A copy (the scan's) vs a fresh load: equal up to ulps.
        assert np.allclose(res.array, direct.array, rtol=1e-12, atol=0)
        assert p_end == res.final['[P]']                  # same run: exact
        assert np.isclose(p_end, direct.final['[P]'], rtol=1e-12, atol=0)
    table = s.pivot('P_end')
    assert table.shape == (2, 3)
    assert list(table.index) == [0.2, 0.5]
    assert list(table.columns) == [30.0, 50.0, 70.0]
    assert table.loc[0.5, 50.0] == m['P_end'][4]
    fig, ax = s.contourf('P_end')
    assert ax.get_xlabel() == 'k' and ax.get_ylabel() == '[S]'
    with pytest.raises(KeyError):
        s.pivot('nope')


def test_scan_modes_placement_and_validation():
    km = _feedspike_model()
    base = nsk.Condition(t_end=10, n_steps=11)
    z = nsk.scan(km, base, {'k': [0.2, 0.4], 'env': [1.0, 2.0]},
                 metrics={'P_end': _p_end}, mode='zip',
                 placement={'env': 'parameters'})
    assert len(z.metrics) == 2 and z.mode == 'zip'
    assert z.placement == {'k': 'parameters', 'env': 'parameters'}
    assert z.conditions[1].parameters == {'k': 0.4, 'env': 2.0}
    auto = nsk.scan(km, base, {'env': [1.0, 2.0], 'n_spk': [0.0]},
                    metrics={'P_end': _p_end})
    assert auto.placement == {'env': 'initial', 'n_spk': 'initial'}
    one = nsk.scan(km, base, {'k': [0.4, 0.2]}, metrics={'P_end': _p_end})
    assert list(one.pivot('P_end').index) == [0.2, 0.4]
    with pytest.raises(ValueError, match='2-variable product'):
        z.contourf('P_end')
    with pytest.raises(ValueError, match='2-variable product'):
        one.contourf('P_end')
    three = nsk.scan(km, base, {'k': [0.2], 'n_max': [1.0], 'threshold': [5.0]},
                     metrics={'P_end': _p_end})
    with pytest.raises(ValueError, match='1- or 2-variable'):
        three.pivot('P_end')
    for bad, error in [
            (dict(grid={'k': [0.1, 0.2], 'env': [1.0]}, mode='zip'),
             ValueError),
            (dict(grid={'k': [0.1, 0.1]}), ValueError),
            (dict(grid={'nope': [0.1]}), ValueError),
            (dict(grid={}), ValueError),
            (dict(grid=[('k', [0.1, 0.2])]), TypeError),     # not a mapping
            (dict(grid={'k': []}), ValueError),
            (dict(grid={'k': 'abc'}), TypeError),
            (dict(grid={'k': [0.1]}, mode='grid'), ValueError),
            (dict(grid={'k': [0.1]}, placement={'env': 'initial'}),
             ValueError),
            (dict(grid={'k': [0.1]}, placement={'k': 'state'}), ValueError),
            (dict(grid={'k': [0.1]}, metrics={'k': _p_end}), ValueError),
            (dict(grid={'[S]': [10.0]}, placement={'[S]': 'parameters'},
                  condition=base.with_initial({'[S]': 5})), ValueError)]:
        args = dict(condition=base, metrics={'P_end': _p_end})
        args.update(bad)
        with pytest.raises(error):
            nsk.scan(km, args.pop('condition'), args.pop('grid'), **args)


    with pytest.raises(ValueError, match='clashes with the batch table'):
        nsk.scan(km, base, {'status': [1.0]}, metrics={'P_end': _p_end})


def _s0(res):
    return res['[S]'][0]


def test_scanning_a_volume_keeps_the_base_concentration():
    # run() sets compartments before concentrations whatever the key order,
    # so a base '[S]' holds in every scanned volume (with_initial appends the
    # grid key after it).
    km = _feedspike_model()
    base = nsk.Condition(initial={'[S]': 50}, t_end=10, n_steps=11)
    s = nsk.scan(km, base, {'env': [1, 2, 4]}, metrics={'S0': _s0})
    assert list(s.conditions[1].initial) == ['[S]', 'env']
    assert list(s.metrics['S0']) == [50.0, 50.0, 50.0]
    # A volume alone still keeps the amounts (concentrations scale).
    alone = nsk.scan(km, nsk.Condition(t_end=10, n_steps=11), {'env': [2]},
                     metrics={'S0': _s0})
    assert list(alone.metrics['S0']) == [25.0]


def test_scan_keeps_a_key_where_the_base_condition_sets_it():
    km = _feedspike_model()
    # 'env' would be auto-placed in initial, 'k' in parameters; a base
    # condition that already sets them elsewhere wins.
    base = nsk.Condition(parameters={'env': 1.5}, initial={'k': 0.2},
                         t_end=10, n_steps=11)
    s = nsk.scan(km, base, {'env': [1.0, 2.0], 'k': [0.3]},
                 metrics={'P_end': _p_end})
    assert s.placement == {'env': 'parameters', 'k': 'initial'}
    for cond, res in zip(s.conditions, s.results):
        assert set(cond.parameters) == {'env'} and set(cond.initial) == {'k'}
        assert cond.initial['k'] == 0.3
        assert np.allclose(res.array, _feedspike_model().run(cond).array,
                           rtol=1e-12, atol=0)
    assert execution._placement_for(km, base, 'env') == 'parameters'
    assert execution._placement_for(km, base, 'k') == 'initial'
    assert execution._placement_for(km, base, 'threshold') == 'parameters'
    assert execution._placement_for(km, base, '[S]') == 'initial'
    # An explicit placement that contradicts the base condition is named.
    with pytest.raises(ValueError, match="contradicts the base condition"):
        nsk.scan(km, base, {'k': [0.3]}, metrics={'P_end': _p_end},
                 placement={'k': 'parameters'})


def test_pivot_rejects_repeated_zip_pairs():
    km = _feedspike_model()
    base = nsk.Condition(t_end=5, n_steps=6)
    z = nsk.scan(km, base, {'k': [0.2, 0.2, 0.4], 'threshold': [5, 5, 6]},
                 metrics={'P_end': _p_end}, mode='zip')
    with pytest.raises(ValueError, match=r'pair once.*\(0\.2, 5\.0\)'):
        z.pivot('P_end')


def test_run_many_rejects_roadrunner_resets_without_touching_the_model():
    km = _feedspike_model()
    km.simulate(0, 5, 6, reset=True)
    before = _model_state(km)
    cond = nsk.Condition(t_end=10, n_steps=11)
    with pytest.raises(ValueError, match="mode='roadrunner'"):
        nsk.run_many(km, [cond, cond.replace(
            reset_kwargs={'mode': 'roadrunner'})])
    after = _model_state(km)
    assert after[:4] == before[:4] and after[4] is before[4]
    assert after[5] == before[5]


def _k_after_a_run(model, k):
    model.run(nsk.Condition(parameters={'k': k}, t_end=1, n_steps=3))
    return (model.get_value('k'), model.get_value('n_spk'))


def test_map_with_model_runs_items_on_one_private_copy():
    km = _feedspike_model()
    km._te.k = 0.45
    before = _model_state(km)
    out = execution._map_with_model(km, _k_after_a_run, [0.1, 0.9],
                                    workers=1)
    assert out == [(0.45, 0.0), (0.45, 0.0)]   # the copy carried the edit
    after = _model_state(km)
    assert after[:4] == before[:4] and after[4] is before[4]
    assert after[5] == before[5]
    assert execution._map_with_model(km, _k_after_a_run, [], workers=2) == []


def _reset_bumping_the_cap(km, bump=True, **kwargs):
    """Like the shipped model's reset_spike_cap: changes a setting (the spike
    cap) only for some reset kwargs."""
    km.reset(mode='baseline')
    if bump:
        km._te['n_max'] = 10


def _n_spk_end(res):
    return res.final['n_spk']


def test_a_setting_changed_by_reset_func_does_not_leak_into_later_rows():
    km = _feedspike_model(reset=_reset_bumping_the_cap)
    km._te['n_max'] = 1
    keep = nsk.Condition(t_end=30, n_steps=31, reset_kwargs={'bump': False},
                         selections=('n_spk',))
    bump = keep.replace(reset_kwargs={'bump': True})
    metrics = {'n': _n_spk_end}
    a = nsk.run_many(km, [keep, bump, keep], metrics=metrics)
    b = nsk.run_many(km, [keep, keep, bump], metrics=metrics)
    assert list(a.metrics['status']) == list(b.metrics['status']) == ['ok'] * 3
    n_a, n_b = list(a.metrics['n']), list(b.metrics['n'])
    assert n_a == [n_b[0], n_b[2], n_b[1]]          # same per condition
    assert n_a[0] == n_a[2] == _n_spk_end(km.run(keep)) == 1.0
    assert n_a[1] > 1.0                              # the bumped cap
    assert km._te['n_max'] == 1                      # the caller's cap
    # The pool path (initializer + per-chunk call), run in this process.
    fn = functools.partial(execution._execute, metrics=metrics,
                           extra_selections=(), keep_results=False,
                           on_error='raise')
    try:
        execution._worker_init(km._blueprint(), fn)
        rows = execution._worker_call([keep, bump, keep])
    finally:
        execution._WORKER.clear()
    assert [row[2]['n'] for row in rows] == n_a


# --- saving batches -------------------------------------------------------------
def _assert_same_batch(back, batch):
    assert type(back) is type(batch)
    pd.testing.assert_frame_equal(back.metrics, batch.metrics,
                                  check_exact=True)
    assert back.conditions == batch.conditions
    assert back.results == [None] * len(batch)


def test_run_batch_csv_round_trip_is_exact(tmp_path):
    km = nsk.KineticModel.from_antimony(BLOW_UP_TOY, units=UNITS)
    conds = [nsk.Condition(t_end=0.5, n_steps=11, label='ok, "quoted" é'),
             nsk.Condition(t_end=2, n_steps=21, label='blow-up'),
             nsk.Condition(parameters={'k': 0.1}, times=[0, 0.3, 1 / 3],
                           selections=['[S]'],
                           reset_kwargs={'preserve': ['k']})]
    metrics = {'S_end': lambda r: r.final['[S]'],
               'inv': lambda r: 1 / (r.final['[S]'] - r.final['[S]']),
               'neg_zero': lambda r: -0.0, 'tiny': lambda r: 1e-300}
    batch = nsk.run_many(km, conds, metrics=metrics)
    assert batch.metrics['inv'][0] == math.inf       # a non-finite value
    assert list(batch.metrics['status']) == ['ok', 'failed', 'ok']
    path = tmp_path / 'batch.csv'
    batch.to_csv(path)
    back = nsk.RunBatch.from_csv(path)
    _assert_same_batch(back, batch)
    assert math.copysign(1, back.metrics['neg_zero'][0]) == -1
    assert back.metrics['error'][0] is None and back.metric_names == list(
        metrics)
    head = pd.read_csv(path, dtype=str, keep_default_na=False)
    assert list(head.columns) == ['row', 'quantity', 'key', 'value']
    assert set(head['quantity']) == {'meta', 'label', 'status', 'error',
                                     'metric', 'condition'}
    # an empty batch keeps its columns and dtypes
    empty = nsk.run_many(km, [], metrics={'S_end': metrics['S_end']})
    empty.to_csv(tmp_path / 'empty.csv')
    _assert_same_batch(nsk.RunBatch.from_csv(tmp_path / 'empty.csv'), empty)
    # files that are not batches, or a batch read as a scan, are rejected
    with pytest.raises(ValueError, match='holds a RunBatch'):
        nsk.ScanResult.from_csv(path)
    batch.results[0].to_csv(tmp_path / 'run.csv')
    with pytest.raises(ValueError, match='not written by RunBatch'):
        nsk.RunBatch.from_csv(tmp_path / 'run.csv')


def test_scan_result_csv_round_trip_is_exact(tmp_path):
    km = _feedspike_model()
    base = nsk.Condition(parameters={'threshold': 8}, t_end=10, n_steps=11,
                         label='FB')
    s = nsk.scan(km, base, {'k': [0.2, 0.5], '[S]': [30, 50]},
                 metrics=[SpikeCount()], keep_results=False)
    s.to_csv(tmp_path / 'scan.csv')
    back = nsk.RunBatch.from_csv(tmp_path / 'scan.csv')   # any class reads it
    _assert_same_batch(back, s)
    assert back.grid == s.grid == {'k': (0.2, 0.5), '[S]': (30.0, 50.0)}
    assert back.mode == 'product' and back.placement == s.placement
    assert back.variables == ['k', '[S]'] and back.metric_names == ['n_spikes']
    pd.testing.assert_frame_equal(back.pivot('n_spikes'), s.pivot('n_spikes'),
                                  check_exact=True)
    back.to_csv(tmp_path / 'again.csv')
    assert ((tmp_path / 'again.csv').read_text()
            == (tmp_path / 'scan.csv').read_text())


# --- (k) process pool -----------------------------------------------------------
def _fork_only():
    if execution._mp_context().get_start_method() != 'fork':
        pytest.skip('needs the fork start method')


def test_two_workers_match_the_serial_batch():
    # A spawn pool (macOS, Windows) re-imports nskinetics in every worker
    # (~20 s), too slow for this suite; and the lambda metric needs fork.
    _fork_only()
    km = nsk.KineticModel.from_antimony(BLOW_UP_TOY, units=UNITS)
    conds = [nsk.Condition(parameters={'k': k}, t_end=2, n_steps=21,
                           label=f'k={k}')
             for k in (0.1, 0.25, 1.0, 0.3, 0.45)]      # k=1 blows up at t=1
    metrics = {'S_end': lambda r: r.final['[S]']}   # lambda: needs fork
    serial = nsk.run_many(km, conds, metrics=metrics)
    parallel = nsk.run_many(km, conds, metrics=metrics, workers=2)
    pd.testing.assert_frame_equal(parallel.metrics, serial.metrics,
                                  check_exact=True)
    assert list(serial.metrics['status']).count('failed') == 1
    for a, b in zip(serial.results, parallel.results):
        assert (a is None) == (b is None)
        if a is not None:
            assert np.array_equal(a.array, b.array)
            assert a.condition == b.condition
    with pytest.raises(KineticSimulationError):
        nsk.run_many(km, conds, metrics=metrics, workers=2, on_error='raise')


def test_workers_reject_unpicklable_models_up_front(monkeypatch):
    km = nsk.KineticModel.from_antimony(BLOW_UP_TOY, units=UNITS,
                                        reset=lambda m, **kw: m._te.reset())
    cond = nsk.Condition(t_end=0.5, n_steps=3)

    def no_pool(*args, **kwargs):
        raise AssertionError('a pool was created')

    monkeypatch.setattr(execution, 'ProcessPoolExecutor', no_pool)
    with pytest.raises(ValueError, match='reset_func'):
        nsk.run_many(km, [cond, cond], workers=2)
    # workers=1 never needs pickling, nor a pool.
    assert nsk.run_many(km, [cond], workers=1).n_failed == 0
    # Under spawn, metrics travel by pickle too.
    km.reset_func = _reset_with_lambda_attribute
    monkeypatch.setattr(execution, '_mp_context',
                        lambda: multiprocessing.get_context('spawn'))
    with pytest.raises(ValueError, match="Metric 'S'"):
        nsk.run_many(km, [cond, cond], workers=2,
                     metrics={'S': lambda r: r.final['[S]']})


def _build_fails(self):
    raise RuntimeError('cannot rebuild this model')


def test_worker_build_failure_surfaces_the_real_error(monkeypatch):
    _fork_only()
    km = nsk.KineticModel.from_antimony(BLOW_UP_TOY, units=UNITS)
    cond = nsk.Condition(t_end=0.5, n_steps=3)
    # Forked workers inherit the patch; the parent never builds a model.
    monkeypatch.setattr(_ModelBlueprint, 'build', _build_fails)
    with pytest.raises(RuntimeError, match='cannot rebuild this model'):
        nsk.run_many(km, [cond, cond], workers=2)


def test_a_crashed_worker_raises_a_runtime_error_with_a_hint():
    _fork_only()
    km = nsk.KineticModel.from_antimony(BLOW_UP_TOY, units=UNITS)
    cond = nsk.Condition(t_end=0.5, n_steps=3)
    with pytest.raises(RuntimeError, match='workers=1'):
        nsk.run_many(km, [cond, cond], workers=2,
                     metrics={'die': lambda r: os._exit(1)})


def test_start_method_prefers_fork_on_linux():
    if 'fork' in multiprocessing.get_all_start_methods() and \
            not execution.sys.platform.startswith(('darwin', 'win')):
        assert execution._mp_context().get_start_method() == 'fork'
