# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Condition, KineticModel.run and SimulationResult (F1 items 4-6; spec
tests d, e, f)."""
import copy
import dataclasses
import json
import pickle

import matplotlib
import numpy as np
import pandas as pd
import pytest
import tellurium as te

import nskinetics as nsk
from nskinetics.exceptions import KineticSimulationError

matplotlib.use('Agg')

# The FeedSpike toy of the verified reset bug: S is consumed into P and
# refilled to `target` (volume grows) whenever it drops below `threshold`, at
# most `n_max` times.
FEEDSPIKE_TOY = """
model fb()
  compartment env = 1; species S in env, P in env;
  S = 50; P = 0; k = 0.5; threshold = 5; target = 50; feed_conc = 500;
  n_max = 3; n_spk = 0; last_vol = 0; tot_vol = 0; curr_env := env;
  J: S => P; k*S*env;
end
"""

# dS/dt = k*S^2 with S(0) = 1 blows up at t = 1/k.
BLOW_UP_TOY = """
model blow()
  species S; S = 1; k = 1;
  J: => S; k*S^2;
end
"""

# The FeedSpike toy plus a boundary species B whose concentration the rate
# law reads; the spikes grow env from 1 to 1.331, diluting [B].
FEEDSPIKE_B_TOY = (FEEDSPIKE_TOY
                   .replace('species S in env, P in env;',
                            'species S in env, P in env; species $B in env;')
                   .replace('S = 50; P = 0;', 'S = 50; P = 0; B = 2;')
                   .replace('J: S => P; k*S*env;', 'J: S => P; k*S*env*B/2;'))

# Initial assignments that read k (k2 a parameter, S a species), and a plain
# parameter k3.
IA_TOY = """
model ia()
  compartment env = 1; species S in env, P in env;
  k = 0.5; k3 = 1; k2 = 2*k; S = 20*k; P = 0;
  J: S => P; k*k3*S;
end
"""

UNITS = {'time': 'h', 'conc': 'g/L'}
SELS = ('time', '[S]', '[P]', 'env', 'n_spk')


def _feedspike_model(reset=None, antimony=FEEDSPIKE_TOY):
    km = nsk.KineticModel(te.loadAntimonyModel(antimony), units=UNITS,
                          reset=reset)
    spike = nsk.FeedSpike('S', 'S <= threshold', 'target', 'feed_conc',
                          volume_var='env', max_count='n_max',
                          count_var='n_spk', last_vol_var='last_vol',
                          tot_vol_var='tot_vol')
    for event in spike.expand():
        km.add_event(event)
    km.compile_events()
    return km


def _cond(**kw):
    kw.setdefault('t_end', 48)
    kw.setdefault('n_steps', 481)
    return nsk.Condition(**kw)


@pytest.fixture
def close_figures():
    """Close every figure after the test, also when an assertion fails."""
    import matplotlib.pyplot as plt
    try:
        yield
    finally:
        plt.close('all')


# --- Condition: construction and validation (d) ---------------------------------
def test_condition_needs_exactly_one_of_t_end_and_times():
    with pytest.raises(ValueError, match='exactly one'):
        nsk.Condition(t_end=10, times=[0, 5, 10])
    with pytest.raises(ValueError, match='exactly one'):
        nsk.Condition()


@pytest.mark.parametrize('times', [[0, 2, 1], [0, 1, 1, 2], [0], [0, np.nan],
                                   [[0, 1], [2, 3]]])
def test_condition_rejects_bad_time_grids(times):
    with pytest.raises(ValueError):
        nsk.Condition(times=times)


@pytest.mark.parametrize('kwargs, error', [
    (dict(t_end=10, n_steps=1), ValueError),
    (dict(t_end=10, n_steps=10.5), TypeError),
    (dict(t_end=10, n_steps=True), TypeError),
    (dict(t0=5, t_end=5), ValueError),
    (dict(t_end=np.inf), ValueError),
    (dict(t_end=10, initial={'[S]': 'fifty'}), TypeError),
    (dict(t_end=10, parameters={'k': np.nan}), ValueError),
    (dict(t_end=10, parameters={'k': True}), TypeError),
    (dict(t_end=10, initial=[('k', 1.0)]), TypeError),
    (dict(t_end=10, initial={'k': 1}, parameters={'k': 2}), ValueError),
    (dict(t_end=10, selections=['time', '[S]', '[S]']), ValueError),
    (dict(t_end=10, selections='[S]'), TypeError),
    (dict(t_end=10, selections=[]), TypeError),
    (dict(t_end=10, reset='yes'), TypeError),
    (dict(t_end=10, reset_kwargs=['mode']), TypeError),
    (dict(t_end=10, label=3), TypeError),
    (dict(t0=1, times=[0, 1, 2]), ValueError),
    # an explicit t0 = 0 before data-only output times is not moved silently
    (dict(t0=0, times=[1, 2]), ValueError),
    (dict(t_end=10, initial={'S': 1, '[S]': 2}), ValueError),
    # Equal, hash-equal conditions would otherwise start at 1 or 2.
    (dict(t_end=10, initial={'init(S)': 1, 'init([S])': 2}), ValueError),
    (dict(t_end=10, initial={'init([S])': 2, 'init(S)': 1}), ValueError),
])
def test_condition_validation(kwargs, error):
    with pytest.raises(error):
        nsk.Condition(**kwargs)


def test_condition_is_immutable_and_normalized():
    initial = {'[S]': 50}
    cond = nsk.Condition(initial=initial, t_end=48, selections=['time', '[S]'],
                         reset_kwargs={'preserve': ['n_max'],
                                       'nested': {'a': [1, 2]}})
    initial['[S]'] = 1.0                    # the caller's dict is copied
    assert cond.initial == {'[S]': 50.0}
    assert isinstance(cond.initial['[S]'], float) and cond.t_end == 48.0
    with pytest.raises(TypeError):
        cond.initial['[S]'] = 2.0
    with pytest.raises(dataclasses.FrozenInstanceError):
        cond.t_end = 10.0
    assert cond.selections == ('time', '[S]')
    # Any sequence of strings works for selections, as for times; numpy
    # strings are stored as plain str.
    arr = nsk.Condition(t_end=48, selections=np.array(['time', '[S]']))
    assert arr.selections == cond.selections and arr == cond.replace(
        initial={}, reset_kwargs={})
    assert all(type(s) is str for s in arr.selections)
    assert json.dumps(arr.to_dict())
    for bad in ('time', {'time': 1}, b'time', 3):
        with pytest.raises(TypeError, match='sequence of strings'):
            nsk.Condition(t_end=48, selections=bad)
    assert cond.reset_kwargs['preserve'] == ('n_max',)
    assert cond.reset_kwargs['nested']['a'] == (1, 2)
    with pytest.raises(TypeError):
        cond.reset_kwargs['nested']['b'] = 1
    # With times, integration starts at times[0].
    grid = nsk.Condition(times=[1, 2.5, 4])
    assert grid.times == (1.0, 2.5, 4.0) and grid.t0 == 1.0
    assert nsk.Condition(t0=1, times=[1, 2]).t0 == 1.0
    assert nsk.Condition(t_end=5).t0 == 0.0          # t0=None -> 0.0
    assert nsk.Condition.from_dict(grid.to_dict()) == grid
    with pytest.raises(TypeError):
        nsk.Condition(0, 10)                # keyword-only


def test_replace_and_merging_helpers():
    base = _cond(initial={'[S]': 50}, parameters={'k': 0.4, 'threshold': 5},
                 label='base')
    merged = base.with_parameters({'k': 0.8}, n_max=2)
    assert merged.parameters == {'k': 0.8, 'threshold': 5.0, 'n_max': 2.0}
    assert merged.initial == base.initial and merged.label == 'base'
    assert base.parameters == {'k': 0.4, 'threshold': 5.0}   # unchanged
    assert base.with_initial({'env': 2}).initial == {'[S]': 50.0, 'env': 2.0}
    # replace() replaces wholesale.
    assert base.replace(parameters={'k': 0.1}).parameters == {'k': 0.1}
    # t_end and times are switched, not combined.
    on_grid = base.replace(times=[0, 1, 2])
    assert on_grid.t_end is None and on_grid.times == (0.0, 1.0, 2.0)
    shifted = on_grid.replace(times=[3, 4])
    assert shifted.t0 == 3.0
    back = shifted.replace(t_end=10)
    assert back.times is None and back.t_end == 10.0 and back.t0 == 3.0
    with pytest.raises(TypeError, match='no field'):
        base.replace(tend=10)
    with pytest.raises(ValueError):
        base.replace(t_end=-1)
    with pytest.raises(ValueError, match='both'):
        base.with_parameters({'[S]': 10})


def test_mappings_may_be_pandas_series():
    """A fit's estimates and sample rows are pandas Series: every entry
    point takes them like a dict, and a non-mapping raises the same clean
    TypeError everywhere."""
    base = _cond(initial={'[S]': 50}, parameters={'k': 0.4})
    estimates = pd.Series({'k': 0.8, 'n_max': 2.0}, name='estimate')
    merged = base.with_parameters(estimates)
    assert merged.parameters == {'k': 0.8, 'n_max': 2.0}
    assert base.with_parameters(estimates, k=0.9).parameters['k'] == 0.9
    assert base.with_initial(pd.Series({'[X]': 0.3})).initial == {
        '[S]': 50.0, '[X]': 0.3}
    # an empty Series is a no-op
    assert base.with_parameters(pd.Series(dtype=float)) == base
    assert base.with_initial(pd.Series(dtype=float)) == base
    cond = nsk.Condition(parameters=pd.Series({'k': 0.5}), t_end=1,
                         initial=pd.Series({'[S]': 2.0}))
    assert cond.parameters == {'k': 0.5} and cond.initial == {'[S]': 2.0}
    assert base.replace(parameters=estimates).parameters == dict(estimates)
    for bad in (5, [('k', 1.0)]):
        with pytest.raises(TypeError, match='must be a mapping'):
            base.with_parameters(bad)
        with pytest.raises(TypeError, match='must be a mapping'):
            base.with_initial(bad)
    with pytest.raises(TypeError, match='keys must be non-empty strings'):
        base.with_parameters(pd.Series([1.0, 2.0]))


def test_condition_serialization_round_trips():
    cond = nsk.Condition(initial={'[S]': 50.5}, parameters={'k': 0.1 + 0.2},
                         times=[0, 0.5, 48], selections=['time', '[P]'],
                         reset_kwargs={'preserve': ['n_max'], 'flag': True},
                         label='FB1 µ\n')
    d = cond.to_dict()
    assert d['times'] == [0.0, 0.5, 48.0]
    assert d['reset_kwargs'] == {'preserve': ['n_max'], 'flag': True}
    again = nsk.Condition.from_dict(json.loads(json.dumps(d)))
    assert again == cond and hash(again) == hash(cond)
    assert pickle.loads(pickle.dumps(cond)) == cond
    assert copy.deepcopy(cond) == cond
    assert cond != cond.with_parameters(k=0.3)
    assert eval(repr(cond), {'Condition': nsk.Condition}) == cond
    with pytest.raises(ValueError, match='Unknown'):
        nsk.Condition.from_dict({**d, 'bogus': 1})
    with pytest.raises(TypeError):
        nsk.Condition(t_end=1, reset_kwargs={'obj': object()}).to_dict()


# --- KineticModel.run (d, e) ----------------------------------------------------
def _snapshot(km):
    r = km._te
    return (r.getCurrentSBML(), km.results_array, list(km.events),
            km._events_compiled)


@pytest.mark.parametrize('bad', [
    dict(parameters={'nope': 1.0}),
    dict(initial={'[nope]': 1.0}),
    dict(initial={'curr_env': 2.0}),          # assignment rule
    dict(initial={'init(curr_env)': 2.0}),    # ... also as init(x)
    dict(parameters={'J': 1.0}),              # reaction rate (read-only)
    dict(initial={"S'": 1.0}),                # rate of change (read-only)
    dict(selections=['time', '[S]', 'nope']),
    dict(reset_kwargs={'mode': 'bogus'}),
    dict(reset_kwargs={'preserve': ['nope']}),
])
def test_run_rejects_bad_ids_before_touching_state(bad):
    km = _feedspike_model()
    km.run(_cond(parameters={'k': 0.3}))     # leave a dirty, non-baseline state
    before = _snapshot(km)
    with pytest.raises(ValueError):
        km.run(_cond(**bad))
    after = _snapshot(km)
    assert after[0] == before[0]              # byte-identical model state
    assert after[1] is before[1] and after[2:] == before[2:]


def test_id_validation_follows_a_replaced_roadrunner():
    km = nsk.KineticModel.from_antimony(BLOW_UP_TOY, units=UNITS)
    km.run(_cond(t_end=0.5, parameters={'k': 0.3}))  # fills the key cache
    km._te = te.loadAntimonyModel(
        'model m(); species A; A = 1; kd = 0.1; J: A => ; kd*A; end')
    res = km.run(_cond(t_end=1, n_steps=3, parameters={'kd': 0.2}))
    assert res.col_names == ['time', '[A]']
    with pytest.raises(ValueError, match="'k'"):
        km.run(_cond(t_end=1, parameters={'k': 0.3}))


def test_run_lists_every_unknown_id():
    km = _feedspike_model()
    with pytest.raises(ValueError, match=r"'nope'.*'\[nada\]'|'\[nada\]'.*'nope'"):
        km.run(_cond(parameters={'nope': 1}, initial={'[nada]': 1}))
    with pytest.raises(TypeError):
        km.run({'t_end': 10})


def test_run_restores_parameters_and_fills_the_legacy_slot():
    km = _feedspike_model()
    km._te.k = 0.45                           # a persistent user edit
    res = km.run(_cond(initial={'[S]': 40, 'n_max': 2}, parameters={'k': 0.8},
                       selections=SELS))
    assert km._te.k == 0.45                   # parameter override restored
    assert km._te.n_max == 3.0                # non-state `initial` restored
    assert res['n_spk'][-1] == 2.0            # ...but applied during the run
    assert res['[S]'][0] == 40.0
    # Legacy slot filled with the same trajectory (but not aliased).
    assert km.results_col_names == list(SELS)
    assert np.array_equal(km.results_array, res.array)
    pd.testing.assert_frame_equal(km.results_df, res.df)
    assert km.results_df is not res.df
    # State is left at the end of the run, as after simulate().
    assert km._te['n_spk'] == 2.0
    # The model's default k is used without an override.
    ref = _feedspike_model()
    ref._te.k = 0.45
    assert np.array_equal(
        km.run(_cond(selections=SELS)).array,
        ref.simulate(0, 48, 481, list(SELS), reset=True))


def test_condition_initial_volume_applies_to_that_run_only():
    # Spec test (b): Condition.initial={'env': 2} is the start of one run; the
    # next run starts from the baseline again.
    km = _feedspike_model()
    base = km.run(_cond(selections=SELS))
    big = km.run(_cond(initial={'env': 2}, selections=SELS))
    assert big['env'][0] == 2.0 and base['env'][0] == 1.0
    # RoadRunner keeps amounts when a volume is set: S is diluted ...
    assert base['[S]'][0] == 50.0 and big['[S]'][0] == 25.0
    # ... unless a concentration is given too: run() sets volumes first, so
    # it holds at t0 whatever the key order.
    both = km.run(_cond(initial={'env': 2, '[S]': 50}, selections=SELS))
    assert both['env'][0] == 2.0 and both['[S]'][0] == 50.0
    assert not np.allclose(big.array, base.array)
    after = km.run(_cond(selections=SELS))
    assert after['env'][0] == 1.0
    assert np.array_equal(after.array, base.array)


def test_initial_key_order_does_not_change_the_run():
    km = _feedspike_model()
    a = _cond(initial={'[S]': 50, 'env': 2}, selections=SELS)
    b = _cond(initial={'env': 2, '[S]': 50}, selections=SELS)
    assert a == b and hash(a) == hash(b)
    ra, rb = km.run(a), km.run(b)
    assert np.array_equal(ra.array, rb.array)
    assert ra['[S]'][0] == 50.0 and ra['env'][0] == 2.0
    # init(x) entries come first, so a current value is never clobbered.
    sbml = km._te.getSBML()
    c = _cond(initial={'[S]': 30, 'init([S])': 7}, selections=SELS)
    d = _cond(initial={'init([S])': 7, '[S]': 30}, selections=SELS)
    rc, rd = km.run(c), km.run(d)
    assert np.array_equal(rc.array, rd.array) and rc['[S]'][0] == 30.0
    assert km._te.getSBML() == sbml
    # parameters follow the same fixed order (compartments before
    # concentrations), so equal, hash-equal conditions run the same.
    e = _cond(parameters={'env': 2, '[S]': 50}, selections=SELS)
    f = _cond(parameters={'[S]': 50, 'env': 2}, selections=SELS)
    assert e == f and hash(e) == hash(f)
    re_, rf = km.run(e), km.run(f)
    assert np.array_equal(re_.array, rf.array)
    assert rf['[S]'][0] == 50.0 and rf['env'][0] == 2.0
    with pytest.raises(ValueError, match='parameters sets'):
        nsk.Condition(t_end=1, parameters={'S': 1, '[S]': 2})


def test_init_of_a_rule_defined_id_is_rejected_everywhere():
    km = _feedspike_model()
    km.run(_cond(parameters={'k': 0.3}))
    before = _snapshot(km)
    for bad in ('init(curr_env)', 'curr_env'):
        with pytest.raises(ValueError, match='assignment rules'):
            km.run(_cond(initial={bad: 2.0}))
        with pytest.raises(ValueError, match='assignment rules'):
            with km.overrides({bad: 2.0}):
                pass
        # preserve takes no init(...) selection at all.
        with pytest.raises(ValueError, match='assignment rules'
                           if bad == 'curr_env' else 'current-value'):
            km.reset(preserve=(bad,))
    after = _snapshot(km)
    assert after[0] == before[0] and after[1] is before[1]


# --- init() forms that would rewrite the SBML for good ---------------------------
# A compartment defined by an initial assignment.
IA_ENV_TOY = """
model iaenv()
  compartment env; species S in env;
  k = 0.5; env = 2*k; S = 10;
  J: S => ; 0.1*S;
end
"""

# A species declared with an initial amount.
SUBSTANCE_TOY = """
model sub()
  compartment env = 1; substanceOnly species S in env;
  S = 5;
  J: S => ; 0.1*S;
end
"""


def _sbml_state(km):
    return km._te.getSBML(), sorted(km._te.getInitialAssignmentIds())


def test_init_of_an_initial_assigned_id_is_rejected():
    # RoadRunner deletes the initial assignment on an init() write, and
    # writing the old value back does not restore it.
    km = nsk.KineticModel.from_antimony(IA_TOY, units=UNITS)
    before = _sbml_state(km)
    assert before[1] == ['S', 'k2']
    for bad in ({'init(S)': 7.0}, {'init([S])': 7.0}):
        with pytest.raises(ValueError, match='initial assignment'):
            km.run(_cond(initial=bad, t_end=1, n_steps=3))
        with pytest.raises(ValueError, match='initial assignment'):
            with km.overrides(bad):
                pass
    # A parameter too: a ValueError up front, not RoadRunner's RuntimeError
    # half-way through the run.
    with pytest.raises(ValueError, match=r"initial assignment \['init\(k2\)'\]"):
        km.run(_cond(parameters={'init(k2)': 3.0}, t_end=1, n_steps=3))
    assert _sbml_state(km) == before
    # init() of the inputs is the supported route; it leaves the SBML alone.
    res = km.run(_cond(initial={'init(k)': 1.0}, t_end=1, n_steps=3))
    assert res['[S]'][0] == 20.0
    assert _sbml_state(km) == before
    # A compartment defined by an initial assignment.
    km = nsk.KineticModel.from_antimony(IA_ENV_TOY, units=UNITS)
    before = _sbml_state(km)
    assert before[1] == ['env']
    with pytest.raises(ValueError, match='initial assignment'):
        km.run(_cond(initial={'init(env)': 3.0}, t_end=1, n_steps=3))
    assert _sbml_state(km) == before


def test_init_of_a_species_must_use_its_declared_form():
    # init(S) of a concentration-declared species would redeclare it as an
    # initial amount, so a later volume change would start it elsewhere.
    km = _feedspike_model()
    r = km._te
    sbml = r.getSBML()
    with pytest.raises(ValueError, match=r"use 'init\(\[S\]\)'"):
        km.run(_cond(initial={'init(S)': 30.0}))
    with pytest.raises(ValueError, match=r"use 'init\(\[S\]\)'"):
        with km.overrides({'init(S)': 30.0}):
            pass
    with pytest.raises(ValueError, match=r"use 'init\(\[S\]\)'"):
        nsk.run_many(km, [_cond(initial={'init(S)': 30.0})])
    assert r.getSBML() == sbml
    # The declared form leaves the SBML alone, so batches do not depend on
    # what ran before.
    first = _cond(initial={'init([S])': 30.0}, selections=SELS)
    second = _cond(initial={'init(env)': 2.0}, selections=SELS)
    assert km.run(first)['[S]'][0] == 30.0
    assert r.getSBML() == sbml
    batch = nsk.run_many(km, [first, second],
                         metrics={'S0': lambda res: res['[S]'][0]})
    assert list(batch.metrics['S0']) == [30.0, 50.0]
    assert _feedspike_model().run(second)['[S]'][0] == 50.0
    # A species declared with an initial amount takes init(S), not init([S]).
    km = nsk.KineticModel.from_antimony(SUBSTANCE_TOY, units=UNITS)
    sbml = km._te.getSBML()
    with pytest.raises(ValueError, match=r"use 'init\(S\)'"):
        km.run(_cond(initial={'init([S])': 3.0}, t_end=1, n_steps=3))
    res = km.run(_cond(initial={'init(S)': 3.0}, t_end=1, n_steps=3,
                       selections=['time', 'S']))
    assert res['S'][0] == 3.0
    assert km._te.getSBML() == sbml


def test_preserve_rejects_init_selections():
    # Writing an init() value back would re-initialize the model and revert
    # every parameter edit; a reset never changes init() values anyway.
    km = _feedspike_model()
    r = km._te
    r['k'] = 0.4
    for preserve in (('init(n_max)',), 'init(n_max)'):
        with pytest.raises(ValueError, match='current-value selections'):
            km.reset(preserve=preserve)
    cond = _cond(t_end=1, n_steps=2,
                 reset_kwargs={'preserve': ['init(threshold)']})
    with pytest.raises(ValueError, match='current-value selections'):
        km.run(cond)
    with pytest.raises(ValueError, match='current-value selections'):
        nsk.run_many(km, [cond])
    assert r['k'] == 0.4


# --- what run() / overrides() restore ------------------------------------------
IA_SELS = ['time', '[S]', '[P]', 'k', 'k2', 'k3']


def _ia_model():
    km = nsk.KineticModel.from_antimony(IA_TOY, units=UNITS)
    km._te['k'] = 0.4                      # user edits that must survive
    km._te['k3'] = 2.0
    return km


@pytest.mark.parametrize('where', ['initial', 'parameters'])
def test_init_entries_define_the_run_and_keep_settings(where):
    km = _ia_model()
    sbml = km._te.getSBML()
    cond = nsk.Condition(**{where: {'init(k)': 0.7}}, t_end=2, n_steps=21,
                         selections=IA_SELS)
    res = km.run(cond)
    # The run starts as a model defined with k = 0.7 would (initial
    # assignments re-evaluated), with the user's other edit (k3) in place.
    first = res.df.iloc[0]
    assert (first['k'], first['k2'], first['[S]'], first['k3']) == (
        0.7, 1.4, 14.0, 2.0)
    fresh = nsk.KineticModel.from_antimony(
        IA_TOY.replace('k = 0.5', 'k = 0.7').replace('k3 = 1', 'k3 = 2'),
        units=UNITS)
    assert np.array_equal(res.array, fresh.run(cond.replace(
        initial={}, parameters={})).array)
    # Nothing leaks: the user's edits, init(k), k2 and the SBML are restored.
    r = km._te
    assert (r['k'], r['k3'], r['init(k)'], r['k2']) == (0.4, 2.0, 0.5, 1.0)
    assert r.getSBML() == sbml
    plain = nsk.Condition(t_end=2, n_steps=21, selections=IA_SELS)
    ref = nsk.KineticModel.from_antimony(IA_TOY, units=UNITS)
    ref._te['k'], ref._te['k3'] = 0.4, 2.0
    assert np.array_equal(km.run(plain).array, ref.run(plain).array)


def test_parameters_hold_alongside_init_entries():
    km = _ia_model()
    res = km.run(nsk.Condition(initial={'init(k)': 0.7},
                               parameters={'k3': 3.0, 'k2': 5.0},
                               t_end=2, n_steps=21, selections=IA_SELS))
    assert res['k3'][0] == 3.0 and res['k2'][0] == 5.0 and res['k'][0] == 0.7
    assert km._te['k3'] == 2.0 and km._te['k2'] == 1.0


def test_init_entries_do_not_make_batches_order_dependent():
    km = _ia_model()
    plain = nsk.Condition(t_end=2, n_steps=21, label='plain')
    initk = plain.with_initial({'init(k)': 0.7}).replace(label='init(k)')
    batch = nsk.run_many(km, [plain, initk, plain],
                         metrics={'P': lambda res: res.final['[P]']})
    p = list(batch.metrics['P'])
    assert p[0] == p[2] != p[1]
    assert p[0] == km.copy().run(plain).final['[P]']


INTEGRATOR_TOY = """
model m()
  S = 50; P = 0; k = 0.5; th = 1;
  J: S => P; k*S;
end
"""


def test_init_entries_keep_the_integrator_settings():
    # RoadRunner 2.10 resets rk45's epsilon (1e-3 -> 1e-12) and gillespie's
    # seed whenever an init(...) write re-initializes the model.
    km = nsk.KineticModel.from_antimony(INTEGRATOR_TOY, units=UNITS)
    km._te.setIntegrator('rk45')
    km._te.getIntegrator().setValue('epsilon', 1e-3)

    def integrator():
        return km._te.getIntegrator()

    plain = nsk.Condition(t_end=10, n_steps=11)
    initc = plain.with_initial({'init(th)': 2})
    km.run(initc)
    assert integrator().getName() == 'rk45'
    assert integrator().getValue('epsilon') == 1e-3
    with km.overrides({'init(th)': 2}):
        assert integrator().getValue('epsilon') == 1e-3
    assert integrator().getValue('epsilon') == 1e-3
    batch = nsk.run_many(km, [plain, initc, plain],
                         metrics={'P': lambda res: res.final['[P]']})
    p = list(batch.metrics['P'])
    assert p[0] == p[1] == p[2]
    assert np.isclose(p[0], 50 * (1 - np.exp(-5)), rtol=1e-2)
    # gillespie's seed survives too.
    km = nsk.KineticModel.from_antimony(INTEGRATOR_TOY, units=UNITS)
    km._te.setIntegrator('gillespie')
    km._te.getIntegrator().setValue('seed', 1234)
    km.run(initc)
    assert km._te.getIntegrator().getValue('seed') == 1234


def test_run_restores_init_entries_but_reinitializes_the_state():
    # Pins the documented behaviour: restoring init(x) re-initializes the
    # model, so it is not left at its end-of-run state.
    km = _feedspike_model()
    r = km._te
    sbml = r.getSBML()
    res = km.run(_cond(initial={'init([S])': 7.0}, t_end=1, n_steps=11,
                       selections=SELS))
    assert res['[S]'][0] == 7.0 and res.final['n_spk'] == 1.0
    assert r['init([S])'] == 50.0
    assert (r['time'], r['[S]'], r['n_spk'], r['env']) == (0.0, 50.0, 0.0, 1.0)
    assert r.getSBML() == sbml


def test_boundary_species_concentration_is_restored_as_an_amount():
    sels = (*SELS, '[B]')
    c = _cond(selections=sels)
    cb = c.with_initial({'[B]': 3.0})
    fresh = _feedspike_model(antimony=FEEDSPIKE_B_TOY).run(c)
    for cond in (cb, c.with_parameters({'[B]': 3.0})):
        km = _feedspike_model(antimony=FEEDSPIKE_B_TOY)
        res = km.run(cond)
        assert res['[B]'][0] == 3.0 and res['env'][-1] > 1.3
        # Restored as the amount B, not as [B] at the end-of-run volume.
        assert km._te['B'] == 2.0
        assert np.array_equal(km.run(c).array, fresh.array)
    km = _feedspike_model(antimony=FEEDSPIKE_B_TOY)
    batch = nsk.run_many(km, [cb, c, c], keep_results=True)
    assert batch.results[0]['[B]'][0] == 3.0
    for res in batch.results[1:]:
        assert np.allclose(res.array, fresh.array, rtol=1e-12, atol=0)
        assert res['[B]'][0] == 2.0


def test_overrides_restore_init_entries_settings_and_amounts():
    km = _ia_model()
    r = km._te
    sbml = r.getSBML()
    with km.overrides({'init(k)': 0.7}):
        assert (r['k'], r['k2'], r['[S]'], r['k3']) == (0.7, 1.4, 14.0, 2.0)
        km.simulate(0, 1, 11)
    assert (r['k'], r['k3'], r['init(k)'], r['k2']) == (0.4, 2.0, 0.5, 1.0)
    with km.overrides({'k': 0.3, 'init(k)': 0.7}):
        assert r['k'] == 0.3 and r['[S]'] == 14.0      # init(...) goes first
    assert (r['k'], r['k3'], r['init(k)']) == (0.4, 2.0, 0.5)
    assert r.getSBML() == sbml
    km = _feedspike_model(antimony=FEEDSPIKE_B_TOY)
    with km.overrides({'[B]': 3.0}):
        km.run(_cond())                    # grows env to 1.331
    assert km._te['B'] == 2.0


def test_flux_summary_after_run_uses_the_values_the_run_used():
    toy = """
    model m()
      compartment env = 1; species S in env, P in env;
      S = 10; P = 0; k = 0.2;
      J: S => P; k*S;
    end"""
    km = nsk.KineticModel.from_antimony(toy, units=UNITS)
    sels = km.state_selections()
    km._te['k'] = 0.8                    # reference: k = 0.8 set by hand
    km.simulate(0, 5, 501, sels, reset=True)
    ref = nsk.compute_flux_summary(km, {}).cumulative_flux['J']
    km._te['k'] = 0.2
    km.simulate(0, 5, 501, sels, reset=True)
    ref_02 = nsk.compute_flux_summary(km, {}).cumulative_flux['J']
    assert not np.isclose(ref, ref_02)
    for cond in (nsk.Condition(parameters={'k': 0.8}, t_end=5, n_steps=501,
                               selections=sels),
                 nsk.Condition(initial={'k': 0.8}, t_end=5, n_steps=501,
                               selections=sels)):
        km.run(cond)
        flux = nsk.compute_flux_summary(km, {}).cumulative_flux['J']
        assert np.isclose(flux, ref, rtol=1e-9, atol=0)
        assert km._te['k'] == 0.2            # restored, also after the replay
    # A plain simulate() at k = 0.2 clears the recorded values.
    km.simulate(0, 5, 501, sels, reset=True)
    assert np.isclose(nsk.compute_flux_summary(km, {}).cumulative_flux['J'],
                      ref_02, rtol=1e-12, atol=0)


def test_run_restores_parameters_when_integration_fails():
    km = nsk.KineticModel.from_antimony(BLOW_UP_TOY, units=UNITS)
    with pytest.raises(KineticSimulationError):
        km.run(nsk.Condition(parameters={'k': 2.0}, t_end=2, n_steps=21))
    assert km._te.k == 1.0
    # The model recovers: the next run equals a fresh model's.
    cond = nsk.Condition(parameters={'k': 0.3}, t_end=2, n_steps=21)
    fresh = nsk.KineticModel.from_antimony(BLOW_UP_TOY, units=UNITS)
    out = km.run(cond)
    assert np.array_equal(out.array, fresh.run(cond).array)
    assert np.isclose(out.final['[S]'], 1 / (1 - 0.3 * 2), rtol=1e-4)


def test_run_returns_a_times_grid_exactly():
    km = _feedspike_model()
    grid = (0.0, 1.5, 3.0, 7.25, 30.0, 48.0)
    res = km.run(nsk.Condition(times=grid, selections=['[P]', 'n_spk']))
    assert tuple(res.time) == grid
    assert res.col_names == ['time', '[P]', 'n_spk']    # 'time' is prepended
    dense = km.run(_cond(n_steps=97, selections=['time', '[P]']))
    assert np.isclose(res.final['[P]'], dense.final['[P]'], rtol=1e-5)


def test_repeated_and_copied_runs_are_identical():
    km = _feedspike_model()
    # No explicit '[S]': the start must come from the baseline reset alone.
    cond = _cond(parameters={'k': 0.4}, selections=SELS)
    r1 = km.run(cond)
    # A dirty run with another k. (Before the baseline reset restored volumes
    # ahead of RoadRunner's reset, this left an ulp in [S]0 and r2 drifted.)
    km.run(cond.with_parameters(k=0.8))
    r2 = km.run(cond)
    assert np.array_equal(r1.array, r2.array)
    assert np.isclose(r1.final['[P]'], 161.9, atol=0.05)
    assert r1.final['n_spk'] == 3.0
    # A fresh model and a copy of the dirty one give the same trajectory.
    assert np.array_equal(_feedspike_model().run(cond).array, r1.array)
    assert np.allclose(km.copy().run(cond).array, r1.array, rtol=1e-12,
                       atol=0)


def test_run_reset_kwargs_and_reset_false():
    calls = []

    def custom_reset(model, **kwargs):
        calls.append(kwargs)
        model._te.reset()

    km = _feedspike_model(reset=custom_reset)
    km.run(_cond(reset_kwargs={'reset_spike_cap': True}))
    assert calls == [{'reset_spike_cap': True}]
    # The reset function gets plain containers, not the frozen ones.
    km.run(_cond(reset_kwargs={'opts': {'caps': [1, 2]}}))
    assert calls[-1] == {'opts': {'caps': [1, 2]}}
    assert type(calls[-1]['opts']) is dict
    assert type(calls[-1]['opts']['caps']) is list
    calls.pop()
    # A mode in reset_kwargs bypasses the custom function.
    km.run(_cond(reset_kwargs={'mode': 'baseline'}))
    assert len(calls) == 1
    # reset=False continues from the current state (end of the last run).
    km.run(_cond(t_end=10, n_steps=11, selections=SELS))
    end = km._te['[P]']
    res = km.run(_cond(t_end=20, n_steps=11, selections=SELS, reset=False,
                       initial={}, parameters={}).replace(t0=10))
    assert len(calls) == 2 and res['[P]'][0] == end and res.time[0] == 10.0


def test_run_does_not_touch_events_or_parameters():
    km = _feedspike_model()
    params = km._user_parameter_values()
    events, compiled = list(km.events), km._events_compiled
    km.run(_cond(parameters={'k': 0.9, 'threshold': 8}))
    assert km._user_parameter_values() == params
    assert km.events == events and km._events_compiled == compiled


# --- SimulationResult (e, f) ------------------------------------------------------
def test_simulation_result_accessors_and_interpolation():
    km = _feedspike_model()
    res = km.run(_cond(selections=SELS, label='FB'))
    assert res.label == 'FB' and res.units == UNITS
    assert res.units is not km._units
    assert res.col_names == list(SELS) and len(res) == 481
    assert res.array.shape == (481, 5)
    assert np.array_equal(res['[P]'], res.array[:, 2])
    # Accessors return new, writable arrays (pandas copy-on-write would
    # otherwise hand out read-only views).
    s_col, first = res['[S]'], res['[S]'][0]
    s_col -= 1.0
    times = res.time
    times[0] = -1.0
    block = res.array
    block[0, 1] = -1.0
    pair = res[['time', '[S]']]
    pair[0, 1] = -1.0
    assert res['[S]'][0] == first and res.time[0] == 0.0
    assert res[['time', 'env']].shape == (481, 2)
    assert '[P]' in res and '[X]' not in res
    with pytest.raises(KeyError, match='available'):
        res['[X]']
    for bad in (0, 1.5, None, {'[S]': 1}):
        with pytest.raises(TypeError, match=r'res\.df\.iloc\[i\]'):
            res[bad]
    assert res[('[S]', 'env')].shape == (481, 2)      # any sequence of names
    assert isinstance(res.final, pd.Series)
    assert res.final['time'] == 48.0 and res.final['n_spk'] == 3.0
    # at(): linear interpolation of every column, exact at grid points.
    t = 12.34
    at = res.at(t)
    for col in SELS[1:]:
        assert at[col] == np.interp(t, res.time, res[col])
    assert at['time'] == t
    assert res.at(res.time[120])['[P]'] == res['[P]'][120]
    table = res.at([0.0, t, 48.0])
    assert table.shape == (3, 5) and table['[P]'].iloc[1] == at['[P]']
    for bad in (-0.1, 48.01, np.nan, [1.0, 49.0]):
        with pytest.raises(ValueError, match='outside'):
            res.at(bad)


def test_simulation_result_validation():
    with pytest.raises(ValueError, match='time'):
        nsk.SimulationResult(df=pd.DataFrame({'[S]': [1.0, 2.0]}))
    with pytest.raises(TypeError):
        nsk.SimulationResult(df={'time': [0, 1]})
    with pytest.raises(TypeError):
        nsk.SimulationResult(df=pd.DataFrame({'time': [0.0, 1.0]}),
                             condition={'t_end': 1})


def test_simulation_result_csv_round_trip_is_exact(tmp_path):
    km = _feedspike_model()
    cond = nsk.Condition(initial={'[S]': 50}, parameters={'k': 0.1 + 0.2},
                         times=[0, 0.5, 1.25, 48], selections=SELS,
                         reset_kwargs={'preserve': ['n_max']},
                         label='run "1", µ')
    res = km.run(cond)
    path = tmp_path / 'run1.csv'
    res.to_csv(path)
    text = path.read_text(encoding='utf-8').splitlines()
    assert text[0] == '# nskinetics: "SimulationResult"'
    assert text[1] == '# units: {"time": "h", "conc": "g/L"}'
    back = nsk.SimulationResult.from_csv(path)
    pd.testing.assert_frame_equal(back.df, res.df, check_exact=True)
    assert back.units == res.units and back.label == res.label
    assert back.condition == cond

    # Awkward values and column names survive too.
    df = pd.DataFrame({'time': [0.0, 1e-300, 1.0],
                       'ec(J, S)': [0.1 + 0.2, -0.0, 1 / 3],
                       "S'": [np.nan, np.inf, 6.02214076e23]})
    raw = nsk.SimulationResult(df=df, units={'time': 'min'})
    raw.to_csv(tmp_path / 'raw.csv')
    back = nsk.SimulationResult.from_csv(tmp_path / 'raw.csv')
    pd.testing.assert_frame_equal(back.df, df, check_exact=True)
    assert back.condition is None and back.label is None
    assert np.signbit(back.df['ec(J, S)'].iloc[1])

    plain = tmp_path / 'plain.csv'
    df.to_csv(plain, index=False)
    with pytest.raises(ValueError, match='to_csv'):
        nsk.SimulationResult.from_csv(plain)


def test_plot_shares_the_plot_simulation_results_core(close_figures):
    km = _feedspike_model()
    res = km.run(_cond(selections=SELS))
    fig1, ax1 = km.plot_simulation_results(markers=[(10, 'a')], harvest=40)
    fig2, ax2 = res.plot(markers=[(10, 'a')], harvest=40)
    assert ax1.get_xlabel() == ax2.get_xlabel()
    assert ax1.get_ylabel() == ax2.get_ylabel()
    assert ax1.get_xlim() == ax2.get_xlim() and ax1.get_ylim() == ax2.get_ylim()
    assert ([l.get_label() for l in ax1.get_lines()]
            == [l.get_label() for l in ax2.get_lines()] == ['[S]', '[P]'])
    for l1, l2 in zip(ax1.get_lines(), ax2.get_lines()):
        assert np.array_equal(l1.get_xydata(), l2.get_xydata())
    _, ax3 = res.plot(variables=['env'], ylabel='Volume', ylabel_unit='L')
    assert 'Volume' in ax3.get_ylabel()
    with pytest.raises(KineticSimulationError):
        nsk.SimulationResult(df=res.df[['time', 'env']]).plot()
