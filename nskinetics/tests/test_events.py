# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
import nskinetics as nsk


def test_exception_hierarchy():
    from nskinetics.exceptions import (
        NSKError, KineticSimulationError, MassBalanceError, EventCompilationError,
    )
    for exc in (KineticSimulationError, MassBalanceError, EventCompilationError):
        assert issubclass(exc, NSKError)
    assert issubclass(NSKError, Exception)
    # re-exported at top level
    assert nsk.NSKError is NSKError
    assert nsk.EventCompilationError is EventCompilationError


import tellurium as te
from nskinetics.exceptions import KineticSimulationError

_SMALL_MODEL = """
model small()
  compartment env; species s in env; s = 10; env = 2; p = 3;
  s' = 0;
end
"""


def _make_trs(units=None):
    r = te.loadAntimonyModel(_SMALL_MODEL)
    return nsk.TelluriumReactionSystem(r, units=units)


def test_get_set_value_and_factors():
    trs = _make_trs(units={'time': 'h', 'conc': 'g/L'})
    assert trs.get_value('p') == 3.0
    trs.set_value('p', 9.0)
    assert trs.get_value('p') == 9.0
    # concentration selection round-trips through the compartment volume
    # NOTE: Antimony's bare `s = 10` sets the *initial concentration* (confirmed
    # via generated SBML: <species ... initialConcentration="10" .../>), not the
    # amount, so [s] == 10.0 directly and the derived amount is 10 * 2 = 20.
    assert trs.get_value('[s]') == 10.0            # concentration set directly
    trs.set_value('[s]', 4.0)
    assert trs.get_value('[s]') == 4.0
    assert trs.time_factor == 1.0
    assert trs.material_indexer == 'imass'
    # bracket-stripped concentration-from-stream helper
    val = trs.set_concentration_from_stream('[s]', amount=8.0, volume=2.0)
    assert val == 4.0
    assert trs.get_value('s') == 4.0              # stripped -> writes amount var 's'


def test_validate_units_rejects_unknown():
    trs = _make_trs(units={'time': 'fortnights', 'conc': 'g/L'})
    import pytest
    with pytest.raises(KineticSimulationError):
        trs.validate_units()


import numpy as np

_DECAY_MODEL = """
model decay()
  species s; s = 100; k = 1; flag = 1;
  s' = -k*flag*s;
end
"""


def test_event_fires():
    r = te.loadAntimonyModel(_DECAY_MODEL)
    ev = nsk.Event(when='time >= 5', do={'flag': '0'}, name='stop_decay')
    ev.compile(r)                      # single event -> regenerate now
    r.reset()
    res = np.array(r.simulate(0, 10, 11, ['time', 's', 'flag']))
    # flag is 1 up to t<5, then 0; decay halts after the event fires
    assert res[4, 2] == 1.0            # t=4 -> flag still 1
    assert res[5, 2] == 0.0            # t=5 -> event fired, flag 0
    assert np.isclose(res[6, 1], res[5, 1], rtol=1e-6)   # s frozen after event


def test_event_autogenerates_name():
    ev = nsk.Event(when='time >= 1', do={'flag': '0'})
    assert isinstance(ev.name, str) and ev.name


def test_compile_events_batches_and_fires():
    r = te.loadAntimonyModel(_DECAY_MODEL)
    trs = nsk.TelluriumReactionSystem(r, units={'time': 'h', 'conc': 'g/L'})
    trs.add_event(nsk.Event(when='time >= 5', do={'flag': '0'}, name='ev'))
    trs.compile_events()
    trs.reset()
    res = np.array(r.simulate(0, 10, 11, ['time', 's', 'flag']))
    assert res[5, 2] == 0.0
    # double-compile is a guarded no-op (no duplicate-event error)
    trs.compile_events()


def test_remove_and_clear_events():
    r = te.loadAntimonyModel(_DECAY_MODEL)
    trs = nsk.TelluriumReactionSystem(r)
    trs.add_event(nsk.Event(when='time >= 5', do={'flag': '0'}, name='a'))
    trs.add_event(nsk.Event(when='time >= 6', do={'flag': '0'}, name='b'))
    trs.remove_event('a')
    assert [e.name for e in trs.events] == ['b']
    trs.clear_events()
    assert trs.events == []


def test_feedspike_expansion_shape():
    fs = nsk.FeedSpike(
        species='s_glu', when='s_glu <= threshold_conc_glu_spike',
        target='target_conc_glu_spike', feed_conc='conc_glu_feed_spike',
        volume_var='env', max_count='max_n_glu_spikes', count_var='n_glu_spikes',
        last_vol_var='last_vol_glu_feed_added', tot_vol_var='tot_vol_glu_feed_added',
        delay='glucose_feed_spikeDelay', priority=5, name='glucose_feed_spike',
    )
    events = fs.expand()
    assert [e.name for e in events] == [
        'glucose_feed_spike_a', 'glucose_feed_spike_b', 'glucose_feed_spike_c']
    a, b, c = events
    # descending priorities enforce a -> b -> c ordering
    assert (a.priority, b.priority, c.priority) == (5, 4, 3)
    assert a.use_trigger_time_values is True
    assert b.use_trigger_time_values is False and c.use_trigger_time_values is False
    # max_count folds into every trigger
    for e in events:
        assert '(n_glu_spikes < max_n_glu_spikes)' in e.when
    # a computes last_vol; b updates volume + tot_vol; c sets target + counter
    assert a.do == {'last_vol_glu_feed_added':
                    'env*(target_conc_glu_spike - s_glu)/(conc_glu_feed_spike - target_conc_glu_spike)'}
    assert b.do == {'env': 'env + last_vol_glu_feed_added',
                    'tot_vol_glu_feed_added': 'tot_vol_glu_feed_added + last_vol_glu_feed_added'}
    assert c.do == {'s_glu': 'target_conc_glu_spike',
                    'n_glu_spikes': 'n_glu_spikes + 1'}


def test_feedspike_fires_and_caps_at_max_count():
    # Model where s_glu decays; each spike resets s_glu to target and bumps count.
    model = """
    model spiker()
      compartment env; species s_glu in env;
      s_glu = 100; env = 1; k = 1;
      threshold = 10; target = 100; feed_conc = 600;
      n_spk = 0; n_max = 2;
      last_vol = 0; tot_vol = 0; dly = 0;
      n_spk has dimensionless; n_max has dimensionless;
      s_glu' = -k*s_glu;
    end
    """
    r = te.loadAntimonyModel(model)
    trs = nsk.TelluriumReactionSystem(r, units={'time': 'h', 'conc': 'g/L'})
    fs = nsk.FeedSpike(
        species='s_glu', when='s_glu <= threshold',
        target='target', feed_conc='feed_conc', volume_var='env',
        max_count='n_max', count_var='n_spk',
        last_vol_var='last_vol', tot_vol_var='tot_vol',
        delay='dly', priority=5, name='spk',
    )
    for e in fs.expand():
        trs.add_event(e)
    trs.compile_events()
    trs.reset()
    res = np.array(r.simulate(0, 40, 401, ['time', 's_glu', 'n_spk']))
    # counter never exceeds the cap
    assert res[:, 2].max() == 2.0
    # after the cap, s_glu is allowed to decay below threshold and stay there
    assert res[-1, 1] < 10.0


import os

_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
_REF_ANTIMONY = os.path.join(_DATA_DIR, 'ibo_antimony_with_events_reference.txt')

# Event-block line prefixes to strip when building the "base" (no-events) model.
_EVENT_NAMES = ('glucose_feed_spike_a', 'glucose_feed_spike_b',
                'glucose_feed_spike_c', 'stage_1_complete_max_time',
                'stage_1_complete_x_target')


def _strip_event_blocks(antimony_text):
    """Remove the five hand-written event definition lines, keeping everything
    else (species, rate laws, parameters, assignment rules)."""
    kept = []
    for line in antimony_text.splitlines():
        stripped = line.strip()
        if any(stripped.startswith(name + ':') for name in _EVENT_NAMES):
            continue
        kept.append(line)
    return '\n'.join(kept)


def _api_events():
    """The isobutanol events, declared via the Python API (mirrors Task 11)."""
    spike = nsk.FeedSpike(
        species='s_glu', when='s_glu <= threshold_conc_glu_spike',
        target='target_conc_glu_spike', feed_conc='conc_glu_feed_spike',
        volume_var='env', max_count='max_n_glu_spikes', count_var='n_glu_spikes',
        last_vol_var='last_vol_glu_feed_added', tot_vol_var='tot_vol_glu_feed_added',
        delay='glucose_feed_spikeDelay', priority=5, name='glucose_feed_spike')
    stage_time = nsk.Event(when='time >= stage_1_max_time', do={'is_aerobic': '0'},
                           name='stage_1_complete_max_time')
    stage_x = nsk.Event(when='x >= stage_1_max_x', do={'is_aerobic': '0'},
                        name='stage_1_complete_x_target')
    return spike, stage_time, stage_x


def _setup_ibo_run(r):
    """Set integrator tolerances and initial conditions matching the example."""
    integ = r.getIntegrator()
    integ.absolute_tolerance = 1e-10
    integ.relative_tolerance = 1e-9
    r.reset()
    r.n_glu_spikes = 0
    r.last_vol_glu_feed_added = 0.
    r.tot_vol_glu_feed_added = 0.
    r.env = 1.
    r.is_aerobic = 1
    r.max_n_glu_spikes = 20
    r.s_glu = 100
    r.x = 1


_COLS = ['time', 's_glu', 's_EtOH', 's_IBO', 'x', 'n_glu_spikes',
         'curr_tot_vol_glu_feed_added']


_RATE_RULE_MODEL = """
model rate_rule_bug()
  compartment env; species s in env; s = 10; env = 1;
  w = 0.5; w has dimensionless;
  w' = 0;
  s' = 0;
  flag = 1;
end
"""


def test_compile_events_preserves_rate_rule_initial_values():
    # Regression test for a roadrunner bug: regenerateModel() (triggered by
    # compiling an event) leaves the model's stored init(X) bookkeeping stale
    # for any rate-rule-governed variable with an explicit initial value; a
    # plain reset() afterward then reads back a corrupted value instead of
    # the SBML-defined initial value. compile_events() guards against this by
    # calling resetToOrigin() right after regenerateModel(). `w` here is a
    # rate-rule-governed variable (w' = 0) with an explicit initial value
    # (w = 0.5), i.e. exactly the kind of variable the bug corrupts.
    r = te.loadAntimonyModel(_RATE_RULE_MODEL)
    trs = nsk.TelluriumReactionSystem(r, units={'time': 'h', 'conc': 'g/L'})
    trs.add_event(nsk.Event(when='time >= 5', do={'flag': '0'}, name='ev'))
    trs.compile_events()
    trs.reset()
    assert r['w'] == 0.5


def test_migration_equivalence_ibo_events():
    with open(_REF_ANTIMONY, encoding='utf-8') as f:
        ref_text = f.read()

    # Reference: hand-written Antimony events.
    r_ref = te.loadAntimonyModel(ref_text)
    _setup_ibo_run(r_ref)
    out_ref = np.array(r_ref.simulate(0, 200, 2001, _COLS))

    # Migrated: base model (events stripped) + events injected via the API.
    r_api = te.loadAntimonyModel(_strip_event_blocks(ref_text))
    trs = nsk.TelluriumReactionSystem(r_api, units={'time': 'h', 'conc': 'g/L'})
    spike, stage_time, stage_x = _api_events()
    for e in spike.expand():
        trs.add_event(e)
    trs.add_event(stage_time)
    trs.add_event(stage_x)
    trs.compile_events()
    _setup_ibo_run(r_api)
    out_api = np.array(r_api.simulate(0, 200, 2001, _COLS))

    assert np.allclose(out_ref, out_api, rtol=1e-3, atol=1e-6)


def test_select_tau_index_policies():
    from nskinetics.units.batch_reactor import select_tau_index
    cols = ['time', 'val']
    results = np.array([[0.0, 1.0], [1.0, 5.0], [2.0, 5.0], [3.0, 2.0]])
    # None -> nearest time index
    idx, ok = select_tau_index(results, cols, tau=2.4, policy=None)
    assert (idx, ok) == (2, True)
    # ('max', var) -> first index achieving the max
    idx, ok = select_tau_index(results, cols, tau=None, policy=('max', 'val'))
    assert (idx, ok) == (1, True)
    # ('min', var) -> first index achieving the min
    idx, ok = select_tau_index(results, cols, tau=None, policy=('min', 'val'))
    assert (idx, ok) == (0, True)
    # ('equals', var, value) -> first index equal (rounded) to value
    idx, ok = select_tau_index(results, cols, tau=None, policy=('equals', 'val', 2.0))
    assert (idx, ok) == (3, True)
    # ('equals', ...) with no match -> success False, index -1
    idx, ok = select_tau_index(results, cols, tau=None, policy=('equals', 'val', 99.0))
    assert ok is False and idx == -1
