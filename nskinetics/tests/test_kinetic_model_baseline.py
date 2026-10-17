# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Baseline reset modes, ``preserve``, ``simulate(times=)`` and the
NSKBatchReactor reset paths (F1 items 1-3)."""
import warnings

import numpy as np
import pytest
import tellurium as te

import nskinetics as nsk

# The FeedSpike toy of the verified reset bug (specs/repro_reset_bug.py): S is
# consumed into P and refilled to `target` (volume grows) whenever it drops
# below `threshold`, at most `n_max` times.
FEEDSPIKE_TOY = """
model fb()
  compartment env = 1; species S in env, P in env;
  S = 50; P = 0; k = 0.5; threshold = 5; target = 50; feed_conc = 500;
  n_max = 3; n_spk = 0; last_vol = 0; tot_vol = 0;
  J: S => P; k*S*env;
end
"""

SELS = ['time', '[S]', '[P]', 'env', 'n_spk']


def _feedspike_model(antimony=FEEDSPIKE_TOY, reset=None):
    km = nsk.KineticModel(te.loadAntimonyModel(antimony),
                          units={'time': 'h', 'conc': 'g/L'}, reset=reset)
    spike = nsk.FeedSpike('S', 'S <= threshold', 'target', 'feed_conc',
                          volume_var='env', max_count='n_max',
                          count_var='n_spk', last_vol_var='last_vol',
                          tot_vol_var='tot_vol')
    for event in spike.expand():
        km.add_event(event)
    km.compile_events()
    return km


def _col(arr, name):
    return arr[:, SELS.index(name)]


# --- (a) the verified bug ---------------------------------------------------
def test_repeated_reset_runs_are_identical():
    km = _feedspike_model()
    a1 = km.simulate(0, 48, 481, SELS, reset=True).copy()
    a2 = km.simulate(0, 48, 481, SELS, reset=True).copy()
    assert np.array_equal(a1, a2)
    for a in (a1, a2):
        assert np.isclose(_col(a, '[P]')[-1], 161.9, atol=0.05)
        assert _col(a, 'n_spk')[-1] == 3.0
        assert _col(a, 'env')[0] == 1.0 and _col(a, 'n_spk')[0] == 0.0


def test_roadrunner_mode_keeps_the_old_stale_start():
    km = _feedspike_model()
    km.simulate(0, 48, 481, SELS, reset=True)          # dirty the state
    km.reset(mode='roadrunner')
    a = km.simulate(0, 48, 481, SELS)
    # The old default: compartment and event-assigned counter are not restored.
    assert np.isclose(_col(a, 'env')[0], 1.331)
    assert _col(a, 'n_spk')[0] == 3.0
    assert np.isclose(_col(a, '[P]')[-1], 37.57, atol=0.05)


def test_baseline_variables_cover_compartments_and_event_targets():
    km = _feedspike_model()
    assert km.baseline_variables == ('env', 'last_vol', 'tot_vol', 'n_spk')


def test_baseline_restores_native_antimony_event_targets():
    # Events written in Antimony (not via nsk.Event) must be restored too.
    km = nsk.KineticModel(te.loadAntimonyModel("""
        model m()
          species S; S = 10; k = 0.5; n = 0;
          J: S => ; k*S;
          E1: at (S < 5): S = S + 5, n = n + 1;
        end"""), units={'time': 'h', 'conc': 'g/L'})
    assert 'n' in km.baseline_variables
    a1 = km.simulate(0, 6, 61, ['time', '[S]', 'n'], reset=True).copy()
    a2 = km.simulate(0, 6, 61, ['time', '[S]', 'n'], reset=True).copy()
    assert a1[-1, 2] == 4.0
    assert np.array_equal(a1, a2)
    # closed form: each event doubles S (5 -> 10), so S(6) = 10*e^(-3)*2^4
    assert np.isclose(a1[-1, 1], 160 * np.exp(-3), rtol=1e-4)



def test_event_target_cache_follows_assignments_added_to_an_event():
    # The cache of event targets must follow an assignment added to an
    # existing event (same event ids), or the baseline reset goes stale.
    km = nsk.KineticModel(te.loadAntimonyModel("""
        model m()
          species S; S = 10; k = 0.5; n = 0;
          J: S => ; k*S;
          E1: at (S < 5): S = S + 5;
        end"""), units={'time': 'h', 'conc': 'g/L'})
    assert 'n' not in km.baseline_variables        # fills the cache
    km._te.addEventAssignment('E1', 'n', 'n + 1', True)
    assert 'n' in km.baseline_variables
    sels = ['time', '[S]', 'n']
    a1 = km.simulate(0, 6, 7, sels, reset=True).copy()
    a2 = km.simulate(0, 6, 7, sels, reset=True).copy()
    assert a1[0, 2] == 0.0 and a1[-1, 2] == 4.0
    assert np.array_equal(a1, a2)
    assert np.allclose(km.copy().simulate(0, 6, 7, sels), a1,
                       rtol=1e-12, atol=0)


def test_flux_summary_needs_native_event_targets():
    # state_selections() includes the targets of events written in Antimony,
    # so compute_flux_summary needs them in the trajectory (a replay without
    # them would use stale values).
    km = nsk.KineticModel(te.loadAntimonyModel("""
        model m()
          compartment env = 1; species S in env, P in env;
          S = 10; P = 0; k = 0.5; n = 0;
          J: S => P; k*S;
          E1: at (S < 5): S = S + 5, n = n + 1;
        end"""), units={'time': 'h', 'conc': 'g/L'})
    assert 'n' in km.state_selections()
    km.simulate(0, 6, 61, ['time', 'env', '[S]', '[P]'], reset=True)
    with pytest.raises(nsk.exceptions.KineticSimulationError,
                       match=r"missing state columns.*\['n'\]"):
        nsk.compute_flux_summary(km, {})
    km.simulate(0, 6, 601, km.state_selections(), reset=True)
    summary = nsk.compute_flux_summary(km, {})
    # P only ever grows through J, so J's cumulative flux is P(6) (up to the
    # trapezoidal rule across the event jumps).
    assert np.isclose(summary.cumulative_mass['J'],
                      km.results_df['[P]'].iloc[-1], rtol=2e-3)


def test_spike_reduce_retry_reads_and_writes_the_cap_by_item():
    # libroadrunner removes a model's r.<id> properties whenever another
    # model is loaded; the retry loop must keep working on such a model.
    r = te.loadAntimonyModel("model a(); n_max = 3; x = 0; end")
    te.loadAntimonyModel("model b(); y = 1; end")
    assert not hasattr(type(r), 'n_max')
    seen = []
    retry = nsk.units.SpikeReduceRetry(
        max_count_var='n_max', stop_when=lambda m: m['n_max'] <= 1)
    retry.run(r, lambda: seen.append(r['n_max']))
    assert seen == [3.0, 2.0] and r['n_max'] == 1.0
    assert 'n_max' not in vars(r)          # no stray instance attribute


def test_baseline_includes_and_restores_rate_rule_variables():
    km = nsk.KineticModel(te.loadAntimonyModel("""
        model m()
          species S; S = 1; z = 2; z' = 0.5;
          J: S => ; 0.1*S;
        end"""), units={'time': 'h', 'conc': 'g/L'})
    # RoadRunner.reset() already restores rate-rule variables; this pins
    # that the baseline path lists z (membership: RoadRunner names the
    # auto-created compartment) and keeps that behaviour.
    assert 'z' in km.baseline_variables
    km.simulate(0, 4, 5, ['time', 'z'])
    assert km._te['z'] == 4.0
    km.reset()
    assert km._te['z'] == 2.0


def test_baseline_reset_restores_volumes_before_roadrunner_reset():
    # RoadRunner's reset() converts species initial values to amounts with
    # the *current* compartment volume. Under a post-run volume this doubled
    # an initial-assigned [S]0 (40 instead of 20) and left an ulp in plain
    # initial concentrations, so runs depended on the previous run.
    km = nsk.KineticModel(te.loadAntimonyModel("""
        model m()
          compartment env = 1; species S in env;
          n = 1; S = 20*n; k = 0.5;
          J: S => ; k*S;
          E1: at (time > 1): n = n + 1, env = 2;
        end"""), units={'time': 'h', 'conc': 'g/L'})
    fresh = km.simulate(0, 3, 31, ['time', '[S]', 'env', 'n']).copy()
    km.reset()
    assert km._te['[S]'] == 20.0 and km._te.env == 1.0 and km._te.n == 1.0
    again = km.simulate(0, 3, 31, ['time', '[S]', 'env', 'n'], reset=True)
    assert np.array_equal(again, fresh)
    # The FeedSpike toy: a run after a dirty run with another k is
    # bit-identical to the same run on a freshly loaded model.
    dirty = _feedspike_model()
    dirty._te.k = 0.8
    dirty.simulate(0, 48, 481, SELS, reset=True)
    dirty._te.k = 0.4
    ref = _feedspike_model()
    ref._te.k = 0.4
    assert np.array_equal(dirty.simulate(0, 48, 481, SELS, reset=True),
                          ref.simulate(0, 48, 481, SELS, reset=True))


# --- (b) user edits persist; per-run initial values do not ------------------
def test_user_parameter_edit_persists_through_baseline_reset():
    km = _feedspike_model()
    km._te.k = 0.4
    km.simulate(0, 48, 481, SELS, reset=True)
    km.reset()
    assert km._te.k == 0.4
    assert km._te['init(k)'] == 0.5        # origin value is untouched
    a = km.simulate(0, 48, 481, SELS, reset=True)
    ref = _feedspike_model()
    ref._te.k = 0.4
    assert np.array_equal(a, ref.simulate(0, 48, 481, SELS, reset=True))


def test_initial_value_set_after_reset_applies_to_that_run_only():
    km = _feedspike_model()
    base = km.simulate(0, 48, 481, SELS, reset=True).copy()
    km.reset()
    km._te.env = 2.0                       # per-run initial volume
    a = km.simulate(0, 48, 481, SELS)
    assert _col(a, 'env')[0] == 2.0
    assert not np.allclose(a, base)
    b = km.simulate(0, 48, 481, SELS, reset=True)
    assert _col(b, 'env')[0] == 1.0
    assert np.array_equal(b, base)


# --- (c) preserve and dispatch -----------------------------------------------
def test_preserve_keeps_a_modified_cap():
    km = _feedspike_model()
    km._te.n_max = 1
    km.reset(preserve=('n_max',))
    assert km._te.n_max == 1
    a = km.simulate(0, 48, 481, SELS)
    assert _col(a, 'n_spk')[-1] == 1.0


def test_preserve_writes_back_a_baseline_variable():
    # env is restored by a baseline reset unless preserved.
    km = _feedspike_model()
    km.simulate(0, 48, 481, SELS, reset=True)
    env_end = km._te.env
    km.reset(preserve='env')               # a bare string is one name
    assert km._te.env == env_end
    assert km._te.n_spk == 0.0
    km.reset()
    assert km._te.env == 1.0


def test_preserve_protects_an_event_assigned_cap():
    # A cap that some event assigns is a baseline variable: a plain baseline
    # reset restores it, preserve keeps it (the NSKBatchReactor retry case).
    km = _feedspike_model(FEEDSPIKE_TOY.replace(
        'end', 'hold: at (time > 1e9): n_max = n_max;\nend'))
    assert 'n_max' in km.baseline_variables
    km._te.n_max = 2
    km.reset()
    assert km._te.n_max == 3.0
    km._te.n_max = 2
    km.reset(preserve=('n_max',))
    assert km._te.n_max == 2.0
    assert _col(km.simulate(0, 48, 481, SELS), 'n_spk')[-1] == 2.0


def test_custom_reset_dispatch_and_preserve_not_forwarded():
    seen = []

    def my_reset(model, **kwargs):
        seen.append(kwargs)
        model._te.reset()
        model._te.n_max = 7

    km = _feedspike_model(reset=my_reset)
    km._te.n_max = 2
    km.reset(preserve=('n_max',), reset_spike_cap=False)
    assert seen == [{'reset_spike_cap': False}]
    assert km._te.n_max == 2.0             # preserve wins over reset_func
    km.reset()
    assert km._te.n_max == 7.0
    km.reset(mode='baseline')              # explicit modes bypass reset_func
    assert len(seen) == 2


def test_reset_argument_errors_touch_nothing():
    km = _feedspike_model()
    km.simulate(0, 48, 481, SELS, reset=True)
    before = {s: km._te[s] for s in km.state_selections()}
    with pytest.raises(ValueError):
        km.reset(mode='bogus')
    with pytest.raises(ValueError):
        km.reset(mode='custom')            # no reset_func
    with pytest.raises(ValueError):
        km.reset(preserve=('not_an_id',))
    # readable but not writable back: a reaction rate and a rate of change
    with pytest.raises(ValueError, match='read-only'):
        km.reset(preserve=('J', "S'"))
    assert {s: km._te[s] for s in km.state_selections()} == before


class _StubTe:
    def __init__(self):
        self.reset_calls = 0

    def reset(self):
        self.reset_calls += 1


def test_baseline_variables_survive_a_missing_event_target_cache():
    # e.g. a model unpickled from before the cache existed; the variables must
    # not silently collapse to () (which degrades to the stale-start reset).
    km = nsk.KineticModel(te.loadAntimonyModel("""
        model m()
          species S; S = 10; k = 0.5; n = 0;
          J: S => ; k*S;
          E1: at (S < 5): S = S + 5, n = n + 1;
        end"""))
    expected = km.baseline_variables
    assert 'n' in expected
    del km._event_targets_cache
    assert km.baseline_variables == expected


def test_stub_te_baseline_falls_back_to_plain_reset():
    stub = _StubTe()
    km = nsk.KineticModel(stub)
    assert km.baseline_variables == ()
    km.reset(mode='baseline')
    km.reset(mode='roadrunner')
    assert stub.reset_calls == 2


# --- (3) simulate(times=) ----------------------------------------------------
def test_simulate_times_grid_is_returned_exactly():
    km = _feedspike_model()
    grid = [0.0, 0.3, 1.7, 5.25, 12.0, 48.0]
    a = km.simulate(selections=SELS, reset=True, times=grid)
    assert a[:, 0].tolist() == grid
    assert km.results_col_names == SELS
    assert km.results_df['time'].tolist() == grid
    # events still fire between output points; endpoint matches a dense run
    dense = km.simulate(0, 48, 481, SELS, reset=True)
    assert np.allclose(a[-1], dense[-1], rtol=1e-6)


@pytest.mark.parametrize('kwargs', [
    {'times': [0.0, 2.0, 1.0]},
    {'times': [0.0, 1.0, 1.0]},
    {'times': [1.0]},
    {'times': [0.0, np.nan]},
    {'times': [0.0, 1.0], 't_end': 1.0},
    {'times': [1.0, 2.0], 't0': 0.5},
    # an explicit t0 = 0 is not silently moved to times[0]
    {'times': [1.0, 2.0], 't0': 0.0},
    {'times': 'abc'},
    {'times': {'a': 1.0}},
])
def test_simulate_invalid_times_raise_before_touching_state(kwargs):
    km = _feedspike_model()
    km.simulate(0, 48, 481, SELS, reset=True)
    before = {s: km._te[s] for s in km.state_selections()}
    results = km.results_array
    # Not a sequence of numbers: TypeError (as for Condition); else ValueError.
    error = TypeError if isinstance(kwargs['times'], (str, dict)) else ValueError
    with pytest.raises(error):
        km.simulate(selections=SELS, reset=True, **kwargs)
    assert {s: km._te[s] for s in km.state_selections()} == before
    assert km.results_array is results


def test_simulate_t0_defaults_follow_the_grid():
    km = _feedspike_model()
    a = km.simulate(selections=SELS, reset=True, times=[1.0, 2.0, 4.0])
    b = km.simulate(t0=1.0, selections=SELS, reset=True, times=[1.0, 2.0, 4.0])
    assert a[0, 0] == 1.0 and np.array_equal(a, b)
    c = km.simulate(None, 4.0, 5, SELS, reset=True)       # t0=None -> 0.0
    assert c[0, 0] == 0.0


def test_simulate_without_t_end_or_times_still_raises():
    km = _feedspike_model()
    with pytest.raises(nsk.KineticSimulationError):
        km.simulate()


# --- (l) NSKBatchReactor, slow -------------------------------------------------
@pytest.fixture
def biosteam_sandbox(request):
    """Run a biosteam test on its own flowsheet and restore the previous
    thermo and flowsheet afterwards (unit IDs never leak into other tests)."""
    import biosteam as bst
    previous_thermo = getattr(bst.settings, '_thermo', None)
    previous_flowsheet = bst.main_flowsheet.get_flowsheet()
    bst.main_flowsheet.set_flowsheet(request.node.name)
    try:
        yield bst
    finally:
        bst.main_flowsheet.set_flowsheet(previous_flowsheet)
        if previous_thermo is None:
            try:
                del bst.settings._thermo
            except AttributeError:
                pass
        else:
            bst.settings._thermo = previous_thermo


def _effluent_and_results(reactor):
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message=r'.*added unit results.*',
                                category=RuntimeWarning)
        reactor.simulate()
    return reactor.outs[1].mol.copy(), np.array(reactor.nsk_results)


@pytest.mark.slow
def test_demo_reactor_simulated_twice_gives_the_same_effluent(
        biosteam_sandbox):
    from nskinetics.models.toys.generic_batch_reactor_demo import (
        build_demo_reactor)
    reactor, te_r = build_demo_reactor()
    assert te_r.reset_func is None
    mol1, res1 = _effluent_and_results(reactor)
    mol2, res2 = _effluent_and_results(reactor)
    assert np.array_equal(res1, res2)
    assert np.array_equal(mol1, mol2)
    assert mol1.sum() > 0


_CAPPED_REACTOR_MODEL = """
model capped_batch()
  compartment env; species S in env, P in env;
  S = 50; P = 0; env = 1; k = 0.5;
  threshold = 5; target = 50; feed_conc = 500;
  n_spk = 0; n_max = 3; last_vol = 0; tot_vol = 0;
  n_spk has dimensionless; n_max has dimensionless;
  J: S => P; k*S*env;
  curr_env := env;
  hold_cap: at (time > 1e9): n_max = n_max;
end
"""


@pytest.mark.slow
def test_spike_retry_reactor_without_reset_func_keeps_its_cap(
        biosteam_sandbox):
    # `hold_cap` never fires but makes the cap event-assigned, i.e. a baseline
    # variable: without preserve, every baseline reset would restore it to 3.
    bst = biosteam_sandbox
    r = te.loadAntimonyModel(_CAPPED_REACTOR_MODEL)
    km = nsk.KineticModel(r, units={'time': 'h', 'conc': 'g/L'})
    spike = nsk.FeedSpike(
        species='S', when='S <= threshold', target='target',
        feed_conc='feed_conc', volume_var='env', max_count='n_max',
        count_var='n_spk', last_vol_var='last_vol', tot_vol_var='tot_vol',
        priority=5, name='S_spike_capped')
    for event in spike.expand():
        km.add_event(event)
    km.compile_events()
    assert km.reset_func is None and 'n_max' in km.baseline_variables

    bst.settings.set_thermo(['Water', 'Glucose', 'Ethanol'], cache=True)
    feed = bst.Stream('feed_capped', Water=1000., Glucose=50., units='kg/hr')
    seed = bst.Stream('seed_capped', Water=10., units='kg/hr')
    spike_feed = bst.Stream('spike_feed_capped', Glucose=500., units='kg/hr')
    checks = []

    def stop_when(rr):
        # Bounded so that a reset undoing the cap (which would make the
        # reduce-then-retry loop spin forever) fails the test, not hangs it.
        checks.append(rr.n_max)
        return rr.n_max <= 1 or len(checks) > 10

    retry = nsk.units.SpikeReduceRetry(
        max_count_var='n_max', count_var='n_spk', stop_when=stop_when)
    reactor = nsk.units.NSKBatchReactor(
        'R_capped', ins=(feed, seed, spike_feed), outs=('vent_capped',
                                                       'effluent_capped'),
        nsk_kinetic_model=km, map_species_to_chemicals={'S': 'Glucose',
                                                        'P': 'Ethanol'},
        track_vars=['[S]', '[P]'], tau=48., tau_max=120.,
        volume_var='curr_env', feed_volume_added_var='tot_vol',
        spike_retry=retry, spike_feed_index=2, N=2)

    def final_count():
        return reactor.nsk_results_dict['n_spk'][-1]

    # Retry path: the full run fires 3 spikes, the cap is lowered to 2 and the
    # retry run must honour it (a baseline reset would restore n_max = 3).
    reactor.simulate_kinetics(feed, 48.)
    assert checks == [2.0, 1.0]
    assert final_count() == 2.0
    assert km._te.n_max == 1.0
    # Frozen-cap (reconciliation) path.
    reactor._frozen_spike_count = 1
    try:
        reactor.simulate_kinetics(feed, 48.)
    finally:
        reactor._frozen_spike_count = None
    assert final_count() == 1.0
    assert km._te.n_max == 1.0
