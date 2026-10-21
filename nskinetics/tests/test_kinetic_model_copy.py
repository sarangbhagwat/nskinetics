# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""KineticModel.copy(), the picklable _ModelBlueprint, from_antimony /
from_sbml and overrides() (F1 items 7-9)."""
import pickle

import numpy as np
import pytest
import tellurium as te

import nskinetics as nsk
from nskinetics.engine.kinetic_model import _ModelBlueprint

FEEDSPIKE_TOY = """
model fb()
  compartment env = 1; species S in env, P in env;
  S = 50; P = 0; k = 0.5; threshold = 5; target = 50; feed_conc = 500;
  n_max = 3; n_spk = 0; last_vol = 0; tot_vol = 0;
  J: S => P; k*S*env;
  curr_env := env;
end
"""

SELS = ['time', '[S]', '[P]', 'env', 'n_spk', 'tot_vol']


def module_level_reset(model, reset_spike_cap=True, **kwargs):
    """Picklable custom reset (module level), mirroring the shipped model's."""
    model._te.reset()
    for var in ('n_spk', 'last_vol', 'tot_vol'):
        model._te[var] = 0.0
    model._te.env = 1.0
    if reset_spike_cap:
        model._te.n_max = model.default_n_max


def _feedspike_model(compile_events=True, wrapper=True, **kwargs):
    if wrapper:     # a KineticModel built around a live RoadRunner
        km = nsk.KineticModel(te.loadAntimonyModel(FEEDSPIKE_TOY),
                              units={'time': 'h', 'conc': 'g/L'}, **kwargs)
    else:
        km = nsk.KineticModel.from_antimony(
            FEEDSPIKE_TOY, units={'time': 'h', 'conc': 'g/L'}, **kwargs)
    spike = nsk.FeedSpike('S', 'S <= threshold', 'target', 'feed_conc',
                          volume_var='env', max_count='n_max',
                          count_var='n_spk', last_vol_var='last_vol',
                          tot_vol_var='tot_vol')
    for event in spike.expand():
        km.add_event(event)
    if compile_events:
        km.compile_events()
    return km


def _state(km):
    return {s: km._te[s] for s in km.state_selections()}


# --- (g) copy equivalence -----------------------------------------------------
@pytest.mark.parametrize('wrapper', [True, False])
def test_copy_after_dirty_run_reproduces_the_run(wrapper):
    km = _feedspike_model(wrapper=wrapper)
    km._te.k = 0.4
    ref = km.simulate(0, 48, 481, SELS, reset=True).copy()
    km.simulate(0, 48, 481, SELS)            # leave the state dirty
    clone = km.copy()
    out = clone.simulate(0, 48, 481, SELS, reset=True)
    assert np.allclose(out, ref, rtol=1e-12, atol=0)
    # a fresh copy starts at the baseline without needing a reset
    fresh = km.copy().simulate(0, 48, 481, SELS)
    assert np.allclose(fresh, ref, rtol=1e-12, atol=0)


def test_copy_carries_settings_and_attributes():
    km = _feedspike_model(reset=module_level_reset)
    km.default_n_max = 2
    km.notes = ['calibrated']
    km._te.k = 0.4
    integ = km._te.getIntegrator()
    integ.absolute_tolerance = 1e-11
    integ.relative_tolerance = 1e-9
    integ.maximum_num_steps = 12345
    km.simulate(0, 48, 481, SELS, reset=True)

    clone = km.copy()
    assert clone._te is not km._te
    assert type(clone) is nsk.KineticModel
    assert clone._te.k == 0.4
    assert clone._units == km._units and clone._units is not km._units
    assert clone.reset_func is module_level_reset
    cinteg = clone._te.getIntegrator()
    assert cinteg.getName() == integ.getName()
    assert cinteg.getValue('absolute_tolerance') == 1e-11
    assert cinteg.getValue('relative_tolerance') == 1e-9
    assert cinteg.getValue('maximum_num_steps') == 12345
    assert clone.default_n_max == 2
    assert clone.notes == ['calibrated'] and clone.notes is not km.notes
    assert [e.name for e in clone.events] == [e.name for e in km.events]
    assert clone.events is not km.events
    assert clone._events_compiled is True
    assert clone.results_df is None and clone.results_array is None
    # the custom reset (which reads the copied attribute) behaves the same
    ref = km.simulate(0, 48, 481, SELS, reset=True)
    out = clone.simulate(0, 48, 481, SELS, reset=True)
    assert np.allclose(out, ref, rtol=1e-12, atol=0)
    assert out[-1, SELS.index('n_spk')] == 2.0


def test_mutating_the_copy_leaves_the_original_unchanged():
    km = _feedspike_model()
    km.flags = {'a': 1}
    km.simulate(0, 48, 481, SELS, reset=True)
    params = km.parameters['value'].to_dict()
    state = _state(km)
    n_events = len(km.events)
    units = dict(km._units)

    clone = km.copy()
    clone._te.k = 9.0
    clone._te.n_max = 0
    clone._units['conc'] = 'M'
    clone.flags['a'] = 2
    clone.add_event(nsk.Event(when='time >= 1', do={'k': '1'}, name='late'))
    clone.simulate(0, 10, 11, SELS, reset=True)

    assert km.parameters['value'].to_dict() == params
    assert _state(km) == state
    assert len(km.events) == n_events and km._events_compiled
    assert km._units == units and km.flags == {'a': 1}


def test_copy_keeps_initial_assignments_live():
    # getCurrentSBML() drops initial assignments; the copy must not.
    km = nsk.KineticModel(te.loadAntimonyModel("""
        model ia()
          compartment env = 1; species S in env, P in env;
          S = 20*k; P = 0; k = 0.5;
          J: S => P; k*S*env;
        end"""), units={'time': 'h', 'conc': 'g/L'})
    km._te.k = 0.4
    km.simulate(0, 5, 6, ['time', '[S]'], reset=True)
    clone = km.copy()
    assert clone._te['[S]'] == 8.0          # baseline with the copied k
    for k in (0.4, 0.2):
        km._te.k = clone._te.k = k
        a = km.simulate(0, 5, 6, ['time', '[S]'], reset=True)
        b = clone.simulate(0, 5, 6, ['time', '[S]'], reset=True)
        assert a[0, 1] == b[0, 1] == 20 * k
        assert np.allclose(a, b, rtol=1e-12, atol=0)


def test_copy_carries_boundary_species_edits():
    # RoadRunner.reset() never restores a boundary species, so an edit
    # persists across resets of the original; the copy must carry it.
    km = nsk.KineticModel.from_antimony("""
        model b()
          compartment c = 1; species S in c, $X in c;
          S = 0; X = 2; k = 0.3;
          J: => S; k*X;
        end""")
    km._te.X = 5.0
    sels = ['time', '[S]', '[X]']
    ref = km.simulate(0, 10, 11, sels, reset=True).copy()
    km.simulate(0, 10, 11, sels, reset=True)
    assert ref[0, 2] == 5.0 and np.isclose(ref[-1, 1], 0.3 * 5.0 * 10)
    clone = km.copy()
    assert clone._te.X == 5.0
    out = clone.simulate(0, 10, 11, sels, reset=True)
    assert np.allclose(out, ref, rtol=1e-12, atol=0)


def test_origin_sbml_is_the_live_get_sbml():
    # The origin SBML is read from the RoadRunner on every use (never a
    # cached copy that could go stale), and a run never changes it.
    for wrapper in (False, True):
        km = _feedspike_model(wrapper=wrapper, compile_events=False)
        at_load = km._get_origin_sbml()
        assert at_load == km._te.getSBML() and '<event' not in at_load
        km.compile_events()
        compiled = km._get_origin_sbml()
        assert compiled == km._te.getSBML()
        assert compiled.count('<event ') == 3
        km.simulate(0, 48, 481, SELS, reset=True)
        km._te.k = 0.3                       # a user edit is not baked in
        assert km._get_origin_sbml() == compiled
        # getCurrentSBML() would bake the run state into the initial values
        assert km._te.getCurrentSBML() != compiled
        assert km._origin_sbml is None       # only a fallback for stubs


def test_copy_of_wrapper_without_compiled_events_keeps_them_pending():
    km = _feedspike_model(compile_events=False)
    assert km._origin_sbml is None
    before = _state(km)
    clone = km.copy()
    assert _state(km) == before               # reading the origin is pure
    assert clone._events_compiled is False and len(clone.events) == 3
    clone.compile_events()
    km.compile_events()
    assert np.allclose(clone.simulate(0, 48, 481, SELS, reset=True),
                       km.simulate(0, 48, 481, SELS, reset=True),
                       rtol=1e-12, atol=0)


COND = nsk.Condition(t_end=48, n_steps=481, selections=SELS)


def _spike():
    return nsk.FeedSpike('S', 'S <= threshold', 'target', 'feed_conc',
                         volume_var='env', max_count='n_max',
                         count_var='n_spk', last_vol_var='last_vol',
                         tot_vol_var='tot_vol')


def test_copy_carries_events_compiled_directly_on_the_roadrunner():
    # FeedSpike.compile / Event.compile edit the RoadRunner without the
    # wrapper knowing; copies (and so run_many) must still see the events.
    km = nsk.KineticModel.from_antimony(FEEDSPIKE_TOY,
                                        units={'time': 'h', 'conc': 'g/L'})
    _spike().compile(km._te)
    ref = km.run(COND)
    assert ref.final['n_spk'] == 3.0
    assert np.isclose(ref.final['[P]'], 161.9, atol=0.05)
    assert np.allclose(km.copy().run(COND).array, ref.array, rtol=1e-12,
                       atol=0)
    batch = nsk.run_many(km, [COND], metrics={'P': lambda r: r.final['[P]']})
    # A copy matches its original to rtol 1e-12, not bit for bit.
    assert np.isclose(batch.metrics['P'][0], ref.final['[P]'],
                      rtol=1e-12, atol=0)
    flags = km.parameters['assigned_by_events']
    assert flags['n_spk'] and flags['last_vol'] and flags['tot_vol']
    assert not flags['k']
    # A wrapper copied before such an edit is copied afresh after it.
    wrapper = nsk.KineticModel(te.loadAntimonyModel(FEEDSPIKE_TOY),
                               units={'time': 'h', 'conc': 'g/L'})
    before = wrapper.copy().run(COND)
    assert before.final['n_spk'] == 0.0
    _spike().compile(wrapper._te)
    out = wrapper.copy().run(COND)
    assert np.allclose(out.array, wrapper.run(COND).array, rtol=1e-12, atol=0)
    assert np.allclose(out.array, ref.array, rtol=1e-12, atol=0)


def test_copy_follows_structural_edits_and_a_replaced_roadrunner():
    km = _feedspike_model(wrapper=False)
    km._te.addEvent('E9', False, 'time > 1', True)
    km._te.addEventAssignment('E9', 'k', '0', True)
    ref = km.run(COND)
    batch = nsk.run_many(km, [COND], metrics={'P': lambda r: r.final['[P]']})
    # A copy matches its original to rtol 1e-12, not bit for bit.
    assert np.isclose(batch.metrics['P'][0], ref.final['[P]'],
                      rtol=1e-12, atol=0)
    assert ref.final['[P]'] < 50                  # the reaction stopped at t=1
    # init(x) values written to the RoadRunner are part of getSBML().
    km._te['init(env)'] = 2.0
    assert np.allclose(km.copy().run(COND).array, km.run(COND).array,
                       rtol=1e-12, atol=0)
    # A replaced RoadRunner is described and copied as the new model.
    km = nsk.KineticModel.from_antimony(
        'model b(); species S; S = 1; k = 0.2; J: S => ; k*S; end')
    assert list(km.parameters.index) == ['k']
    km._te = te.loadAntimonyModel(
        'model m(); species A; A = 1; kd = 0.1; J: A => ; kd*A; end')
    assert list(km.parameters.index) == ['kd']
    assert list(km.species.index) == ['A']
    cond = nsk.Condition(t_end=1, n_steps=3)
    # A copy matches its original to rtol 1e-12, not necessarily bit for bit.
    assert np.allclose(km.copy().run(cond).array, km.run(cond).array,
                       rtol=1e-12, atol=0)
    assert nsk.run_many(km, [cond]).n_failed == 0


def test_copy_starts_at_the_baseline_with_a_compartment_initial_assignment():
    km = nsk.KineticModel.from_antimony("""
        model ia()
          compartment env; env = 2*k; species S in env; S = 10; k = 0.5;
          n = 0;
          J: S => ; k*S;
          E1: at (time > 1): env = env*1.5, n = n + 1;
        end""")
    km._te['k'] = 0.4
    km.reset()
    clone = km.copy()
    assert clone._te['env'] == km._te['env'] == 0.8
    assert clone._te['[S]'] == km._te['[S]'] == 12.5
    cond = nsk.Condition(t_end=3, n_steps=7, selections=['time', '[S]', 'env',
                                                        'n'])
    assert np.allclose(clone.run(cond).array, km.run(cond).array, rtol=1e-12,
                       atol=0)


def test_copy_of_a_gillespie_model_runs():
    km = nsk.KineticModel.from_antimony(
        'model a(); species S, P; S = 100; P = 0; k = 0.1; '
        'J: S => P; k*S; end')
    km._te.setIntegrator('gillespie')
    km._te.getIntegrator().seed = 42
    clone = km.copy()
    integ = clone._te.getIntegrator()
    assert integ.getName() == 'gillespie' and integ.getValue('seed') == 42
    cond = nsk.Condition(t_end=10, n_steps=11)
    out = clone.run(cond)            # failed with a variant error before
    assert np.array_equal(out.array, km.run(cond).array)


def test_copy_of_a_conserved_moiety_model_runs():
    km = nsk.KineticModel.from_antimony(
        'model c(); species A, B; A = 10; B = 0; k1 = 0.3; k2 = 0.1; '
        'J1: A => B; k1*A; J2: B => A; k2*B; end')
    km._te.conservedMoietyAnalysis = True
    km._te['k1'] = 0.5
    assert '_CSUM0' not in km._user_parameter_values()   # not in the SBML
    clone = km.copy()
    assert clone._te.conservedMoietyAnalysis and clone._te['k1'] == 0.5
    cond = nsk.Condition(t_end=10, n_steps=11)
    assert np.allclose(clone.run(cond).array, km.run(cond).array,
                       rtol=1e-12, atol=0)


def test_copies_keep_attribute_access_on_other_live_models():
    # libroadrunner keeps r.<id> properties on the RoadRunner class and drops
    # the previous model's on every load. The loads copy()/run_many make on
    # the user's behalf must not break `r.kb = ...` on another live model.
    toy_a = nsk.KineticModel.from_antimony(
        'model a(); species A; A = 1; kd = 0.1; J: A => ; kd*A; end')
    toy_b = nsk.KineticModel.from_antimony(      # loaded last: owns r.<id>
        'model b(); species B; B = 1; kb = 0.1; J: B => ; kb*B; end')
    nsk.run_many(toy_a, [nsk.Condition(t_end=1)])
    toy_a.copy()
    toy_b._te.kb = 0.7
    assert toy_b._te['kb'] == 0.7 and 'kb' not in vars(toy_b._te)
    toy_b._te.B = 3.0                        # a species (concentration)
    assert toy_b._te['[B]'] == 3.0 and 'B' not in vars(toy_b._te)


# --- blueprint ------------------------------------------------------------------
def test_blueprint_pickles_and_builds_an_equivalent_model():
    km = _feedspike_model(reset=module_level_reset)
    km.default_n_max = 3
    km._te.k = 0.3
    bp = km._blueprint()
    assert isinstance(bp, _ModelBlueprint)
    bp.check_picklable()
    built = pickle.loads(pickle.dumps(bp)).build()
    assert built._te.k == 0.3 and built.default_n_max == 3
    assert built.reset_func is module_level_reset
    assert np.allclose(built.simulate(0, 48, 481, SELS, reset=True),
                       km.simulate(0, 48, 481, SELS, reset=True),
                       rtol=1e-12, atol=0)
    # building twice gives independent models
    a, b = bp.build(), bp.build()
    assert a._te is not b._te


def test_blueprint_rejects_unpicklable_reset_func_and_attributes():
    km = _feedspike_model(reset=lambda m, **kw: m._te.reset())
    with pytest.raises(ValueError, match='reset_func'):
        km._blueprint().check_picklable()
    km.reset_func = module_level_reset
    km.default_n_max = 3
    km.handle = (lambda: None)
    with pytest.raises(ValueError, match='handle'):
        km._blueprint().check_picklable()
    km.copy()                                  # in-process copy still works


# --- (9) constructors ---------------------------------------------------------------
def test_from_antimony_accepts_text_and_paths(tmp_path):
    path = tmp_path / 'toy.ant'
    path.write_text(FEEDSPIKE_TOY, encoding='utf-8')
    for source in (FEEDSPIKE_TOY, path, str(path)):
        km = nsk.KineticModel.from_antimony(source)
        assert km._units == {'time': 'min', 'conc': 'M'}   # default unchanged
        assert km._te.k == 0.5
    with pytest.raises(TypeError):
        nsk.KineticModel.from_antimony(42)
    # A mistyped path is reported as such, not as an Antimony parse error.
    for missing in ('missing_model.ant', str(tmp_path / 'nope' / 'm.txt'),
                    tmp_path / 'missing.ant'):
        with pytest.raises(FileNotFoundError):
            nsk.KineticModel.from_antimony(missing)
    one_line = nsk.KineticModel.from_antimony('J: S -> P; k*S; S = 1; k = 0.1')
    assert one_line._te.k == 0.1


def test_from_sbml_takes_units_and_reset(tmp_path):
    path = tmp_path / 'toy.xml'
    path.write_text(te.antimonyToSBML(FEEDSPIKE_TOY), encoding='utf-8')
    km = nsk.KineticModel.from_sbml(str(path), units={'time': 'h',
                                                      'conc': 'g/L'},
                                    reset=module_level_reset)
    assert km._units == {'time': 'h', 'conc': 'g/L'}
    assert km.reset_func is module_level_reset
    assert nsk.KineticModel.from_sbml(path)._units == {'time': 'min',
                                                       'conc': 'M'}


# --- (7) overrides ------------------------------------------------------------------
def test_overrides_apply_and_restore():
    km = _feedspike_model()
    with km.overrides({'threshold': 6.0}, k=0.1) as m:
        assert m is km
        assert km._te.k == 0.1 and km._te.threshold == 6.0
    assert km._te.k == 0.5 and km._te.threshold == 5.0


def test_overrides_restore_when_the_body_raises():
    km = _feedspike_model()
    with pytest.raises(RuntimeError, match='boom'):
        with km.overrides(k=0.1):
            raise RuntimeError('boom')
    assert km._te.k == 0.5


@pytest.mark.parametrize('bad', [{'nope': 1.0}, {'curr_env': 2.0},
                                 {'init(curr_env)': 2.0}, {'[nope]': 1.0}])
def test_overrides_reject_bad_ids_before_mutating(bad):
    km = _feedspike_model()
    params = km.parameters['value'].to_dict()
    with pytest.raises(ValueError):
        with km.overrides(dict(bad), k=0.1):
            pass
    assert km.parameters['value'].to_dict() == params


# --- (l) shipped model, slow ----------------------------------------------------------
@pytest.mark.slow
def test_copy_equivalence_on_the_shipped_fed_batch_model():
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
        te_r, reset_nsk_kinetic_model)
    cols = ['time', '[s_glu]', '[x]', '[s_EtOH]', '[s_IBO]', 'env',
            'n_glu_spikes', 'tot_vol_glu_feed_added', 'is_aerobic']
    r = te_r._te

    def run(km):
        km.reset(reset_spike_cap=True)
        # Item access: attribute properties (r.s_glu) belong to whichever
        # model was loaded last.
        km._te['[s_glu]'] = 100.
        km._te['[x]'] = 1.
        return km.simulate(0, 60, 601, cols)

    # te_r is shared by the whole test session and other tests (e.g. the
    # spike-feed reconciliation system) leave parameter edits on it; run from
    # its as-loaded values (init(p): edits never change them) so the result
    # does not depend on test order. overrides() puts the found values back.
    before = te_r._user_parameter_values()
    pristine = {p: r[f'init({p})'] for p in before}
    kla = pristine['kLa'] * 1.25             # a user parameter edit
    try:
        with te_r.overrides({**pristine, 'kLa': kla}):
            ref = run(te_r).copy()
            run(te_r)                        # leave the state dirty
            clone = te_r.copy()
        # overrides restored on the original (a boundary species is written
        # back as an amount in the grown vessel: equal up to rounding)
        assert te_r._user_parameter_values() == pytest.approx(before,
                                                              rel=1e-12)
        assert clone._te['kLa'] == kla
        assert clone.reset_func is reset_nsk_kinetic_model
        assert clone.default_max_n_glu_spikes == te_r.default_max_n_glu_spikes
        assert (clone._te.getIntegrator().getValue('relative_tolerance')
                == r.getIntegrator().getValue('relative_tolerance'))
        out = run(clone)
        assert ref[-1, cols.index('n_glu_spikes')] > 0
        assert np.allclose(out, ref, rtol=1e-12, atol=0)
    finally:
        reset_nsk_kinetic_model(te_r)


@pytest.mark.slow
def test_shipped_model_repeated_runs_do_not_alternate():
    # The shipped reset function used to call RoadRunner's reset() at the
    # previous run's end volume, so runs that do not set [s_glu] themselves
    # alternated between two starting states (15 and 229.65 g/L).
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import te_r
    km = te_r.copy()
    # te_r is shared by the session; other tests (e.g. the spike-feed system)
    # leave parameter edits on it, which the copy carries. Start this private
    # copy from the as-loaded values (init(p)) so the test does not depend on
    # test order.
    for p in km._user_parameter_values():
        km._te[p] = km._te[f'init({p})']
    cond = nsk.Condition(t_end=60, n_steps=121)
    runs = [km.run(cond) for _ in range(3)]
    assert runs[0]['[s_glu]'][0] == 15.0
    for res in runs[1:]:
        assert np.array_equal(res.array, runs[0].array)
    batch = nsk.run_many(te_r, [cond] * 3,
                         metrics={'EtOH': lambda r: r.final['[s_EtOH]']})
    assert batch.n_failed == 0
    assert batch.metrics['EtOH'].nunique() == 1
