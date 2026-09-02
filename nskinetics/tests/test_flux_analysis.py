# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import os

import numpy as np
import tellurium as te
import nskinetics as nsk

# A small model with: an exp product-inhibition term, a saturable denominator
# inhibition term, a threshold-gated (piecewise) enhancement term, a rate-rule
# variable, and two events (a parameter switch and a compartment change). This
# exercises every feature compute_flux_summary must handle, without importing
# the heavy shipped te_r (which would pollute sys.modules and break the
# import-guard tests under `-m "not slow"`).
#
# `T` and `Ptot` are assignment-rule variables, mirroring the shipped model's
# `cell` / `a` / `AcDH`: `T` is a *species* (roadrunner reports it among the
# floating species even though it cannot be set), `Ptot` a parameter. Neither
# appears in any rate law, so they do not perturb the r1/r2/r3 dynamics.
_TOY = """
model flux_toy()
  compartment env; species S in env, P in env, Q in env, T in env;
  r1: S -> P; k1*S*Xa*exp(-ki*P)*exp(-kq*Q)*env;
  r2: P -> Q; sw*k2*P/(P + K2 + Ke*Q)*env;
  r3: S -> ; k3*S*piecewise(exp(kd*(P - Pth)), P > Pth, 1)*env;
  Xa' = -0.01*Xa;
  T := S + P;
  Ptot := P + Q;
  S = 10; P = 0; Q = 0; env = 1; Xa = 1; sw = 1;
  k1 = 0.5; ki = 0.2; kq = 0.1; k2 = 0.3; K2 = 0.1; Ke = 0.5;
  k3 = 0.05; kd = 0.3; Pth = 2;
end
"""


def _make_toy():
    r = te.loadAntimonyModel(_TOY)
    km = nsk.KineticModel(r, units={'time': 'h', 'conc': 'g/L'})
    km.add_event(nsk.Event(when='time >= 5', do={'sw': '0'}, name='switch'))
    km.add_event(nsk.Event(when='time >= 3', do={'env': '1.5'}, name='dilute'))
    km.compile_events()
    return km


def test_state_selections_covers_state_and_excludes_assignment_rules():
    km = _make_toy()
    sels = km.state_selections()
    assert sels[0] == 'time'
    assert 'env' in sels                 # settable compartment
    for s in ('[S]', '[P]', '[Q]'):
        assert s in sels                 # floating species as concentrations
    assert 'Xa' in sels                  # rate-rule variable
    assert 'sw' in sels                  # event-assigned parameter
    # no assignment-rule-defined names, in either bare or bracketed form (the
    # shipped model likewise defines `cell`, `a`, `AcDH` by assignment rules).
    # roadrunner reports an assignment-rule species among the floating species
    # even though writing that selection back raises a RuntimeError:
    assert 'T' not in sels and '[T]' not in sels          # a-rule species
    assert 'Ptot' not in sels and '[Ptot]' not in sels    # a-rule parameter
    # and no duplicates:
    assert len(sels) == len(set(sels))


# NOTE: the public API (``compute_flux_summary``, ``FluxSummary``,
# ``KineticSimulationError``) is imported inside the test bodies, not here:
# ``nskinetics/__init__.py`` imports ``nskinetics.tests``, which imports this
# module, so a module-level ``from nskinetics import ...`` would make the
# package's own initialization order load-bearing.

# Toy inhibition map: r1 is inhibited by P (via ki, "self") and Q (via kq);
# r2 is inhibited by Q via the denominator term Ke; r3 is ENHANCED by P (kd).
_TOY_MAP = {
    'ki': ('r1', 'P'),
    'kq': ('r1', 'Q'),
    'Ke': ('r2', 'Q'),
    'kd': ('r3', 'P'),
}


def _simulate_toy():
    km = _make_toy()
    km.reset()
    sels = km.state_selections() + ['r1', 'r2', 'r3']
    km.simulate(0, 10, 2001, sels)
    return km


def test_reeval_matches_recorded_rates():
    from nskinetics import compute_flux_summary
    km = _simulate_toy()
    s = compute_flux_summary(km, _TOY_MAP, reactions=['r1', 'r2', 'r3'])
    df = km.results_df
    t = df['time'].to_numpy()
    for rid in ('r1', 'r2', 'r3'):
        assert np.isclose(s.cumulative_mass[rid],
                          np.trapezoid(df[rid].to_numpy(), t), rtol=1e-9)


# --- time write-back ---------------------------------------------------------

# A second toy whose rate law DOES read the simulation clock, so `time` must be
# written back on every row for the re-evaluation to reproduce the recorded
# rate. No events, so nothing here depends on the time write-back being safe.
_TOY_TIME = """
model flux_toy_time()
  compartment env; species S in env, P in env;
  r1: S -> P; k1*S*(1 + 0.01*time)*exp(-ki*P)*env;
  S = 10; P = 0; env = 1; k1 = 0.3; ki = 0.2;
end
"""


def _simulate_time_toy():
    km = nsk.KineticModel(te.loadAntimonyModel(_TOY_TIME),
                          units={'time': 'h', 'conc': 'g/L'})
    km.reset()
    km.simulate(0, 10, 2001, km.state_selections() + ['r1'])
    return km


def test_time_is_not_written_back_when_no_rate_law_reads_it():
    # The toy references `time` only in event triggers, so writing it back per
    # row would re-arm compiled native events outside an integration step.
    from nskinetics import compute_flux_summary
    from nskinetics.engine import flux_analysis as fa
    km = _simulate_toy()
    assert 'time' in km.state_selections()      # still recorded and restored
    assert fa._kinetic_laws_use_time(km._te) is False
    seen = []
    original = fa._rates_along

    def _spy(r_, arrs, cols, idx_of, n):
        seen.append(list(cols))
        return original(r_, arrs, cols, idx_of, n)

    fa._rates_along = _spy
    try:
        s = compute_flux_summary(km, _TOY_MAP, reactions=['r1', 'r2', 'r3'])
    finally:
        fa._rates_along = original
    assert seen, 'the write-back loop never ran'
    assert all('time' not in cols for cols in seen)
    assert all('env' in cols and '[S]' in cols for cols in seen)
    # ... and the numbers are unchanged: still the recorded rate integrals.
    df = km.results_df
    t = df['time'].to_numpy()
    for rid in ('r1', 'r2', 'r3'):
        assert np.isclose(s.cumulative_mass[rid],
                          np.trapezoid(df[rid].to_numpy(), t), rtol=1e-9)


def test_time_is_written_back_when_a_rate_law_reads_it():
    from nskinetics import compute_flux_summary
    from nskinetics.engine import flux_analysis as fa
    km = _simulate_time_toy()
    assert fa._kinetic_laws_use_time(km._te) is True
    s = compute_flux_summary(km, {'ki': ('r1', 'P')}, reactions=['r1'])
    df = km.results_df
    # would fail if `time` were held at its restored value during replay
    assert np.isclose(s.cumulative_mass['r1'],
                      np.trapezoid(df['r1'].to_numpy(),
                                   df['time'].to_numpy()), rtol=1e-9)
    assert 0.0 < s.fraction_lost['r1']['P'] < 1.0


def test_mass_balance_closes_on_Q():
    from nskinetics import compute_flux_summary
    km = _simulate_toy()
    s = compute_flux_summary(km, _TOY_MAP, reactions=['r1', 'r2', 'r3'])
    df = km.results_df
    # Q accumulates from r2 only (no Q outflow): cumulative r2 mass == Q*env final
    q_final = df['[Q]'].iloc[-1] * df['env'].iloc[-1]
    assert abs(s.cumulative_mass['r2'] - q_final) / q_final < 1e-3


def test_fraction_lost_signs_and_zero():
    from nskinetics import compute_flux_summary
    km = _simulate_toy()
    s = compute_flux_summary(km, _TOY_MAP, reactions=['r1', 'r2', 'r3'])
    assert 0.0 < s.fraction_lost['r1']['P'] < 1.0     # real inhibition
    assert 0.0 < s.fraction_lost['r1']['Q'] < 1.0
    assert 0.0 < s.fraction_lost['r2']['Q'] < 1.0     # denominator inhibition
    assert s.fraction_lost['r3']['P'] < 0.0           # enhancement -> negative
    assert 0.0 < s.fraction_lost_all['r1'] < 1.0      # joint


def test_state_restored_after_compute():
    from nskinetics import compute_flux_summary
    km = _simulate_toy()
    r = km._te
    # every selection, 'time' included: it is part of the snapshot/restore
    # contract whether or not it is written back per row.
    before = {c: r[c] for c in km.state_selections()}
    assert 'time' in before
    before_params = {p: r[p] for p in _TOY_MAP}
    compute_flux_summary(km, _TOY_MAP, reactions=['r1', 'r2', 'r3'])
    for c, v in before.items():
        assert r[c] == v, f'{c} not restored'
    for p, v in before_params.items():
        assert r[p] == v, f'{p} not restored'


def test_state_restored_when_computation_raises():
    from nskinetics import compute_flux_summary
    # The restore lives in a `finally`, so a failure part-way through must not
    # leave the model parked at some arbitrary trajectory row.
    from nskinetics.engine import flux_analysis as fa
    km = _simulate_toy()
    r = km._te
    before = {c: r[c] for c in km.state_selections()}
    before_params = {p: r[p] for p in _TOY_MAP}

    def _boom(r_, arrs, cols, idx_of, n):
        fa._apply_row(r_, arrs, cols, n // 2)   # move the model off t_end
        raise RuntimeError('boom')

    original = fa._rates_along
    fa._rates_along = _boom
    try:
        compute_flux_summary(km, _TOY_MAP, reactions=['r1', 'r2', 'r3'])
        assert False, 'expected RuntimeError'
    except RuntimeError as e:
        assert 'boom' in str(e)
    finally:
        fa._rates_along = original

    for c, v in before.items():
        assert r[c] == v, f'{c} not restored after failure'
    for p, v in before_params.items():
        assert r[p] == v, f'{p} not restored after failure'


def test_missing_columns_raise():
    from nskinetics import compute_flux_summary
    from nskinetics.exceptions import KineticSimulationError
    km = _make_toy()
    km.reset()
    km.simulate(0, 10, 101, ['time', 'r1'])   # no state columns recorded
    try:
        compute_flux_summary(km, _TOY_MAP, reactions=['r1'])
        assert False, 'expected KineticSimulationError'
    except KineticSimulationError as e:
        assert '[S]' in str(e) or 'env' in str(e)


def test_bad_mapping_raises_before_touching_state():
    from nskinetics import compute_flux_summary
    km = _simulate_toy()
    r = km._te
    s_before = r['[S]']
    try:
        compute_flux_summary(km, {'not_a_param': ('r1', 'X')}, reactions=['r1'])
        assert False, 'expected ValueError'
    except ValueError:
        pass
    assert r['[S]'] == s_before


def test_reactions_subset_ignores_unrequested_mapping_entries():
    from nskinetics import compute_flux_summary
    # A mapping may name reactions outside `reactions`; those entries are
    # validated but not summarized (they must not raise).
    km = _simulate_toy()
    s = compute_flux_summary(km, _TOY_MAP, reactions=['r1'])
    assert s.reaction_ids == ['r1']
    assert set(s.cumulative_mass) == {'r1'}
    assert set(s.fraction_lost) == {'r1'}
    assert set(s.fraction_lost_all) == {'r1'}


def test_csv_roundtrip(tmp_path):
    from nskinetics import compute_flux_summary, FluxSummary
    km = _simulate_toy()
    s = compute_flux_summary(km, _TOY_MAP, reactions=['r1', 'r2', 'r3'],
                             label='toy')
    p = os.path.join(tmp_path, 'summary.csv')
    s.to_csv(p)
    s2 = FluxSummary.from_csv(p)
    assert s2.label == 'toy'
    assert s2.reaction_ids == s.reaction_ids
    assert np.isclose(s2.cumulative_flux['r1'], s.cumulative_flux['r1'])
    assert np.isclose(s2.fraction_lost['r1']['P'], s.fraction_lost['r1']['P'])
    assert np.isclose(s2.fraction_lost_all['r1'], s.fraction_lost_all['r1'])
    assert np.isclose(s2.final_volume, s.final_volume)
    assert np.isclose(s2.t_end, s.t_end)
    assert s2.inhibitors == s.inhibitors


def test_csv_roundtrip_partial_map(tmp_path):
    from nskinetics import compute_flux_summary, FluxSummary
    # A reaction with no entry in the inhibition map is absent from
    # fraction_lost / fraction_lost_all; the CSV must round-trip that absence
    # rather than inventing a nan entry. Also covers label=None, and the two
    # fields the first round-trip test does not assert on.
    km = _simulate_toy()
    s = compute_flux_summary(km, {'ki': ('r1', 'P')}, reactions=['r1', 'r2'])
    assert 'r2' not in s.fraction_lost          # producer-side precondition
    assert 'r2' not in s.fraction_lost_all
    p = os.path.join(tmp_path, 'partial.csv')
    s.to_csv(p)
    s2 = FluxSummary.from_csv(p)
    assert 'r2' not in s2.fraction_lost         # absent key, not nan
    assert 'r2' not in s2.fraction_lost_all
    assert s2.label is None
    assert s2.reaction_ids == s.reaction_ids    # order preserved
    # values are written at full repr, so these compare exactly
    assert s2.cumulative_mass == s.cumulative_mass
    assert s2.final_concentrations == s.final_concentrations
    assert np.isclose(s2.fraction_lost['r1']['P'], s.fraction_lost['r1']['P'])
    assert np.isclose(s2.fraction_lost_all['r1'], s.fraction_lost_all['r1'])
    assert s2.inhibitors == s.inhibitors == ['P']


# --- integration window (t_end) ---------------------------------------------

def test_t_end_truncates_the_integration_window():
    from nskinetics import compute_flux_summary
    # Q accumulates from r2 only, so the cumulative r2 mass up to t = 5 must
    # close against Q*env at the t = 5 row -- not at the end of the run.
    km = _simulate_toy()
    s = compute_flux_summary(km, _TOY_MAP, reactions=['r1', 'r2', 'r3'],
                             t_end=5.0)
    assert np.isclose(s.t_end, 5.0) and s.t_end <= 5.0
    df = km.results_df
    i = int(np.flatnonzero(df['time'].to_numpy() <= 5.0)[-1])
    q_at_t_end = df['[Q]'].iat[i] * df['env'].iat[i]
    assert abs(s.cumulative_mass['r2'] - q_at_t_end) / q_at_t_end < 1e-3
    # the window really is shorter than the run
    assert s.final_volume == df['env'].iat[i]
    assert s.final_concentrations['Q'] == df['[Q]'].iat[i]


def test_t_end_default_is_the_whole_trajectory():
    from nskinetics import compute_flux_summary
    km = _simulate_toy()
    full = compute_flux_summary(km, _TOY_MAP, reactions=['r1', 'r2', 'r3'])
    explicit = compute_flux_summary(km, _TOY_MAP, reactions=['r1', 'r2', 'r3'],
                                    t_end=None)
    assert full.t_end == explicit.t_end == km.results_df['time'].iat[-1]
    assert full.cumulative_mass == explicit.cumulative_mass
    # and a t_end at/after the last row changes nothing
    late = compute_flux_summary(km, _TOY_MAP, reactions=['r1', 'r2', 'r3'],
                                t_end=1e6)
    assert late.cumulative_mass == full.cumulative_mass


def test_t_end_before_the_first_row_raises():
    from nskinetics import compute_flux_summary
    km = _simulate_toy()
    r = km._te
    s_before = r['[S]']
    try:
        compute_flux_summary(km, _TOY_MAP, reactions=['r1'], t_end=-1.0)
        assert False, 'expected ValueError'
    except ValueError as e:
        assert 't_end' in str(e)
    assert r['[S]'] == s_before          # raised before touching state
