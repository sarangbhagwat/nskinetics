# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

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


from nskinetics import compute_flux_summary, FluxSummary
from nskinetics.exceptions import KineticSimulationError

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
    km = _simulate_toy()
    s = compute_flux_summary(km, _TOY_MAP, reactions=['r1', 'r2', 'r3'])
    df = km.results_df
    t = df['time'].to_numpy()
    for rid in ('r1', 'r2', 'r3'):
        assert np.isclose(s.cumulative_mass[rid],
                          np.trapezoid(df[rid].to_numpy(), t), rtol=1e-9)


def test_mass_balance_closes_on_Q():
    km = _simulate_toy()
    s = compute_flux_summary(km, _TOY_MAP, reactions=['r1', 'r2', 'r3'])
    df = km.results_df
    # Q accumulates from r2 only (no Q outflow): cumulative r2 mass == Q*env final
    q_final = df['[Q]'].iloc[-1] * df['env'].iloc[-1]
    assert abs(s.cumulative_mass['r2'] - q_final) / q_final < 1e-3


def test_fraction_lost_signs_and_zero():
    km = _simulate_toy()
    s = compute_flux_summary(km, _TOY_MAP, reactions=['r1', 'r2', 'r3'])
    assert 0.0 < s.fraction_lost['r1']['P'] < 1.0     # real inhibition
    assert 0.0 < s.fraction_lost['r1']['Q'] < 1.0
    assert 0.0 < s.fraction_lost['r2']['Q'] < 1.0     # denominator inhibition
    assert s.fraction_lost['r3']['P'] < 0.0           # enhancement -> negative
    assert 0.0 < s.fraction_lost_all['r1'] < 1.0      # joint


def test_state_restored_after_compute():
    km = _simulate_toy()
    r = km._te
    before = {c: r[c] for c in km.state_selections() if c != 'time'}
    before_params = {p: r[p] for p in _TOY_MAP}
    compute_flux_summary(km, _TOY_MAP, reactions=['r1', 'r2', 'r3'])
    for c, v in before.items():
        assert r[c] == v, f'{c} not restored'
    for p, v in before_params.items():
        assert r[p] == v, f'{p} not restored'


def test_missing_columns_raise():
    km = _make_toy()
    km.reset()
    km.simulate(0, 10, 101, ['time', 'r1'])   # no state columns recorded
    try:
        compute_flux_summary(km, _TOY_MAP, reactions=['r1'])
        assert False, 'expected KineticSimulationError'
    except KineticSimulationError as e:
        assert '[S]' in str(e) or 'env' in str(e)


def test_bad_mapping_raises_before_touching_state():
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
    # A mapping may name reactions outside `reactions`; those entries are
    # validated but not summarized (they must not raise).
    km = _simulate_toy()
    s = compute_flux_summary(km, _TOY_MAP, reactions=['r1'])
    assert s.reaction_ids == ['r1']
    assert set(s.cumulative_mass) == {'r1'}
    assert set(s.fraction_lost) == {'r1'}
    assert set(s.fraction_lost_all) == {'r1'}
