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
