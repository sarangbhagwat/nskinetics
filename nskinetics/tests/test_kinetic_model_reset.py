# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import nskinetics as nsk


class _StubTe:
    """Minimal stand-in for a Tellurium RoadRunner: counts reset() calls.

    KineticModel.__init__ only stores the object; nothing else here
    touches it, so a stub suffices to unit-test the reset dispatch without
    loading a real model.
    """
    def __init__(self):
        self.reset_calls = 0

    def reset(self):
        self.reset_calls += 1


def test_default_reset_calls_underlying_te_reset():
    stub = _StubTe()
    km = nsk.KineticModel(stub, units={'time': 'h', 'conc': 'g/L'})
    km.reset()
    assert stub.reset_calls == 1


def test_default_reset_ignores_extra_kwargs():
    stub = _StubTe()
    km = nsk.KineticModel(stub)
    km.reset(reset_spike_cap=True, anything=123)  # must not raise
    assert stub.reset_calls == 1


def test_custom_reset_receives_trs_and_replaces_default():
    stub = _StubTe()
    seen = {}

    def my_reset(km, **kwargs):
        seen['km'] = km
        seen['kwargs'] = kwargs

    km = nsk.KineticModel(stub, reset=my_reset)
    km.reset(reset_spike_cap=False)
    assert seen['km'] is km
    assert seen['kwargs'] == {'reset_spike_cap': False}
    # Full replacement: the default self._te.reset() is NOT auto-called.
    assert stub.reset_calls == 0


def test_reset_func_settable_after_construction():
    stub = _StubTe()
    km = nsk.KineticModel(stub)
    seen = []
    km.reset_func = lambda t, **kw: seen.append(t)
    km.reset()
    assert seen == [km]
    assert stub.reset_calls == 0
    assert km.reset_func is not None
