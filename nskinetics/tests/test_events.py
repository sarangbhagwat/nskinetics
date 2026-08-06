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
