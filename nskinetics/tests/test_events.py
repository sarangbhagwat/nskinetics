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
