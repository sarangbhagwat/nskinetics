# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

__all__ = ('NSKError', 'KineticSimulationError', 'MassBalanceError',
           'EventCompilationError')


class NSKError(Exception):
    """Base class for all NSKinetics-raised errors."""


class KineticSimulationError(NSKError):
    """Raised when a kinetic simulation fails or produces an invalid state."""


class MassBalanceError(NSKError):
    """Raised when a post-simulation mass/atom balance check fails."""


class EventCompilationError(NSKError):
    """Raised when compiling a Python Event into a native SBML event fails."""
