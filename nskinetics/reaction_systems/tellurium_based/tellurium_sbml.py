# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import tellurium as te

from ...exceptions import KineticSimulationError, EventCompilationError
from .events import _regenerate_and_resync

__all__ = ('TelluriumReactionSystem',)

_TIME_FACTORS = {'h': 1.0, 'hr': 1.0, 'min': 60.0, 'm': 60.0, 's': 3600.0, 'sec': 3600.0}
# 'kg/m3' (== 'g/L') is a *mass* concentration; it must live only in
# _MASS_CONC_UNITS. material_indexer checks the molar tuple first, so listing it
# in both would mis-resolve a mass unit to 'imol'.
_MOLAR_CONC_UNITS = ('M', 'mol/L')
_MASS_CONC_UNITS = ('g/L', 'kg/m3', 'kg/m^3')


class TelluriumReactionSystem():
    """
    Thin wrapper around a Tellurium extended RoadRunner model, adding a Python
    event API and unit-aware value access.

    Parameters
    ----------
    te : tellurium.roadrunner.extended_roadrunner.ExtendedRoadRunner
        A Tellurium extended RoadRunner object.
    units : dict, optional
        Mapping with keys ``'time'`` and ``'conc'``. Defaults to
        ``{'time': 'min', 'conc': 'M'}``.
    """
    def __init__(self, te, units=None):
        self._te = te
        if units is not None:
            self._units = units
        else:
            self._units = {'time': 'min', 'conc': 'M'}
        self.events = []
        self._events_compiled = False

    @classmethod
    def from_sbml(cls, filepath):
        return cls(te.loadSBMLModel(filepath))

    # --- value access -------------------------------------------------------
    def get_value(self, selection):
        """Return the current value of a roadrunner selection (plain name or
        ``'[conc]'``)."""
        return self._te[selection]

    def set_value(self, selection, value):
        """Set a roadrunner selection to ``value``."""
        self._te[selection] = value

    def set_concentration_from_stream(self, var, amount, volume):
        """Set the bracket-stripped model variable ``var`` to ``amount/volume``
        and return the value. Replicates the legacy concentration assignment."""
        plain = var.replace('[', '').replace(']', '')
        value = amount / volume
        self._te[plain] = value
        return value

    # --- units --------------------------------------------------------------
    @property
    def time_factor(self):
        """Hours per model time unit."""
        time = self._units.get('time', '')
        try:
            return _TIME_FACTORS[time.lower()]
        except (AttributeError, KeyError):
            raise KineticSimulationError(
                f'Unrecognized time units {time!r}; '
                f'expected one of {sorted(_TIME_FACTORS)}.')

    @property
    def material_indexer(self):
        """biosteam stream indexer matching the model's concentration units."""
        conc = self._units.get('conc')
        if conc in _MOLAR_CONC_UNITS:
            return 'imol'
        elif conc in _MASS_CONC_UNITS:
            return 'imass'
        raise KineticSimulationError(f'Unrecognized concentration units {conc!r}.')

    def validate_units(self):
        """Raise ``KineticSimulationError`` if time/conc units are unrecognized."""
        time_u = self._units.get('time', '').lower()
        if time_u not in _TIME_FACTORS:
            raise KineticSimulationError(
                f'Unrecognized time units {self._units.get("time")!r}; '
                f'expected one of {sorted(_TIME_FACTORS)}.')
        conc_u = self._units.get('conc')
        if conc_u not in _MOLAR_CONC_UNITS and conc_u not in _MASS_CONC_UNITS:
            raise KineticSimulationError(
                f'Unrecognized concentration units {conc_u!r}.')

    # --- events -------------------------------------------------------------
    def add_event(self, event):
        """Append an :class:`Event` to the pending list (compiled by
        :meth:`compile_events`). Returns the event."""
        self.events.append(event)
        self._events_compiled = False
        return event

    def remove_event(self, name):
        """Remove a pending event by name."""
        self.events = [e for e in self.events if e.name != name]

    def clear_events(self):
        """Drop all pending events."""
        self.events = []

    def compile_events(self):
        """Inject all pending events into the model, then regenerate once.

        Guarded so repeated calls do not double-inject. Raises
        :class:`EventCompilationError` with context on failure.

        Notes
        -----
        Compiling events regenerates the model and resets all model state to
        its SBML-defined origin values as a side effect. Set any custom
        initial conditions or parameter overrides AFTER compiling events, not
        before, or they will be discarded.
        """
        if self._events_compiled:
            return
        r = self._te
        try:
            for event in self.events:
                event.compile(r, force_regenerate=False)
            _regenerate_and_resync(r)
        except EventCompilationError:
            raise
        except Exception as e:
            raise EventCompilationError(
                f'Failed to compile {len(self.events)} event(s): {e}') from e
        self._events_compiled = True

    # --- reset --------------------------------------------------------------
    def reset(self, **kwargs):
        """Reset model state (not structure); extra kwargs are ignored so
        callers may pass bookkeeping flags."""
        self._te.reset()
