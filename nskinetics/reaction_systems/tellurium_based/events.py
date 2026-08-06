# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
import itertools

from ...exceptions import EventCompilationError

__all__ = ('Event', 'FeedSpike')

_event_counter = itertools.count()


def _next_event_name(prefix='nsk_event'):
    return f'{prefix}_{next(_event_counter)}'


class Event:
    """
    A single SBML event, compiled into a native roadrunner event.

    Parameters
    ----------
    when : str
        Trigger expression, e.g. ``'time >= stage_1_max_time'`` or
        ``'s_glu <= threshold'``.
    do : dict
        Mapping of ``{var: expr}`` event assignments, e.g. ``{'is_aerobic': '0'}``.
        Every ``var`` must be a **non-constant** model variable.
    delay : str or float, optional
        Event delay; a number or the name of a model parameter.
    priority : str or float, optional
        Event priority; a number or the name of a model parameter. Higher fires
        first when multiple events trigger simultaneously.
    name : str, optional
        Event id. Auto-generated if ``None``.
    use_trigger_time_values : bool
        SBML ``useValuesFromTriggerTime``. Default ``False``.
    """
    def __init__(self, when, do, *, delay=None, priority=None, name=None,
                 use_trigger_time_values=False):
        self.when = when
        self.do = dict(do)
        self.delay = delay
        self.priority = priority
        self.name = name if name is not None else _next_event_name()
        self.use_trigger_time_values = use_trigger_time_values

    def compile(self, r, force_regenerate=True):
        """Inject this event into roadrunner model ``r``.

        Parameters
        ----------
        r : roadrunner.RoadRunner
            The model to inject into.
        force_regenerate : bool
            If ``True``, regenerate the model immediately. Pass ``False`` to
            batch multiple events and regenerate once afterward.
        """
        name = self.name
        try:
            r.addEvent(name, self.use_trigger_time_values, self.when, False)
            for var, expr in self.do.items():
                r.addEventAssignment(name, var, str(expr), False)
            if self.delay is not None:
                r.addDelay(name, str(self.delay), False)
            if self.priority is not None:
                r.addPriority(name, str(self.priority), False)
        except Exception as e:
            raise EventCompilationError(
                f'Failed to compile event {name!r} (when={self.when!r}, '
                f'do={self.do!r}): {e}') from e
        if force_regenerate:
            r.regenerateModel()


class FeedSpike:
    """Placeholder; implemented in a later task."""
    def __init__(self, *args, **kwargs):
        raise NotImplementedError('FeedSpike not yet implemented.')
