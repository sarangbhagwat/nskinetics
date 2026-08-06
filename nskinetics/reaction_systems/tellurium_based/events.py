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

        Notes
        -----
        When ``force_regenerate`` is ``True``, compiling this event regenerates
        the model and resets all model state to its SBML-defined origin values
        as a side effect. Set any custom initial conditions or parameter
        overrides AFTER compiling events, not before, or they will be
        discarded.
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
            # roadrunner's regenerateModel() leaves the model's stored
            # init(X) bookkeeping stale for any rate-rule-governed variable
            # with an explicit initial value; a plain r.reset() afterward
            # then reads back corrupted values. resetToOrigin() re-syncs
            # init(X) to the SBML-defined values so later reset() calls are
            # correct.
            r.resetToOrigin()


class FeedSpike:
    """
    High-level fed-batch spike: expands to three coordinated SBML events that
    add feed to hit a target concentration, update the working volume, and
    increment a spike counter.

    Parameters
    ----------
    species : str
        Model variable whose concentration is spiked (e.g. ``'s_glu'``).
    when : str
        Base trigger expression (e.g. ``'s_glu <= threshold_conc_glu_spike'``).
    target : str or float
        Target concentration after the spike; a number or model-parameter name.
    feed_conc : str or float
        Concentration of the feed being added; a number or parameter name.
    volume_var : str
        Working-volume variable the reactor reads to scale effluent (default
        ``'env'``).
    max_count : str, float, or None
        Optional cap on the number of spikes; a number or parameter name. When
        given, ``&& (count_var < max_count)`` is folded into every trigger.
    count_var : str, optional
        Spike-counter variable. Defaults to ``f'{species}_spike_count'``.
    last_vol_var : str, optional
        Bookkeeping var for the last volume added. Defaults to
        ``f'last_vol_{species}_added'``.
    tot_vol_var : str, optional
        Bookkeeping var for total volume added. Defaults to
        ``f'tot_vol_{species}_added'``.
    delay : str or float, optional
        Event delay; a number or parameter name.
    priority : int
        Priority of the first (``_a``) event; ``_b`` and ``_c`` get
        ``priority-1`` and ``priority-2`` to enforce ordering. Default ``5``.
    name : str, optional
        Base event name; the three events are ``{name}_a/_b/_c``. Auto-generated
        if ``None``.
    """
    def __init__(self, species, when, target, feed_conc, *, volume_var='env',
                 max_count=None, count_var=None, last_vol_var=None,
                 tot_vol_var=None, delay=None, priority=5, name=None):
        self.species = species
        self.when = when
        self.target = target
        self.feed_conc = feed_conc
        self.volume_var = volume_var
        self.max_count = max_count
        self.count_var = count_var if count_var is not None else f'{species}_spike_count'
        self.last_vol_var = last_vol_var if last_vol_var is not None else f'last_vol_{species}_added'
        self.tot_vol_var = tot_vol_var if tot_vol_var is not None else f'tot_vol_{species}_added'
        self.delay = delay
        self.priority = priority
        self.name = name if name is not None else _next_event_name('feed_spike')

    def _trigger(self):
        if self.max_count is not None:
            return f'({self.when}) && ({self.count_var} < {self.max_count})'
        return self.when

    def expand(self):
        """Return the three coordinated :class:`Event` objects."""
        s, vol = self.species, self.volume_var
        tgt, feed = self.target, self.feed_conc
        last_vol, tot_vol, count = self.last_vol_var, self.tot_vol_var, self.count_var
        trig = self._trigger()
        common = dict(when=trig, delay=self.delay)
        a = Event(
            do={last_vol: f'{vol}*({tgt} - {s})/({feed} - {tgt})'},
            priority=self.priority, name=f'{self.name}_a',
            use_trigger_time_values=True, **common)
        b = Event(
            do={vol: f'{vol} + {last_vol}',
                tot_vol: f'{tot_vol} + {last_vol}'},
            priority=self.priority - 1, name=f'{self.name}_b',
            use_trigger_time_values=False, **common)
        c = Event(
            do={s: f'{tgt}', count: f'{count} + 1'},
            priority=self.priority - 2, name=f'{self.name}_c',
            use_trigger_time_values=False, **common)
        return [a, b, c]

    def compile(self, r, force_regenerate=True):
        """Expand and compile all three events into roadrunner model ``r``.

        Notes
        -----
        When ``force_regenerate`` is ``True``, compiling this spike regenerates
        the model and resets all model state to its SBML-defined origin values
        as a side effect. Set any custom initial conditions or parameter
        overrides AFTER compiling events, not before, or they will be
        discarded.
        """
        events = self.expand()
        for event in events:
            event.compile(r, force_regenerate=False)
        if force_regenerate:
            r.regenerateModel()
            # See the comment in Event.compile() for why this is needed.
            r.resetToOrigin()
