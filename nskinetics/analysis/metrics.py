# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Kinetics-level metrics: scalar summaries of one simulated trajectory.

A metric maps a :class:`~nskinetics.SimulationResult` to one float and
declares the columns it reads, so it plugs directly into
:func:`nskinetics.run_many` / :func:`nskinetics.scan` (which record the
``selections`` it declares), :func:`~nskinetics.analysis.local_sensitivity`
(``metrics=``) and the screening and calibration tools. Every metric states
its basis explicitly -- which columns, which time points, amounts or
concentrations -- and nothing is hard-coded: no species names, no ``'env'``
volume, no yields.

All arithmetic runs under ``numpy.errstate(all='ignore')`` (importing
nskinetics sets numpy's divide/invalid handling to ``'raise'``), and a zero
denominator gives ``NaN`` rather than an error or an infinity.
"""
import math

import numpy as np
import pandas as pd

from ..engine.conditions import SimulationResult
from ..engine.execution import _scalar
from ._common import _finite

__all__ = ('Metric', 'Final', 'At', 'Max', 'Min', 'TimeAtMax', 'Integral',
           'EventCount', 'Productivity', 'Yield', 'CallableMetric')

_YIELD_BASES = ('consumed', 'added')


# --- validation helpers ---------------------------------------------------------
def _selection(arg, value):
    """``value`` as a selection string, or ``TypeError``/``ValueError``."""
    if not isinstance(value, str):
        raise TypeError(f'{arg} must be a selection string (e.g. "[P]"); got '
                        f'{type(value).__name__}.')
    if not value or value != value.strip():
        raise ValueError(f'{arg} must be a non-empty selection without '
                         f'surrounding spaces; got {value!r}.')
    if value == 'time':
        raise ValueError(f"{arg} cannot be 'time'.")
    return value


def _ratio(numerator, denominator):
    """``numerator / denominator`` as a float; ``NaN`` for a zero or NaN
    denominator (never an infinity or a ``ZeroDivisionError``)."""
    numerator, denominator = float(numerator), float(denominator)
    if denominator == 0.0 or math.isnan(denominator):
        return math.nan
    with np.errstate(all='ignore'):
        return float(np.float64(numerator) / np.float64(denominator))


def _fmt(value):
    """Short, exact-enough text for a number inside a default metric name."""
    return f'{value:g}'


# --- base class -----------------------------------------------------------------
class Metric:
    """Base class of the metric protocol: a named scalar summary of a run.

    A metric is any object with

    - ``name`` (str): the column name it gets in batch tables;
    - ``selections`` (tuple of str): the columns of the trajectory it reads
      (:func:`nskinetics.run_many` records them in addition to the
      condition's own selections; ``'time'`` is always recorded);
    - ``__call__(result) -> float`` for a
      :class:`~nskinetics.SimulationResult`.

    The concrete metrics in this module subclass :class:`Metric`, which
    evaluates them under ``numpy.errstate(all='ignore')`` and returns a
    Python float. Use :class:`CallableMetric` (or :meth:`from_function`) to
    wrap your own function.

    Parameters
    ----------
    name : str
        Metric name (non-empty).
    selections : sequence of str
        Columns the metric reads.

    Raises
    ------
    TypeError, ValueError
        For an empty or non-string name or selection.

    Notes
    -----
    Metrics read the recorded output grid only. Values between output
    times are linearly interpolated (:class:`At`), integrals use the
    trapezoid rule over the grid (:class:`Integral`), extrema are taken
    over the grid points (:class:`Max`, :class:`Min`) and the time of the
    maximum is refined by a parabola through the grid maximum and its two
    neighbours (:class:`TimeAtMax`); request a finer grid (``n_steps`` or
    ``times``) where that matters, e.g. around feed spikes. Two metrics are
    equal when they are of the same type with the same arguments.

    Local sensitivities difference a metric as computed, so a metric that
    is piecewise constant in the parameters -- :class:`EventCount`, or a
    :class:`CallableMetric` that returns a grid time -- has a derivative of
    0, or a jump; screen such metrics with
    :func:`~nskinetics.analysis.morris` / :func:`~nskinetics.analysis.sobol`
    instead.
    """
    def __init__(self, name, selections=()):
        if not isinstance(name, str) or not name:
            raise TypeError(f'Metric name must be a non-empty string; got '
                            f'{name!r}.')
        if isinstance(selections, str):
            selections = (selections,)
        sels = []
        for s in selections:
            _selection('selections', s)
            if s not in sels:
                sels.append(s)
        self.name = name
        self.selections = tuple(sels)

    def __call__(self, result):
        """Evaluate the metric on one trajectory.

        Parameters
        ----------
        result : SimulationResult or pandas.DataFrame
            The run (a DataFrame must have a ``'time'`` column).

        Returns
        -------
        float
            The metric value; ``NaN`` where a denominator is zero.

        Raises
        ------
        TypeError
            If ``result`` is neither a SimulationResult nor a DataFrame.
        KeyError
            If a column the metric reads was not recorded.
        """
        if isinstance(result, pd.DataFrame):
            result = SimulationResult(result)
        elif not isinstance(result, SimulationResult):
            raise TypeError(f'{type(self).__name__} metrics evaluate a '
                            f'SimulationResult; got {type(result).__name__}.')
        with np.errstate(all='ignore'):
            return _scalar(self.name, self._evaluate(result))

    def _evaluate(self, result):        # abstract
        raise NotImplementedError(
            f'{type(self).__name__} does not compute a value; use a concrete '
            'metric (Final, Yield, ...) or wrap a function with '
            'CallableMetric / Metric.from_function.')

    def _check_condition(self, condition):
        """Raise ``ValueError`` if no run of ``condition`` can be evaluated
        (e.g. :class:`At` outside its time span). Called by the analyses
        before any run; the base class accepts every condition."""

    #: Names of the constructor arguments stored as attributes (used by
    #: ``==``, ``hash`` and ``repr``).
    _fields = ()

    def _args(self):
        """The arguments that define the metric."""
        return (self.name, self.selections,
                *(getattr(self, f) for f in self._fields))

    def __eq__(self, other):
        if type(other) is not type(self):
            return NotImplemented
        return self._args() == other._args()

    def __hash__(self):
        return hash((type(self).__name__, self._args()))

    def __repr__(self):
        args = ''.join(f'{f}={getattr(self, f)!r}, ' for f in self._fields)
        return f'{type(self).__name__}({args}name={self.name!r})'

    @classmethod
    def from_function(cls, fn, name, selections=()):
        """Wrap a function as a metric (same as :class:`CallableMetric`).

        Parameters
        ----------
        fn : callable
            ``fn(result) -> float`` for a :class:`~nskinetics.SimulationResult`.
        name : str
            Metric name.
        selections : sequence of str, optional
            Columns ``fn`` reads, recorded by batch runs.

        Returns
        -------
        CallableMetric

        Raises
        ------
        TypeError, ValueError
            As for :class:`CallableMetric`.
        """
        return CallableMetric(fn, name, selections)


def _time_span(condition):
    """``(first, last)`` output time of a run of ``condition`` (model time
    units): its ``times``, else ``t0 .. t_end`` (a Condition has exactly
    one of them)."""
    if condition.times is not None:
        return float(condition.times[0]), float(condition.times[-1])
    return float(condition.t0), float(condition.t_end)


def _check_metrics_condition(metrics, condition):
    """Call each metric's ``_check_condition`` hook (if it has one) on
    ``condition``; ``ValueError`` prefixed with ``'metrics: '``.

    Parameters
    ----------
    metrics : mapping
        ``{name: metric}`` (normalized, e.g. by ``run_many``'s validation).
    condition : Condition
    """
    for metric in metrics.values():
        check = getattr(metric, '_check_condition', None)
        if check is None:
            continue
        try:
            check(condition)
        except ValueError as e:
            raise ValueError(f'metrics: {e}') from e


def _init_single(metric, selection, name):
    """Shared constructor of the metrics of one ``selection``."""
    metric.selection = _selection('selection', selection)
    Metric.__init__(metric, name if name is not None
                    else f'{type(metric).__name__}({selection})',
                    (selection,))


# --- concrete metrics -----------------------------------------------------------
class Final(Metric):
    """Value of a selection at the last output time.

    Parameters
    ----------
    selection : str
        Column to read, e.g. ``'[P]'`` (concentration, model units).
    name : str, optional
        Defaults to ``'Final(<selection>)'``.

    Raises
    ------
    TypeError, ValueError
        For a selection that is not a non-empty string, or is ``'time'``.

    Notes
    -----
    Basis: ``y(t_end)``, the last recorded row, in the selection's own units.
    """
    _fields = ('selection',)

    def __init__(self, selection, *, name=None):
        _init_single(self, selection, name)

    def _evaluate(self, result):
        return result[self.selection][-1]


class Max(Metric):
    """Largest value of a selection over the output grid.

    Parameters
    ----------
    selection : str
        Column to read.
    name : str, optional
        Defaults to ``'Max(<selection>)'``.

    Raises
    ------
    TypeError, ValueError
        For a selection that is not a non-empty string, or is ``'time'``.

    Notes
    -----
    Basis: the maximum over the recorded rows (not between them), in the
    selection's units; ``NaN`` if any value is ``NaN``.
    """
    _fields = ('selection',)

    def __init__(self, selection, *, name=None):
        _init_single(self, selection, name)

    def _evaluate(self, result):
        return np.max(result[self.selection])


class Min(Metric):
    """Smallest value of a selection over the output grid.

    Parameters
    ----------
    selection : str
        Column to read.
    name : str, optional
        Defaults to ``'Min(<selection>)'``.

    Raises
    ------
    TypeError, ValueError
        For a selection that is not a non-empty string, or is ``'time'``.

    Notes
    -----
    Basis: the minimum over the recorded rows, in the selection's units;
    ``NaN`` if any value is ``NaN``.
    """
    _fields = ('selection',)

    def __init__(self, selection, *, name=None):
        _init_single(self, selection, name)

    def _evaluate(self, result):
        return np.min(result[self.selection])


class TimeAtMax(Metric):
    """Time at which a selection first reaches its maximum.

    Parameters
    ----------
    selection : str
        Column to read.
    name : str, optional
        Defaults to ``'TimeAtMax(<selection>)'``.

    Raises
    ------
    TypeError, ValueError
        For a selection that is not a non-empty string, or is ``'time'``.

    Notes
    -----
    Basis (model time units): the vertex of the parabola through the first
    grid maximum and its two neighbouring rows, which lies within one
    output interval of it; the grid time itself when the maximum is the
    first or last row or is tied with the next row (a plateau). ``NaN`` if
    any value is ``NaN``. The vertex locates a smooth peak to ``O(dt^2)``
    in the output spacing ``dt`` and, unlike the grid time, moves
    continuously with the parameters, so finite-difference sensitivities
    of it are meaningful. Its parameter derivative, however, converges
    only to first order in ``dt``: it oscillates with the peak's position
    within a grid cell. For an A -> B -> C intermediate (``k1 = 1``,
    ``k2 = 0.5``, peak near 1.4) the worst relative error of
    ``d TimeAtMax / d k1`` over ``k1`` in ``[0.9, 1.1]`` is about 7% at
    ``dt = 0.1``, 3% at ``dt = 0.05`` and 0.6% at ``dt = 0.01``, so use a
    dense output grid around the peak when its sensitivities matter. A
    maximum at a kink, e.g. just before a feed spike, is only located to
    within the grid.
    """
    _fields = ('selection',)

    def __init__(self, selection, *, name=None):
        _init_single(self, selection, name)

    def _evaluate(self, result):
        y, t = result[self.selection], result.time
        if np.isnan(y).any():
            return math.nan
        i = int(np.argmax(y))
        # The first argmax has y[i-1] < y[i]; a tie with y[i+1] is a plateau.
        if i == 0 or i == len(y) - 1 or y[i + 1] == y[i]:
            return t[i]
        t0, t1, t2 = t[i - 1:i + 2]
        y0, y1, y2 = y[i - 1:i + 2]
        with np.errstate(all='ignore'):
            num = (t1 - t0) ** 2 * (y1 - y2) - (t1 - t2) ** 2 * (y1 - y0)
            den = (t1 - t0) * (y1 - y2) - (t1 - t2) * (y1 - y0)
            vertex = t1 - 0.5 * num / den
        if den == 0 or not np.isfinite(vertex):
            return t[i]
        return min(max(vertex, t0), t2)


class Integral(Metric):
    """Time integral of a selection over the run (trapezoid rule).

    Parameters
    ----------
    selection : str
        Column to integrate.
    name : str, optional
        Defaults to ``'Integral(<selection>)'``.

    Raises
    ------
    TypeError, ValueError
        For a selection that is not a non-empty string, or is ``'time'``.

    Notes
    -----
    Basis: ``∫ y dt`` from the first to the last output time by the
    trapezoid rule on the recorded grid, in the selection's units times
    model time units (e.g. g/L·h). Jumps at events are interpolated across,
    so resolve them with a fine grid.
    """
    _fields = ('selection',)

    def __init__(self, selection, *, name=None):
        _init_single(self, selection, name)

    def _evaluate(self, result):
        return np.trapezoid(result[self.selection], result.time)


class EventCount(Metric):
    """Number of increments of an event counter during the run.

    Parameters
    ----------
    counter : str
        A counter variable, e.g. the ``count_var`` of a
        :class:`~nskinetics.FeedSpike` (``'n_spk'``) or any variable an event
        increments by one (``n = n + 1``).
    name : str, optional
        Defaults to ``'EventCount(<counter>)'``.

    Raises
    ------
    TypeError, ValueError
        For a counter that is not a non-empty string, or is ``'time'``.

    Notes
    -----
    Basis: final minus initial recorded value of ``counter``, i.e. the
    number of firings of an event that adds 1 per firing, independent of
    the counter's starting value. The count is piecewise constant in the
    parameters: its local derivative is 0 unless a perturbation changes the
    count, which :func:`~nskinetics.analysis.local_sensitivity` flags
    (``event_count_changed``) when it watches the counter. Use the
    screening functions for its dependence on parameters.
    """
    _fields = ('counter',)

    def __init__(self, counter, *, name=None):
        self.counter = _selection('counter', counter)
        super().__init__(name if name is not None
                         else f'EventCount({counter})', (counter,))

    def _evaluate(self, result):
        y = result[self.counter]
        return y[-1] - y[0]


class At(Metric):
    """Value of a selection at a given time.

    Parameters
    ----------
    selection : str
        Column to read.
    t : float
        Time in model time units, within the simulated range.
    name : str, optional
        Defaults to ``'At(<selection>, <t>)'``.

    Raises
    ------
    TypeError, ValueError
        For a non-numeric or non-finite ``t``. Evaluating a result whose
        range does not contain ``t`` raises ``ValueError``; the analyses
        (:func:`~nskinetics.analysis.local_sensitivity`,
        :func:`~nskinetics.analysis.morris`, :func:`~nskinetics.analysis.sobol`)
        check ``t`` against the condition's time span before any run.

    Notes
    -----
    Basis: linear interpolation between the two neighbouring output rows
    (exact at an output time), in the selection's units. An event jump
    between those rows is interpolated across; put ``t`` on the grid (e.g.
    with ``Condition(times=...)``) when that matters.
    """
    _fields = ('selection', 't')

    def __init__(self, selection, t, *, name=None):
        self.selection = _selection('selection', selection)
        self.t = _finite('t', t)
        super().__init__(name if name is not None
                         else f'At({selection}, {_fmt(self.t)})',
                         (selection,))

    def _check_condition(self, condition):
        """``ValueError`` if ``t`` is outside the time span of
        ``condition`` (its ``times``, else ``t0 .. t_end``)."""
        span = _time_span(condition)
        if not span[0] <= self.t <= span[1]:
            raise ValueError(f'{self.name}: t={self.t!r} is outside the '
                             f'simulated range [{span[0]!r}, {span[1]!r}] of '
                             'the condition.')

    def _evaluate(self, result):
        time = result.time
        if not time[0] <= self.t <= time[-1]:
            raise ValueError(f'{self.name}: t={self.t!r} is outside the '
                             f'simulated range [{time[0]!r}, {time[-1]!r}].')
        return np.interp(self.t, time, result[self.selection])


class Productivity(Metric):
    """Volumetric productivity: product formed per unit of process time.

    Parameters
    ----------
    product : str
        Product concentration selection, e.g. ``'[P]'``.
    turnaround : float, optional
        Non-productive time added to the run length (cleaning, filling...),
        in model time units. Defaults to ``0.0``.
    name : str, optional
        Defaults to ``'Productivity(<product>)'``.

    Raises
    ------
    TypeError, ValueError
        For an invalid ``product`` selection (see :class:`Metric`), or a
        non-numeric, non-finite or negative ``turnaround``.

    Notes
    -----
    Basis: ``(P_end - P_0) / (t_end - t_0 + turnaround)`` from the first and
    last output rows, in concentration units per model time unit (e.g.
    g/L/h). It uses concentrations, so dilution by fed volume lowers it;
    ``NaN`` if the denominator is zero.
    """
    _fields = ('product', 'turnaround')

    def __init__(self, product, *, turnaround=0.0, name=None):
        self.product = _selection('product', product)
        self.turnaround = _finite('turnaround', turnaround)
        if self.turnaround < 0:
            raise ValueError(f'turnaround must be >= 0; got {turnaround!r}.')
        super().__init__(name if name is not None
                         else f'Productivity({product})', (product,))

    def _evaluate(self, result):
        p, t = result[self.product], result.time
        return _ratio(p[-1] - p[0], t[-1] - t[0] + self.turnaround)


class Yield(Metric):
    """Product formed per substrate consumed or added, on an amount basis.

    Parameters
    ----------
    product : str
        Product concentration selection, e.g. ``'[P]'``.
    substrate : str
        Substrate concentration selection, e.g. ``'[S]'``.
    volume : str, optional
        Volume selection (e.g. a compartment such as ``'env'``) by which
        concentrations are multiplied to give amounts. ``None`` (default)
        means a constant volume: concentrations are used as they are.
    fed : tuple, optional
        ``(fed_volume, feed_conc)`` for substrate added during the run:
        ``fed_volume`` is the selection of the *total* volume fed so far
        (e.g. a :class:`~nskinetics.FeedSpike` ``tot_vol_var``) and
        ``feed_conc`` the substrate concentration of the feed, as a
        selection (read at the last row) or a number. Requires ``volume``.
    basis : {'consumed', 'added'}
        Denominator: substrate consumed (default) or substrate added.
    name : str, optional
        Defaults to ``'Yield(<product>/<substrate>, <basis>)'``.

    Raises
    ------
    TypeError, ValueError
        For invalid selections, an unknown ``basis``, a malformed ``fed``,
        or ``fed`` without ``volume``.

    Notes
    -----
    With amounts ``A_x = [x] * V`` (``V = 1`` without ``volume``) and the fed
    amount ``F = (fed_volume_end - fed_volume_0) * feed_conc``:

    - formed product ``A_P,end - A_P,0``;
    - ``basis='added'``: ``(A_P,end - A_P,0) / (A_S,0 + F)``;
    - ``basis='consumed'``: ``(A_P,end - A_P,0) / (A_S,0 + F - A_S,end)``.

    The value is in product amount per substrate amount in the model's
    concentration basis (g/g for mass concentrations, mol/mol for molar
    ones); no stoichiometric or theoretical-yield factor is applied.
    ``NaN`` if the denominator is zero.
    """
    _fields = ('product', 'substrate', 'volume', 'fed', 'basis')

    def __init__(self, product, substrate, *, volume=None, fed=None,
                 basis='consumed', name=None):
        self.product = _selection('product', product)
        self.substrate = _selection('substrate', substrate)
        self.volume = None if volume is None else _selection('volume', volume)
        if basis not in _YIELD_BASES:
            raise ValueError(f'basis must be one of {_YIELD_BASES}; got '
                             f'{basis!r}.')
        self.basis = basis
        if fed is not None:
            if (isinstance(fed, (str, bytes)) or not hasattr(fed, '__len__')
                    or len(fed) != 2):
                raise TypeError('fed must be a (fed_volume, feed_conc) pair; '
                                f'got {fed!r}.')
            fed_volume, feed_conc = fed
            fed_volume = _selection('fed[0] (fed volume)', fed_volume)
            if isinstance(feed_conc, str):
                feed_conc = _selection('fed[1] (feed concentration)',
                                       feed_conc)
            else:
                feed_conc = _finite('fed[1] (feed concentration)', feed_conc)
            if self.volume is None:
                raise ValueError('fed= needs volume=: fed substrate is an '
                                 'amount (fed volume x feed concentration), '
                                 'so the broth amounts need the volume '
                                 'selection too.')
            fed = (fed_volume, feed_conc)
        self.fed = fed
        sels = [product, substrate]
        if self.volume is not None:
            sels.append(self.volume)
        if fed is not None:
            sels.append(fed[0])
            if isinstance(fed[1], str):
                sels.append(fed[1])
        super().__init__(name if name is not None
                         else f'Yield({product}/{substrate}, {basis})', sels)

    def _evaluate(self, result):
        v = (result[self.volume] if self.volume is not None
             else np.ones(len(result)))
        p = result[self.product] * v
        s = result[self.substrate] * v
        fed_amount = 0.0
        if self.fed is not None:
            fed_volume, feed_conc = self.fed
            fv = result[fed_volume]
            conc = (result[feed_conc][-1] if isinstance(feed_conc, str)
                    else feed_conc)
            fed_amount = (fv[-1] - fv[0]) * conc
        supplied = s[0] + fed_amount
        denominator = supplied if self.basis == 'added' else supplied - s[-1]
        return _ratio(p[-1] - p[0], denominator)


class CallableMetric(Metric):
    """A user function as a metric.

    Parameters
    ----------
    fn : callable
        ``fn(result) -> float`` for a :class:`~nskinetics.SimulationResult`
        (e.g. ``lambda res: res.final['[P]'] / res.final['[X]']``). It is
        evaluated under ``numpy.errstate(all='ignore')``.
    name : str
        Metric name.
    selections : sequence of str, optional
        Columns ``fn`` reads, recorded by batch runs in addition to the
        condition's selections.

    Raises
    ------
    TypeError, ValueError
        If ``fn`` is not callable, or for an invalid ``name`` or selection
        (see :class:`Metric`).

    Notes
    -----
    Named ``CallableMetric`` rather than ``Callable`` so that
    ``from nskinetics.analysis import *`` cannot shadow
    ``typing.Callable`` / ``collections.abc.Callable``;
    :meth:`Metric.from_function` is an equivalent constructor. A function
    that returns a non-scalar raises ``TypeError`` when evaluated. For
    worker processes started with ``spawn`` (macOS, Windows), ``fn`` must
    be a module-level function.
    """
    def __init__(self, fn, name, selections=()):
        if not callable(fn):
            raise TypeError(f'fn must be callable; got {type(fn).__name__}.')
        self.fn = fn
        super().__init__(name, selections)

    _fields = ('fn',)

    def __repr__(self):
        return (f'CallableMetric({self.fn!r}, name={self.name!r}, '
                f'selections={self.selections!r})')

    def _evaluate(self, result):
        return self.fn(result)
