# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Run descriptions and stored trajectories.

:class:`Condition` is the single, immutable description of one simulation run
(initial values, parameter overrides, output grid, reset behaviour) and
:class:`SimulationResult` the single trajectory type that
:meth:`nskinetics.KineticModel.run`, :func:`nskinetics.run_many` and the
analysis metrics exchange.

This module deliberately does not import :mod:`.kinetic_model` (which imports
it), so there is no import cycle.
"""
import json
import math
import numbers
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field, fields
from types import MappingProxyType
from typing import Any

import numpy as np
import pandas as pd

from ._plotting import _plot_time_course

__all__ = ('Condition', 'SimulationResult')


# --- helpers ------------------------------------------------------------------
def _freeze(value):
    """Deep-freeze plain containers: mappings become read-only
    ``MappingProxyType`` views over copies, lists and tuples become tuples."""
    if isinstance(value, Mapping):
        return MappingProxyType({k: _freeze(v) for k, v in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(v) for v in value)
    return value


def _thaw(value):
    """Inverse of :func:`_freeze` with JSON-friendly containers: mappings
    become dicts and tuples become lists (recursively)."""
    if isinstance(value, Mapping):
        return {k: _thaw(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_thaw(v) for v in value]
    return value


def _hashable(value):
    """Hashable stand-in for a frozen value (mappings -> frozensets)."""
    if isinstance(value, Mapping):
        return frozenset((k, _hashable(v)) for k, v in value.items())
    if isinstance(value, tuple):
        return tuple(_hashable(v) for v in value)
    return value


def _is_json_safe(value):
    """Whether ``value`` (thawed) serializes to JSON and reads back equal."""
    if value is None or isinstance(value, (bool, str)):
        return True
    if isinstance(value, numbers.Integral):
        return isinstance(value, int)
    if isinstance(value, float):
        return math.isfinite(value)
    if isinstance(value, list):
        return all(_is_json_safe(v) for v in value)
    if isinstance(value, dict):
        return all(isinstance(k, str) and _is_json_safe(v)
                   for k, v in value.items())
    return False


def _real(name, value):
    """Return ``value`` as a finite float, or raise TypeError/ValueError."""
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise TypeError(f'{name} must be a real number; got '
                        f'{type(value).__name__} {value!r}.')
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f'{name} must be finite; got {value!r}.')
    return value


def _value_map(name, mapping):
    """Validate a ``{selection: value}`` mapping; return a plain dict of
    finite floats."""
    if not isinstance(mapping, Mapping):
        raise TypeError(f'{name} must be a mapping of model ids to numbers; '
                        f'got {type(mapping).__name__}.')
    out = {}
    for key, value in mapping.items():
        if not isinstance(key, str) or not key:
            raise TypeError(f'{name} keys must be non-empty strings; got '
                            f'{key!r}.')
        out[key] = _real(f'{name}[{key!r}]', value)
    return out


def _time_grid(times):
    """Validate an explicit output grid; return it as a tuple of floats."""
    if isinstance(times, (str, bytes)) or isinstance(times, Mapping):
        raise TypeError('times must be a sequence of numbers.')
    try:
        grid = np.asarray(times, dtype=float)
    except (TypeError, ValueError) as e:
        raise TypeError(f'times must be a sequence of numbers: {e}') from e
    if grid.ndim != 1 or grid.size < 2:
        raise ValueError('times must be a 1-D sequence of at least two '
                         'output times.')
    if not np.all(np.isfinite(grid)):
        raise ValueError('times must be finite.')
    if not np.all(np.diff(grid) > 0):
        raise ValueError('times must be strictly increasing.')
    return tuple(float(t) for t in grid)


# --- Condition ----------------------------------------------------------------
@dataclass(frozen=True, kw_only=True, repr=False)
class Condition:
    """Immutable description of one simulation run.

    A condition says how to start, perturb and record a run of a
    :class:`~nskinetics.KineticModel`: :meth:`KineticModel.run` resets the
    model (if ``reset``), applies ``parameters``, then ``initial``, integrates
    over ``[t0, t_end]`` (or the explicit ``times`` grid) and returns a
    :class:`SimulationResult`. The same object describes an experiment in
    calibration, a perturbation in a sensitivity analysis and a sample in a
    screening or scan (see :func:`nskinetics.run_many`).

    All arguments are keyword-only. Every field is validated on construction
    and the object cannot be modified afterwards: mappings are stored as
    read-only views over copies and sequences as tuples. Use
    :meth:`replace`, :meth:`with_parameters` and :meth:`with_initial` to
    derive new conditions.

    Parameters
    ----------
    initial : mapping of str to float, optional
        Values set *after* the reset and the parameters, i.e. the initial
        state of this run, as RoadRunner selections in the model's own units:
        ``{'[S]': 50.0}`` (concentration), ``'S'`` (amount), compartments
        (``'env'``) or any other settable id. The order of the keys does not
        matter: ``init(x)`` selections are set first, then compartment
        volumes, then the other entries, so every value you give holds at
        ``t0``. A volume alone keeps species *amounts*, so ``{'env': 2}``
        halves the starting concentrations of a 1 L baseline, while
        ``{'env': 2, '[S]': 50}`` starts 2 L at ``[S] = 50``. A species may
        be given as an amount or a concentration, not both. ``init(x)``
        entries re-initialize the model; see Notes. Defaults to ``{}``.
    parameters : mapping of str to float, optional
        Values set right after the reset (before ``initial``) and restored to
        their pre-run values when the run ends, e.g. ``{'k': 0.4}``. They are
        applied in the same fixed order as ``initial`` (``init()`` keys, then
        compartments, then the rest), so their key order does not matter,
        and a species may be given as an amount or a concentration, not
        both. Plain ids do not re-evaluate initial assignments (see Notes).
        Defaults to ``{}``. An id may not appear in both ``initial`` and
        ``parameters``.
    t0 : float, optional
        Start time, in model time units. ``None`` (default) means ``0.0``
        with ``t_end`` and ``times[0]`` with ``times``; it is stored as that
        float. Integration with ``times`` starts at ``times[0]``, so an
        explicit ``t0`` must then equal ``times[0]`` (include ``t0`` in
        ``times`` to start earlier).
    t_end : float, optional
        End time, in model time units; must exceed ``t0``. Exactly one of
        ``t_end`` and ``times`` must be given.
    n_steps : int
        Number of output points over ``[t0, t_end]`` (at least 2). Defaults
        to ``1000``. Ignored with ``times``.
    times : sequence of float, optional
        Explicit output grid in model time units: at least two finite,
        strictly increasing points. Rows are reported exactly at these times
        (events still fire in between).
    selections : sequence of str, optional
        Columns to record. ``'time'`` is always recorded (first, if not
        listed). ``None`` (default) records ``'time'`` and the concentration
        of every floating species.
    reset : bool
        Whether to reset the model before applying the condition. Defaults to
        ``True``. :func:`nskinetics.run_many` requires ``True``.
    reset_kwargs : mapping, optional
        Keyword arguments for :meth:`KineticModel.reset`: a ``mode``, a
        ``preserve`` sequence, or arguments for a custom reset function such
        as ``{'reset_spike_cap': True}``. Nested lists/dicts are frozen to
        tuples/read-only mappings. Defaults to ``{}``.
    label : str, optional
        Name of the run, carried by its :class:`SimulationResult` and the
        ``label`` column of batch tables.

    Raises
    ------
    TypeError
        If a field has the wrong type (e.g. a non-numeric value, a non-integer
        ``n_steps`` or a non-bool ``reset``).
    ValueError
        If neither or both of ``t_end`` and ``times`` are given, ``times`` is
        not strictly increasing, an explicit ``t0`` differs from
        ``times[0]``, ``t_end <= t0``, ``n_steps < 2``, a value is not
        finite, ``selections`` has duplicates, an id is listed in both
        ``initial`` and ``parameters``, or ``initial`` or ``parameters`` sets
        a species both as an amount and as a concentration (``'S'`` and
        ``'[S]'``, or ``'init(S)'`` and ``'init([S])'``).

    Notes
    -----
    Model ids are checked against a model only when the condition is run
    (:meth:`KineticModel.run` raises ``ValueError`` listing unknown ids
    before changing any state).

    ``parameters`` (plain ids) are applied after the reset, so they do not
    re-evaluate SBML initial assignments (``S = 20*k`` keeps the ``S``
    computed from the model's own ``k``); set such initial values through
    ``initial``, or give ``init(k)`` to run the model as if it were defined
    with that value.

    ``init(x)`` keys, in either mapping, are written before the reset,
    because RoadRunner re-initializes the model on every ``init(...)``
    write; ``{'init(k)': v}`` runs the model as if it were defined with
    ``k = v``. ``init()`` of a species must use the form its SBML declares:
    ``init([S])`` for the usual Antimony ``S = 50`` (an initial
    concentration), ``init(S)`` for one declared with an initial amount.
    :meth:`KineticModel.run` rejects the other form, and ``init()`` of an
    id defined by an initial assignment (give ``init()`` of its inputs,
    e.g. ``init(k)``), because RoadRunner would rewrite the model for good.
    See the Notes of :meth:`KineticModel.run` and the concepts page
    ("``Condition.parameters`` vs ``Condition.initial``").

    Examples
    --------
    >>> cond = Condition(initial={'[S]': 50}, parameters={'k': 0.4},
    ...                  t_end=48, n_steps=481)
    >>> cond.with_parameters(k=0.8).parameters['k']
    0.8
    """
    initial: Mapping[str, float] = field(default_factory=dict)
    parameters: Mapping[str, float] = field(default_factory=dict)
    t0: float | None = None
    t_end: float | None = None
    n_steps: int = 1000
    times: Sequence[float] | None = None
    selections: Sequence[str] | None = None
    reset: bool = True
    reset_kwargs: Mapping[str, Any] = field(default_factory=dict)
    label: str | None = None

    def __post_init__(self):
        put = object.__setattr__
        initial = _value_map('initial', self.initial)
        parameters = _value_map('parameters', self.parameters)
        both = [k for k in parameters if k in initial]
        if both:
            raise ValueError(f'id(s) {both} appear in both initial and '
                             'parameters; give each id once.')
        # The only conflicts the fixed apply order of run() cannot settle:
        # 'S' with '[S]', and 'init(S)' with 'init([S])'.
        for name, values in (('initial', initial),
                             ('parameters', parameters)):
            twice = [k for k in values if f'[{k}]' in values]
            twice += [k for k in values
                      if k.startswith('init(') and k.endswith(')')
                      and f'init([{k[5:-1]}])' in values]
            if twice:
                raise ValueError(f'{name} sets {twice} both as an amount and '
                                 'as a concentration ([...]); give one of '
                                 'them.')
        t0 = None if self.t0 is None else _real('t0', self.t0)
        n_steps = self.n_steps
        if isinstance(n_steps, bool) or not isinstance(n_steps,
                                                       numbers.Integral):
            raise TypeError(f'n_steps must be an integer; got {n_steps!r}.')
        n_steps = int(n_steps)
        if n_steps < 2:
            raise ValueError(f'n_steps must be at least 2; got {n_steps}.')
        if (self.t_end is None) == (self.times is None):
            raise ValueError('Give exactly one of t_end and times.')
        t_end, times = None, None
        if self.times is not None:
            times = _time_grid(self.times)
            if t0 is None:
                t0 = times[0]
            elif t0 != times[0]:
                raise ValueError(f't0={t0!r} conflicts with times[0]='
                                 f'{times[0]!r}; integration starts at '
                                 'times[0] (include t0 in times to start '
                                 'there).')
        else:
            if t0 is None:
                t0 = 0.0
            t_end = _real('t_end', self.t_end)
            if t_end <= t0:
                raise ValueError(f't_end={t_end!r} must exceed t0={t0!r}.')
        selections = self.selections
        if selections is not None:
            if isinstance(selections, (str, bytes, Mapping)) or not isinstance(
                    selections, (Sequence, np.ndarray)):
                raise TypeError('selections must be a sequence of strings.')
            selections = tuple(selections)
            if not selections or not all(isinstance(s, str) and s
                                         for s in selections):
                raise TypeError('selections must be non-empty strings.')
            # Plain str (not np.str_), so to_dict, equality and hashing agree.
            selections = tuple(str(s) for s in selections)
            dupes = sorted({s for s in selections if selections.count(s) > 1})
            if dupes:
                raise ValueError(f'Duplicate selections {dupes}.')
        reset = self.reset
        if not isinstance(reset, (bool, np.bool_)):
            raise TypeError(f'reset must be a bool; got {reset!r}.')
        if not isinstance(self.reset_kwargs, Mapping):
            raise TypeError('reset_kwargs must be a mapping.')
        bad = [k for k in self.reset_kwargs if not isinstance(k, str)]
        if bad:
            raise TypeError(f'reset_kwargs keys must be strings; got {bad}.')
        if self.label is not None and not isinstance(self.label, str):
            raise TypeError(f'label must be a string or None; got '
                            f'{self.label!r}.')
        put(self, 'initial', MappingProxyType(initial))
        put(self, 'parameters', MappingProxyType(parameters))
        put(self, 't0', t0)
        put(self, 't_end', t_end)
        put(self, 'n_steps', n_steps)
        put(self, 'times', times)
        put(self, 'selections', selections)
        put(self, 'reset', bool(reset))
        put(self, 'reset_kwargs', _freeze(dict(self.reset_kwargs)))

    # --- derivation ---------------------------------------------------------
    def replace(self, **changes):
        """Return a copy with the given fields replaced, then re-validated.

        Fields are replaced wholesale (``dataclasses.replace`` semantics); use
        :meth:`with_parameters` / :meth:`with_initial` to *merge* new values
        into the existing mappings.

        Parameters
        ----------
        **changes
            New field values, e.g. ``times=[0, 1, 2]`` or ``label='FB2'``.

        Returns
        -------
        Condition

        Raises
        ------
        TypeError
            If a keyword is not a field of :class:`Condition`.
        ValueError
            If the resulting condition is invalid.

        Notes
        -----
        ``t_end`` and ``times`` are mutually exclusive, so setting one of them
        to a value clears the other unless it is also passed. Setting
        ``times`` without ``t0`` makes the run start at the new ``times[0]``.
        Setting ``t_end`` keeps the current ``t0``; for a condition built
        from ``times``, that is ``times[0]``. Pass ``t0`` as well (e.g.
        ``replace(t_end=10, t0=0.0)``) to start elsewhere.
        """
        names = [f.name for f in fields(self)]
        unknown = [k for k in changes if k not in names]
        if unknown:
            raise TypeError(f'Condition has no field(s) {unknown}; expected '
                            f'some of {names}.')
        if changes.get('times') is not None:
            changes.setdefault('t_end', None)
            changes.setdefault('t0', None)
        elif changes.get('t_end') is not None:
            changes.setdefault('times', None)
        kwargs = {name: getattr(self, name) for name in names}
        kwargs.update(changes)
        return type(self)(**kwargs)

    def with_parameters(self, mapping=None, **values):
        """Return a copy whose ``parameters`` are updated with new entries.

        Parameters
        ----------
        mapping : mapping of str to float, optional
            ``{id: value}`` entries to add or override.
        **values : float
            More ``id=value`` entries (plain ids); applied after ``mapping``.

        Returns
        -------
        Condition
            The existing parameters, updated with the new entries; all other
            fields unchanged.

        Raises
        ------
        TypeError, ValueError
            As for the constructor (e.g. an id that is also in ``initial``).
        """
        merged = dict(self.parameters)
        merged.update(mapping or {})
        merged.update(values)
        return self.replace(parameters=merged)

    def with_initial(self, mapping=None, **values):
        """Return a copy whose ``initial`` values are updated with new entries.

        Parameters
        ----------
        mapping : mapping of str to float, optional
            ``{selection: value}`` entries to add or override, e.g.
            ``{'[S]': 30.0}``.
        **values : float
            More ``id=value`` entries (plain ids); applied after ``mapping``.

        Returns
        -------
        Condition
            The existing initial values, updated with the new entries; all
            other fields unchanged.

        Raises
        ------
        TypeError, ValueError
            As for the constructor (e.g. an id that is also in
            ``parameters``).
        """
        merged = dict(self.initial)
        merged.update(mapping or {})
        merged.update(values)
        return self.replace(initial=merged)

    # --- serialization ------------------------------------------------------
    def to_dict(self):
        """Return a JSON-safe ``dict`` of every field.

        Mappings become dicts and sequences lists, so
        ``Condition.from_dict(json.loads(json.dumps(c.to_dict()))) == c``.

        Returns
        -------
        dict

        Raises
        ------
        TypeError
            If ``reset_kwargs`` holds a value JSON cannot represent exactly
            (anything but None, bool, int, finite float, str, list, dict).
        """
        out = {f.name: _thaw(getattr(self, f.name)) for f in fields(self)}
        if not _is_json_safe(out['reset_kwargs']):
            raise TypeError('reset_kwargs holds values that are not '
                            'JSON-safe; cannot serialize the condition: '
                            f'{out["reset_kwargs"]!r}.')
        return out

    @classmethod
    def from_dict(cls, d):
        """Build a condition from a mapping such as :meth:`to_dict` returns.

        Parameters
        ----------
        d : mapping
            Field names to values; missing fields take their defaults.

        Returns
        -------
        Condition

        Raises
        ------
        ValueError
            If ``d`` has keys that are not fields, or the condition is
            invalid.
        TypeError
            If ``d`` is not a mapping, or a field has the wrong type.
        """
        if not isinstance(d, Mapping):
            raise TypeError(f'from_dict expects a mapping; got '
                            f'{type(d).__name__}.')
        names = {f.name for f in fields(cls)}
        unknown = sorted(k for k in d if k not in names)
        if unknown:
            raise ValueError(f'Unknown Condition field(s) {unknown}.')
        return cls(**d)

    # --- dunder -------------------------------------------------------------
    def __repr__(self):
        body = ', '.join(f'{f.name}={_thaw_repr(getattr(self, f.name))}'
                         for f in fields(self))
        return f'{type(self).__name__}({body})'

    def __hash__(self):
        return hash(tuple(_hashable(getattr(self, f.name))
                          for f in fields(self)))

    def __getstate__(self):
        # MappingProxyType cannot be pickled: ship plain containers and
        # re-validate/re-freeze on the way back in.
        return {f.name: _thaw(getattr(self, f.name)) for f in fields(self)}

    def __setstate__(self, state):
        self.__init__(**state)


def _thaw_repr(value):
    """repr of a frozen field with mappings shown as dicts."""
    if isinstance(value, Mapping):
        return repr({k: _thaw(v) for k, v in value.items()})
    return repr(value)


# --- SimulationResult ---------------------------------------------------------
_CSV_FORMAT = 'SimulationResult'


@dataclass(eq=False, repr=False)
class SimulationResult:
    """One simulated trajectory with its units and the condition that made it.

    Returned by :meth:`nskinetics.KineticModel.run` and collected by
    :func:`nskinetics.run_many`; the analysis metrics read it.

    Parameters
    ----------
    df : pandas.DataFrame
        The trajectory, one row per output time; must have a ``'time'``
        column (model time units). Other columns are RoadRunner selections in
        the model's own units (``'[S]'`` concentration, ``'env'`` volume...).
    units : dict, optional
        The model's units, e.g. ``{'time': 'h', 'conc': 'g/L'}`` (copied).
    condition : Condition, optional
        The condition that produced the run.
    label : str, optional
        Name of the run. Defaults to ``condition.label``.

    Raises
    ------
    TypeError
        If ``df`` is not a DataFrame or ``condition`` not a :class:`Condition`.
    ValueError
        If ``df`` has no ``'time'`` column.

    Notes
    -----
    Indexing returns numpy arrays: ``res['[P]']`` is the ``[P]`` column, so
    ``res['[P]'][-1]`` is its final value (as is ``res.final['[P]']``).
    """
    df: pd.DataFrame
    units: dict = field(default_factory=dict)
    condition: Condition | None = None
    label: str | None = None

    def __post_init__(self):
        if not isinstance(self.df, pd.DataFrame):
            raise TypeError(f'df must be a pandas DataFrame; got '
                            f'{type(self.df).__name__}.')
        if 'time' not in self.df.columns:
            raise ValueError("A SimulationResult needs a 'time' column.")
        if self.condition is not None and not isinstance(self.condition,
                                                         Condition):
            raise TypeError(f'condition must be a Condition; got '
                            f'{type(self.condition).__name__}.')
        self.units = dict(self.units) if self.units is not None else {}
        if self.label is None and self.condition is not None:
            self.label = self.condition.label

    # --- accessors ----------------------------------------------------------
    @property
    def time(self):
        """Output times (model time units), as a new, writable 1-D
        ``numpy.ndarray``."""
        return self.df['time'].to_numpy(dtype=float, copy=True)

    @property
    def array(self):
        """The trajectory as a new, writable 2-D ``numpy.ndarray`` (rows =
        times, columns = :attr:`col_names`)."""
        return self.df.to_numpy(dtype=float, copy=True)

    @property
    def col_names(self):
        """Column names (selections), in order, as a list."""
        return list(self.df.columns)

    @property
    def final(self):
        """Values at the last output time, as a ``pandas.Series`` indexed by
        column name."""
        return self.df.iloc[-1].copy()

    def __getitem__(self, key):
        """Column ``key`` as a 1-D array, or columns ``[k1, k2, ...]`` as a
        2-D array (new, writable arrays that do not alias :attr:`df`).

        Raises
        ------
        KeyError
            If a column is not recorded, listing the recorded ones.
        TypeError
            If ``key`` is neither a column name nor a list of names (e.g. a
            row number; use ``res.df.iloc[i]`` for rows).
        """
        if not isinstance(key, str) and (isinstance(key, (bytes, Mapping))
                                         or not isinstance(key, Iterable)):
            raise TypeError(
                "index a SimulationResult with a column name or a list of "
                "names (e.g. res['[S]']); use res.df.iloc[i] for rows; got "
                f'{type(key).__name__}.')
        keys = [key] if isinstance(key, str) else list(key)
        missing = [k for k in keys if k not in self.df.columns]
        if missing:
            raise KeyError(f'Column(s) {missing} not recorded; available: '
                           f'{self.col_names}.')
        if isinstance(key, str):
            return self.df[key].to_numpy(dtype=float, copy=True)
        return self.df[keys].to_numpy(dtype=float, copy=True)

    def __contains__(self, key):
        return key in self.df.columns

    def __len__(self):
        return len(self.df)

    def __repr__(self):
        t = self.time
        return (f'{type(self).__name__}(label={self.label!r}, rows={len(t)}, '
                f't=[{t[0]!r}, {t[-1]!r}], columns={self.col_names}, '
                f'units={self.units})')

    def at(self, t):
        """Linearly interpolate every column at time(s) ``t``.

        Parameters
        ----------
        t : float or sequence of float
            Time(s) in model time units, within ``[time[0], time[-1]]``.

        Returns
        -------
        pandas.Series or pandas.DataFrame
            For a scalar ``t``, a Series indexed by column name (its
            ``'time'`` entry is ``t``); for a sequence, a DataFrame with one
            row per time. Event-driven jumps between output points are
            interpolated across, so request a grid that resolves them (or
            use ``times=``) when that matters.

        Raises
        ------
        ValueError
            If any ``t`` is not finite or lies outside the simulated range.
        """
        time = self.time
        tt = np.asarray(t, dtype=float)
        if tt.ndim > 1:
            raise ValueError('at() takes a scalar or a 1-D sequence of times.')
        flat = np.atleast_1d(tt)
        if (not np.all(np.isfinite(flat)) or np.any(flat < time[0])
                or np.any(flat > time[-1])):
            raise ValueError(f'at({t!r}) is outside the simulated range '
                             f'[{time[0]!r}, {time[-1]!r}].')
        cols = self.col_names
        values = np.column_stack(
            [flat if c == 'time' else
             np.interp(flat, time, self.df[c].to_numpy(dtype=float))
             for c in cols])
        if tt.ndim == 0:
            return pd.Series(values[0], index=cols, name=float(tt))
        return pd.DataFrame(values, columns=cols)

    # --- IO -----------------------------------------------------------------
    def to_csv(self, path):
        """Write the trajectory to a wide CSV with a JSON metadata header.

        The file starts with ``# key: <json>`` lines (``nskinetics``,
        ``units``, ``label`` and ``condition``, the latter from
        :meth:`Condition.to_dict`), followed by the trajectory with one column
        per selection. Floats are written with their shortest round-trip
        representation, so :meth:`from_csv` restores the result exactly.

        Parameters
        ----------
        path : str or path-like
            Destination; an existing file is overwritten.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If the condition (its ``reset_kwargs``) or units are not
            JSON-serializable.
        """
        meta = {'nskinetics': _CSV_FORMAT,
                'units': self.units,
                'label': self.label,
                'condition': (None if self.condition is None
                              else self.condition.to_dict())}
        with open(path, 'w', encoding='utf-8', newline='') as f:
            for key, value in meta.items():
                f.write(f'# {key}: {json.dumps(value, allow_nan=False)}\n')
            self.df.to_csv(f, index=False, lineterminator='\n')

    @classmethod
    def from_csv(cls, path):
        """Read a result written by :meth:`to_csv`.

        Parameters
        ----------
        path : str or path-like
            CSV written by :meth:`to_csv`.

        Returns
        -------
        SimulationResult
            Equivalent to the one written: same columns, exact float values,
            units, label and condition (:class:`SimulationResult` does not
            define ``==``; compare e.g. ``back.df.equals(res.df)``).

        Raises
        ------
        ValueError
            If the file lacks the metadata header.
        """
        meta, n_header = {}, 0
        with open(path, encoding='utf-8', newline='') as f:
            for line in f:
                if not line.startswith('#'):
                    break
                n_header += 1
                key, sep, text = line[1:].strip().partition(': ')
                if not sep:
                    raise ValueError(f'Malformed header line {line!r} in '
                                     f'{path!s}.')
                meta[key] = json.loads(text)
        if meta.get('nskinetics') != _CSV_FORMAT:
            raise ValueError(f'{path!s} was not written by '
                             'SimulationResult.to_csv().')
        df = pd.read_csv(path, skiprows=n_header, encoding='utf-8',
                         float_precision='round_trip')
        cond = meta.get('condition')
        return cls(df=df, units=meta.get('units') or {},
                   condition=None if cond is None else Condition.from_dict(cond),
                   label=meta.get('label'))

    # --- plotting -----------------------------------------------------------
    def plot(self, variables=None, labels=None, xlim=None, ylim=None,
             markers=None, save_fig=False, filename=None, figwidth=3.9,
             ylabel=None, ylabel_unit=None, **kwargs):
        """Plot columns against time with the package's standard styling.

        Same arguments and output as
        :meth:`nskinetics.KineticModel.plot_simulation_results`, for this
        stored trajectory: by default every concentration column (names with
        ``'['`` and ``']'``) is plotted, with axis labels in :attr:`units`.

        Parameters
        ----------
        variables, labels, xlim, ylim, markers, save_fig, filename, figwidth, ylabel, ylabel_unit, **kwargs
            See :meth:`nskinetics.KineticModel.plot_simulation_results`.

        Returns
        -------
        (matplotlib.figure.Figure, matplotlib.axes.Axes)

        Raises
        ------
        KineticSimulationError
            If there is no concentration column and ``variables`` is not
            given.
        """
        return _plot_time_course(
            self.df, self.col_names, self.units, variables=variables,
            labels=labels, xlim=xlim, ylim=ylim, markers=markers,
            save_fig=save_fig, filename=filename, figwidth=figwidth,
            ylabel=ylabel, ylabel_unit=ylabel_unit, **kwargs)
