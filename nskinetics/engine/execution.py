# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Batch execution of :class:`~nskinetics.Condition` objects.

:func:`run_many` and :func:`scan` are the package's only parallel executor;
sensitivity screening, calibration multistarts and prediction bands all run
through them. Every batch runs on a private model -- a :meth:`KineticModel.copy`
in the calling process, or one model per worker process built from a picklable
blueprint -- so the caller's model is never touched.
"""
import functools
import itertools
import json
import math
import multiprocessing
import numbers
import pickle
import reprlib
import sys
from collections.abc import Mapping
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt

from .conditions import Condition
from .kinetic_model import KineticModel
from ._plotting import _set_font, _style_axes

__all__ = ('RunBatch', 'ScanResult', 'run_many', 'scan')

# Columns every batch table carries ahead of the metric columns.
_RESERVED_COLUMNS = ('label', 'status', 'error')
_ON_ERROR = ('record', 'raise')
_SCAN_MODES = ('product', 'zip')
_PLACEMENTS = ('initial', 'parameters')
#: Columns of the long-format CSV written by :meth:`RunBatch.to_csv`.
_CSV_COLUMNS = ('row', 'quantity', 'key', 'value')


# --- results ------------------------------------------------------------------
@dataclass(eq=False, repr=False)
class RunBatch:
    """Outcome of :func:`run_many`: one row of metrics per condition.

    Parameters
    ----------
    metrics : pandas.DataFrame
        One row per condition, in order, with the columns ``label`` (the
        condition's label or ``None``), ``status`` (``'ok'`` or
        ``'failed'``), ``error`` (``None``, or the error message of a failed
        run or metric) and one float column per metric (``NaN`` where the run
        or that metric failed).
    results : list of SimulationResult or None
        Trajectories aligned with ``conditions``; ``None`` where the
        integration failed (a row whose run succeeded but a metric failed
        keeps its result), and for every run when ``keep_results=False``.
    conditions : list of Condition
        The conditions, in the order run.
    """
    metrics: pd.DataFrame
    results: list
    conditions: list

    def __len__(self):
        return len(self.conditions)

    def __repr__(self):
        return (f'{type(self).__name__}(n={len(self)}, n_failed='
                f'{self.n_failed}, metrics={self.metric_names})')

    @property
    def metric_names(self):
        """Names of the metric columns, in order (list of str)."""
        return [c for c in self.metrics.columns
                if c not in _RESERVED_COLUMNS and c not in self._input_columns()]

    def _input_columns(self):
        return ()

    @property
    def ok(self):
        """Boolean ``numpy.ndarray``: ``True`` where the run and every metric
        succeeded."""
        return (self.metrics['status'] == 'ok').to_numpy(dtype=bool)

    @property
    def n_failed(self):
        """Number of failed rows (int)."""
        return int(len(self) - self.ok.sum())

    @property
    def inputs(self):
        """The condition values as a table (``pandas.DataFrame``).

        One row per condition and one column per id set in any condition's
        ``parameters`` or ``initial`` (first-seen order), in the model's own
        units; ``NaN`` where a condition does not set that id.
        """
        rows = [{**c.parameters, **c.initial} for c in self.conditions]
        cols = list(dict.fromkeys(k for row in rows for k in row))
        return pd.DataFrame(rows, columns=cols, dtype=float)

    # --- IO -----------------------------------------------------------------
    def _csv_meta(self):
        """Extra ``meta`` rows a subclass writes (JSON-safe values)."""
        return {}

    def _quantity(self, column):
        """The ``quantity`` of a metrics-table column in the CSV."""
        if column in _RESERVED_COLUMNS:
            return column
        return 'input' if column in self._input_columns() else 'metric'

    def to_csv(self, path):
        """Write the batch to a long-format CSV.

        The file has the columns ``row, quantity, key, value``. Each batch row
        ``i`` contributes one CSV row per metrics-table column (``quantity``
        ``label``, ``status``, ``error``, ``input`` for a grid variable or
        ``metric``; ``key`` is the column name) and one ``condition`` row (the
        :meth:`Condition.to_dict` of that run). Rows with ``quantity`` ``meta``
        and an empty ``row`` record the class, the row count, the column order
        and, for a :class:`ScanResult`, the grid, mode and placement. Every
        ``value`` is JSON (``NaN`` and ``Infinity`` for non-finite floats),
        so :meth:`from_csv` restores the metrics table, the conditions and
        the grid exactly.

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
            If a condition's ``reset_kwargs`` are not JSON-safe.

        Notes
        -----
        Trajectories (:attr:`results`) are not written; they read back as
        ``None``. Save the ones you need with :meth:`SimulationResult.to_csv`.
        """
        columns = list(self.metrics.columns)
        meta = {'nskinetics': type(self).__name__, 'n': len(self),
                'columns': columns, **self._csv_meta()}
        rows = [('', 'meta', key, json.dumps(value, allow_nan=False))
                for key, value in meta.items()]
        values = {c: self.metrics[c].tolist() for c in columns}
        quantity = {c: self._quantity(c) for c in columns}
        for i, condition in enumerate(self.conditions):
            for c in columns:
                rows.append((i, quantity[c], c, json.dumps(values[c][i])))
            rows.append((i, 'condition', '',
                         json.dumps(condition.to_dict(), allow_nan=False)))
        pd.DataFrame(rows, columns=list(_CSV_COLUMNS)).to_csv(
            path, index=False, lineterminator='\n', encoding='utf-8')

    @classmethod
    def from_csv(cls, path):
        """Read a batch written by :meth:`to_csv`.

        Parameters
        ----------
        path : str or path-like
            CSV written by :meth:`RunBatch.to_csv` (also for a
            :class:`ScanResult`).

        Returns
        -------
        RunBatch or ScanResult
            The class that wrote the file, with the same metrics table
            (columns, order, dtypes and exact values), conditions and, for a
            scan, grid, mode and placement. :attr:`results` are ``None``.

        Raises
        ------
        ValueError
            If the file was not written by :meth:`to_csv`, or holds a
            :class:`RunBatch` and this is called on :class:`ScanResult`.
        """
        with open(path, encoding='utf-8', newline='') as f:
            header = f.readline().rstrip('\r\n')
        if header != ','.join(_CSV_COLUMNS):
            raise ValueError(f'{path!s} was not written by RunBatch.to_csv().')
        df = pd.read_csv(path, dtype=str, keep_default_na=False,
                         encoding='utf-8')
        meta, cells, conditions = {}, {}, {}
        for row, quantity, key, value in df.itertuples(index=False):
            if quantity == 'meta':
                meta[key] = json.loads(value)
            elif quantity == 'condition':
                conditions[int(row)] = Condition.from_dict(json.loads(value))
            else:
                cells[int(row), key] = json.loads(value)
        kinds = {c.__name__: c for c in (RunBatch, ScanResult)}
        target = kinds.get(meta.get('nskinetics'))
        if target is None:
            raise ValueError(f'{path!s} was not written by RunBatch.to_csv().')
        if not issubclass(target, cls):
            raise ValueError(f'{path!s} holds a {target.__name__}, not a '
                             f'{cls.__name__}.')
        n, columns = meta['n'], meta['columns']
        data = {}
        for c in columns:
            column = [cells[i, c] for i in range(n)]
            data[c] = pd.Series(column, dtype=object if c in _RESERVED_COLUMNS
                                else float)
        return target._from_csv_parts(
            metrics=pd.DataFrame(data, columns=columns),
            results=[None] * n,
            conditions=[conditions[i] for i in range(n)], meta=meta)

    @classmethod
    def _from_csv_parts(cls, metrics, results, conditions, meta):
        return cls(metrics=metrics, results=results, conditions=conditions)


@dataclass(eq=False, repr=False)
class ScanResult(RunBatch):
    """Outcome of :func:`scan`: a :class:`RunBatch` over a grid of values.

    Parameters
    ----------
    metrics : pandas.DataFrame
        As for :class:`RunBatch`, plus one column per grid variable (after
        ``label``) holding that row's value.
    results, conditions
        As for :class:`RunBatch`.
    grid : dict
        ``{variable: tuple of values}`` as scanned.
    mode : {'product', 'zip'}
        How the grid was combined.
    placement : dict
        ``{variable: 'initial' | 'parameters'}``: where each variable was put
        in the generated conditions.
    """
    grid: dict = field(default_factory=dict)
    mode: str = 'product'
    placement: dict = field(default_factory=dict)

    @property
    def variables(self):
        """Grid variable names, in grid order (list of str)."""
        return list(self.grid)

    def _input_columns(self):
        return tuple(self.grid)

    def _csv_meta(self):
        return {'grid': {k: list(v) for k, v in self.grid.items()},
                'mode': self.mode, 'placement': dict(self.placement)}

    @classmethod
    def _from_csv_parts(cls, metrics, results, conditions, meta):
        return cls(metrics=metrics, results=results, conditions=conditions,
                   grid={k: tuple(v) for k, v in meta['grid'].items()},
                   mode=meta['mode'], placement=dict(meta['placement']))

    def pivot(self, metric):
        """Arrange one metric over the grid.

        Parameters
        ----------
        metric : str
            A metric column name.

        Returns
        -------
        pandas.Series or pandas.DataFrame
            For a 1-variable grid, a Series indexed by that variable; for a
            2-variable grid, a DataFrame indexed by the first variable with
            one column per value of the second (sorted ascending; ``NaN`` for
            failed runs and, in ``'zip'`` mode, for combinations not run).

        Raises
        ------
        KeyError
            If ``metric`` is not a metric column.
        ValueError
            If the grid has more than two variables, or (a ``'zip'`` grid) a
            pair of values of a 2-variable grid repeats, so a cell would hold
            several runs; aggregate ``metrics`` instead (e.g.
            ``.metrics.groupby([var0, var1])[metric].mean()``).
        """
        if metric not in self.metric_names:
            raise KeyError(f'Unknown metric {metric!r}; available: '
                           f'{self.metric_names}.')
        v = self.variables
        if len(v) == 1:
            return self.metrics.set_index(v[0])[metric].sort_index()
        if len(v) == 2:
            pairs = self.metrics[v]
            repeated = pairs[pairs.duplicated(keep=False)]
            if len(repeated):
                pairs = sorted(set(repeated.itertuples(index=False,
                                                       name=None)))
                raise ValueError(
                    f'pivot() needs each ({v[0]}, {v[1]}) pair once; '
                    f'{pairs} repeat. Aggregate .metrics instead, e.g. '
                    f'.metrics.groupby({v})[{metric!r}].mean().')
            return self.metrics.pivot(index=v[0], columns=v[1], values=metric)
        raise ValueError(f'pivot() needs a 1- or 2-variable grid; this one has '
                         f'{len(v)} ({v}). Select rows of .metrics instead.')

    def contourf(self, metric, ax=None, levels=10, colorbar=True, **kwargs):
        """Filled contour plot of one metric over a 2-variable product grid.

        Parameters
        ----------
        metric : str
            A metric column name.
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure is created by default.
        levels : int or sequence of float
            Contour levels, passed to ``Axes.contourf``. Defaults to ``10``.
        colorbar : bool
            Whether to add a colorbar labelled with ``metric``. Defaults to
            ``True``.
        **kwargs
            Passed to ``Axes.contourf`` (e.g. ``cmap``).

        Returns
        -------
        (matplotlib.figure.Figure, matplotlib.axes.Axes)
            The first grid variable is on the x-axis and the second on the
            y-axis, in the model's own units. Failed runs (``NaN``) are left
            blank.

        Raises
        ------
        KeyError
            If ``metric`` is not a metric column.
        ValueError
            If the grid is not a 2-variable ``'product'`` grid with at least
            two values per variable.
        """
        v = self.variables
        if len(v) != 2 or self.mode != 'product':
            raise ValueError('contourf() needs a 2-variable product grid.')
        table = self.pivot(metric)
        if table.shape[0] < 2 or table.shape[1] < 2:
            raise ValueError('contourf() needs at least two values per grid '
                             'variable.')
        x = table.index.to_numpy(dtype=float)
        y = table.columns.to_numpy(dtype=float)
        z = np.ma.masked_invalid(table.to_numpy(dtype=float).T)
        _set_font()
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        cs = ax.contourf(x, y, z, levels=levels, **kwargs)
        if colorbar:
            fig.colorbar(cs, ax=ax).set_label(metric)
        ax.set_xlabel(v[0])
        ax.set_ylabel(v[1])
        _style_axes(ax)
        return fig, ax


# --- helpers ------------------------------------------------------------------
def _normalize_metrics(metrics):
    """Return ``{name: callable}`` from a mapping or a sequence of objects
    with a ``name`` attribute."""
    if metrics is None:
        return {}
    if isinstance(metrics, Mapping):
        items = list(metrics.items())
    elif isinstance(metrics, (list, tuple)):
        items = []
        for m in metrics:
            name = getattr(m, 'name', None)
            if not isinstance(name, str):
                raise TypeError('metrics given as a sequence need a str '
                                f'`name` attribute; {m!r} has none. Pass a '
                                "mapping instead, e.g. metrics={'P_end': fn}.")
            items.append((name, m))
    else:
        raise TypeError('metrics must be a mapping {name: callable} or a '
                        f'sequence of named metrics; got '
                        f'{type(metrics).__name__}.')
    out = {}
    for name, fn in items:
        if not isinstance(name, str) or not name:
            raise TypeError(f'Metric names must be non-empty strings; got '
                            f'{name!r}.')
        if not callable(fn):
            raise TypeError(f'Metric {name!r} is not callable.')
        if name in out:
            raise ValueError(f'Duplicate metric name {name!r}.')
        out[name] = fn
    return out


def _metric_selections(metrics):
    """Selections the metrics declare via a ``selections`` attribute (the
    analysis Metric protocol), de-duplicated in order."""
    sels = []
    for fn in metrics.values():
        own = getattr(fn, 'selections', None)
        if own is None:
            continue
        if isinstance(own, str):
            own = (own,)
        for s in own:
            if not isinstance(s, str):
                raise TypeError(f'Metric selections must be strings; got '
                                f'{s!r}.')
            if s not in sels:
                sels.append(s)
    return sels


def _check_workers(workers):
    if isinstance(workers, bool) or not isinstance(workers, numbers.Integral):
        raise TypeError(f'workers must be an integer; got {workers!r}.')
    if workers < 1:
        raise ValueError(f'workers must be at least 1; got {workers}.')
    return int(workers)


def _describe(error):
    return f'{type(error).__name__}: {error}'


def _unrecorded_hint(error, result):
    """For a metric's ``KeyError`` on a column the run did not record (e.g.
    ``res.final['n_spk']``, which bypasses :class:`SimulationResult`'s own
    message), how to record it; else ``''``."""
    if not (isinstance(error, KeyError) and len(error.args) == 1
            and isinstance(error.args[0], str)):
        return ''
    key = error.args[0]
    if key in result.df.columns or 'not recorded' in key:
        return ''
    return (f' (column {key!r} was not recorded; recorded: '
            f'{result.col_names}. Add it to Condition.selections or give the '
            'metric a `selections` attribute)')


def _scalar(name, value):
    """Return a metric value as a float, or raise TypeError.

    Real numbers are accepted, and booleans (Python or numpy, e.g.
    ``res.final['[S]'] < 5``) count as ``1.0`` / ``0.0``. The error for
    anything else names the type and shape, never the full repr, so a metric
    that returns a whole column does not put a huge array into every row's
    ``error`` cell.
    """
    arr = np.asarray(value)
    if arr.ndim == 0 and arr.dtype == bool:
        return float(bool(arr))
    if arr.ndim != 0 or not (np.issubdtype(arr.dtype, np.integer)
                             or np.issubdtype(arr.dtype, np.floating)):
        if arr.ndim == 0:
            got = f'{type(value).__name__} {reprlib.repr(value)}'
        else:
            got = (f'{type(value).__name__} with shape {arr.shape} and dtype '
                   f'{arr.dtype}')
        raise TypeError(f'metric {name!r} must return a real scalar; got '
                        f"{got} (use e.g. res.final['x'] or res['x'][-1]).")
    return float(arr)


def _execute(model, condition, *, metrics, extra_selections, keep_results,
             on_error):
    """Run one condition and evaluate the metrics.

    Returns
    -------
    tuple
        ``(status, error, {metric: float}, SimulationResult or None)``.
        With ``on_error='record'`` every failure (of the run or of a metric)
        is recorded here; with ``'raise'`` it propagates.
    """
    try:
        result = model.run(condition, extra_selections=extra_selections)
    except Exception as e:  # recorded as a failed row (on_error='record')
        if on_error == 'raise':
            raise
        return ('failed', _describe(e), {name: math.nan for name in metrics},
                None)
    values, errors = {}, []
    for name, fn in metrics.items():
        try:
            with np.errstate(divide='ignore', invalid='ignore'):
                values[name] = _scalar(name, fn(result))
        except Exception as e:  # recorded; the other metrics still count
            if on_error == 'raise':
                raise
            values[name] = math.nan
            errors.append(f'metric {name!r}: {_describe(e)}'
                          + _unrecorded_hint(e, result))
    if errors:
        return ('failed', '; '.join(errors), values,
                result if keep_results else None)
    return ('ok', None, values, result if keep_results else None)


def _mp_context():
    """'fork' where available (and not on macOS/Windows), else 'spawn'."""
    if ('fork' in multiprocessing.get_all_start_methods()
            and sys.platform not in ('darwin', 'win32')):
        return multiprocessing.get_context('fork')
    return multiprocessing.get_context('spawn')


def _restore_settings(model, snapshot):
    """Write back the settings in ``snapshot`` that differ on ``model``.

    ``snapshot`` is :meth:`KineticModel._user_parameter_values` of the
    private model as built (the caller's settings). Writing it back before
    each item keeps a setting that an earlier item changed -- e.g. a custom
    ``reset_func`` that writes a spike cap only for some ``reset_kwargs`` --
    from leaking into later items. Only differing values are written.
    """
    r = model._te
    for pid, value in snapshot.items():
        if r[pid] != value:
            r[pid] = value


# Per-process state of a pool worker: one model built from the blueprint, its
# settings as built (see _restore_settings) and the function to apply (or the
# error that building the model raised).
_WORKER = {}


def _worker_init(blueprint, fn):
    """Pool initializer: build this worker's private model once.

    A failure is kept and re-raised by every task this worker receives (see
    :func:`_worker_call`): an exception escaping an initializer only breaks
    the pool, and the caller would see an opaque ``BrokenProcessPool``
    instead of the real error.
    """
    _WORKER.clear()
    try:
        model = blueprint.build()
        _WORKER.update(model=model, fn=fn,
                       settings=model._user_parameter_values())
    except Exception as e:  # re-raised by each task (see above)
        _WORKER.clear()
        _WORKER['error'] = e


def _worker_call(items):
    """Apply the worker's function to a chunk of items on its model, each
    item starting from the model's settings as built."""
    error = _WORKER.get('error')
    if error is not None:
        raise error
    model, fn, settings = _WORKER['model'], _WORKER['fn'], _WORKER['settings']
    out = []
    for item in items:
        _restore_settings(model, settings)
        out.append(fn(model, item))
    return out


def _map_with_model(model, fn, items, *, workers, chunksize=None,
                    spawn_check=None):
    """Return ``[fn(private_model, item) for item in items]``.

    The generic engine behind :func:`run_many` and :func:`scan`, for any
    per-item work on a model (e.g. a calibration multistart, one fit per
    item). The caller's ``model`` is never touched: ``workers=1`` runs on one
    ``model.copy()`` in this process (no pool is created); more workers run
    on a ``ProcessPoolExecutor`` (start method from :func:`_mp_context`)
    whose initializer builds one model per worker from the model's picklable
    blueprint. Before each item the private model's settings (global
    parameters and boundary species; :meth:`KineticModel._user_parameter_values`)
    are written back to their values as built -- the caller's -- so a
    setting an earlier item changed cannot leak into a later one.

    Parameters
    ----------
    model : KineticModel
        The model to copy.
    fn : callable
        ``fn(private_model, item)``; return values are collected in order.
        Under the ``spawn`` start method it is pickled, so it must be a
        module-level function (or a ``functools.partial`` of one) with
        picklable arguments; under ``fork`` it is inherited.
    items : sequence
        One task per item (pickled to the workers when ``workers > 1``).
    workers : int
        Number of processes (validated by the caller); ``1`` runs serially.
    chunksize : int, optional
        Items per task sent to a worker. Defaults to
        ``ceil(len(items) / (4 * workers))``: a few chunks per worker balance
        the load while keeping inter-process traffic low for many short
        items. Pass ``1`` for long, uneven items.
    spawn_check : callable, optional
        Called without arguments before the pool starts when the start
        method is not ``fork``, to raise a specific ``ValueError`` for an
        unpicklable part of ``fn`` (e.g. naming a metric) ahead of the
        generic check.

    Returns
    -------
    list
        ``fn``'s return values, in the order of ``items``.

    Raises
    ------
    ValueError
        With ``workers > 1``, before any process starts: if the model's
        blueprint (its ``reset_func`` or an extra attribute) cannot be
        pickled, or, without ``fork``, if ``fn`` cannot.
    RuntimeError
        If a worker process dies (``BrokenProcessPool``), with a hint to run
        with ``workers=1`` to see the underlying error.
    Exception
        Whatever ``fn`` -- or building a worker's model -- raises propagates
        (the real exception, also from a worker); remaining tasks are
        cancelled.
    """
    items = list(items)
    if not items:
        return []
    if workers == 1:
        private = model.copy()
        # Taken on the copy as built: the caller's settings, at baseline.
        settings = private._user_parameter_values()
        out = []
        for item in items:
            _restore_settings(private, settings)
            out.append(fn(private, item))
        return out
    blueprint = model._blueprint()
    blueprint.check_picklable()
    ctx = _mp_context()
    if ctx.get_start_method() != 'fork':
        # Without fork, the initializer arguments are pickled too.
        if spawn_check is not None:
            spawn_check()
        try:
            pickle.dumps(fn)
        except (pickle.PicklingError, TypeError, AttributeError) as e:
            raise ValueError(
                f'{fn!r} is not picklable ({e}); worker processes started '
                'with spawn need module-level functions.') from e
    n = min(workers, len(items))
    size = chunksize or max(1, math.ceil(len(items) / (4 * n)))
    chunks = [items[i:i + size] for i in range(0, len(items), size)]
    with ProcessPoolExecutor(max_workers=n, mp_context=ctx,
                             initializer=_worker_init,
                             initargs=(blueprint, fn)) as pool:
        futures = [pool.submit(_worker_call, chunk) for chunk in chunks]
        try:
            return [out for f in futures for out in f.result()]
        except BrokenProcessPool as e:
            pool.shutdown(wait=True, cancel_futures=True)
            raise RuntimeError(
                'A worker process failed to start or crashed. Run with '
                'workers=1 to see the underlying error; under the spawn '
                'start method (macOS, Windows), reset_func and metrics must '
                'be importable module-level functions and scripts must call '
                "run_many/scan under `if __name__ == '__main__':`.") from e
        except BaseException:
            pool.shutdown(wait=True, cancel_futures=True)
            raise


def _validate_batch(model, conditions, metrics, workers, keep_results,
                    on_error, reserved=()):
    """Validate everything a batch needs before running anything.

    Returns
    -------
    (list of Condition, dict, list of str, int)
        Conditions, normalized metrics, extra selections, workers.
    """
    if not isinstance(model, KineticModel):
        raise TypeError(f'model must be a KineticModel; got '
                        f'{type(model).__name__}.')
    if isinstance(conditions, (Condition, str, bytes, Mapping)):
        raise TypeError('conditions must be a sequence of Condition objects.')
    conditions = list(conditions)
    for i, c in enumerate(conditions):
        if not isinstance(c, Condition):
            raise TypeError(f'conditions[{i}] is a {type(c).__name__}, not a '
                            'Condition.')
    no_reset = [i for i, c in enumerate(conditions) if not c.reset]
    if no_reset:
        raise ValueError(
            f'conditions {no_reset} have reset=False. A batch runs every '
            'condition from the baseline so results do not depend on run '
            'order or worker assignment; set the starting state through '
            'Condition.initial instead.')
    resumed = [i for i, c in enumerate(conditions)
               if c.reset_kwargs.get('mode') == 'roadrunner']
    if resumed:
        raise ValueError(
            f"conditions {resumed} use reset_kwargs mode='roadrunner', which "
            "resumes from the previous run's event and compartment state, so "
            'batch results would depend on run order or worker assignment; '
            'use the default baseline reset or a reset_func.')
    if on_error not in _ON_ERROR:
        raise ValueError(f'on_error must be one of {_ON_ERROR}; got '
                         f'{on_error!r}.')
    if not isinstance(keep_results, (bool, np.bool_)):
        raise TypeError(f'keep_results must be a bool; got {keep_results!r}.')
    workers = _check_workers(workers)
    metrics = _normalize_metrics(metrics)
    clash = [n for n in metrics if n in _RESERVED_COLUMNS or n in reserved]
    if clash:
        raise ValueError(f'Metric name(s) {clash} clash with the batch table '
                         f'columns {list(_RESERVED_COLUMNS) + list(reserved)}.')
    extra = _metric_selections(metrics)
    # Validate the union of ids, each distinct selection list and each
    # distinct reset mode once (batches can hold thousands of conditions).
    ids = dict.fromkeys(i for c in conditions
                        for i in (*c.parameters, *c.initial))
    model._validate_ids(ids)
    seen_sels = set()
    for c in conditions:
        if c.selections not in seen_sels:
            seen_sels.add(c.selections)
            model._validate_ids(model._run_selections(c, extra),
                                settable=False)
    modes, preserve = set(), {}
    for i, c in enumerate(conditions):
        mode = c.reset_kwargs.get('mode')
        if not (mode is None or isinstance(mode, str)) or mode not in modes:
            # raises for an invalid mode
            model._check_reset_kwargs({'mode': mode})
            modes.add(mode)
        kept = c.reset_kwargs.get('preserve', ())
        for sel in ((kept,) if isinstance(kept, str) else kept):
            preserve.setdefault(sel, []).append(i)
    model._check_reset_kwargs({'preserve': tuple(preserve)})
    if preserve:
        state = model._state_ids()
        carried = {sel: idx for sel, idx in preserve.items() if sel in state}
        if carried:
            where = sorted({i for idx in carried.values() for i in idx})
            raise ValueError(
                f'conditions {where} preserve state variable(s) '
                f'{list(carried)} across the reset, which carries the '
                "previous run's end state into the next run, so batch results "
                'would depend on run order or worker assignment; preserve '
                'settings only (e.g. a spike cap) and set starting values '
                'through Condition.initial.')
    return conditions, metrics, extra, workers


def _check_metrics_picklable(metrics):
    """Raise ``ValueError`` naming the first metric that cannot be pickled
    (needed by worker processes started without ``fork``)."""
    for name, fn in metrics.items():
        try:
            pickle.dumps(fn)
        except (pickle.PicklingError, TypeError, AttributeError) as e:
            raise ValueError(
                f'Metric {name!r} is not picklable ({e}); worker '
                'processes started with spawn need module-level '
                'functions or metric objects.') from e


def _run_records(model, conditions, metrics, extra, workers, keep_results,
                 on_error):
    """Execute a validated batch serially or on a process pool."""
    fn = functools.partial(_execute, metrics=metrics, extra_selections=extra,
                           keep_results=keep_results, on_error=on_error)
    return _map_with_model(
        model, fn, conditions, workers=workers,
        spawn_check=functools.partial(_check_metrics_picklable, metrics))


def _table(conditions, records, metrics, inputs=None):
    """Assemble the metrics DataFrame of a batch."""
    data = {'label': pd.Series([c.label for c in conditions], dtype=object)}
    for name, values in (inputs or {}).items():
        data[name] = pd.Series(values, dtype=float)
    data['status'] = pd.Series([rec[0] for rec in records], dtype=object)
    data['error'] = pd.Series([rec[1] for rec in records], dtype=object)
    for name in metrics:
        data[name] = pd.Series([rec[2][name] for rec in records], dtype=float)
    return pd.DataFrame(data)


# --- public API ---------------------------------------------------------------
def run_many(model, conditions, *, metrics=None, workers=1, keep_results=True,
             on_error='record'):
    """Run many conditions on a private copy of a model.

    Parameters
    ----------
    model : KineticModel
        The model to run. It is never modified: the batch runs on
        ``model.copy()`` (``workers=1``) or on one model per worker process
        built from the same picklable blueprint.
    conditions : sequence of Condition
        The runs, all with ``reset=True``, without
        ``reset_kwargs={'mode': 'roadrunner'}`` and with a ``preserve`` (if
        any) that names settings only, not state (each run starts from the
        model's baseline, so results do not depend on order or worker).
    metrics : mapping or sequence, optional
        ``{name: callable}`` or a sequence of callables with a ``name``
        attribute (e.g. :mod:`nskinetics.analysis` metrics). Each is called
        with the run's :class:`SimulationResult` and must return a real
        scalar. A metric with a ``selections`` attribute has those columns
        recorded in addition to the condition's. Evaluated with numpy
        ``divide``/``invalid`` warnings silenced, so a zero denominator gives
        ``inf``/``NaN`` instead of raising.
    workers : int
        ``1`` (default) runs serially in this process and never creates a
        pool. More starts a ``ProcessPoolExecutor`` of up to that many
        processes (``fork`` start method where available, else ``spawn``),
        each building its model once.
    keep_results : bool
        Keep every :class:`SimulationResult` in :attr:`RunBatch.results`
        (default). ``False`` keeps only the metrics, saving memory and
        inter-process traffic.
    on_error : {'record', 'raise'}
        ``'record'`` (default): a failed integration or metric becomes a row
        with ``status='failed'``, the error message and ``NaN`` for the
        metric(s) that could not be computed. ``'raise'``: the first failure
        is raised (a :class:`~nskinetics.exceptions.KineticSimulationError`
        for a failed integration) and the remaining runs are cancelled.

    Returns
    -------
    RunBatch
        Metrics table, results and conditions, in the order given.

    Raises
    ------
    TypeError
        If ``model`` is not a :class:`KineticModel`, an element of
        ``conditions`` is not a :class:`Condition`, a metric is not callable,
        or ``workers`` / ``keep_results`` have the wrong type.
    ValueError
        Before anything runs: if a condition has ``reset=False``,
        ``reset_kwargs`` ``mode='roadrunner'`` or a ``preserve`` naming a
        state variable (a species, given as an amount or a concentration, a
        compartment, or a rate-rule or event-assigned variable), uses an id
        the model does not have (or cannot set, see
        :meth:`KineticModel.run`), has invalid ``reset_kwargs``; if ``on_error``
        or ``workers`` is invalid; if metric names clash with the table's
        columns; or, with ``workers > 1``, if the model (its ``reset_func``
        or an extra attribute) cannot be pickled, or a metric cannot be
        pickled under the ``spawn`` start method.
    KineticSimulationError
        With ``on_error='raise'``, if a run fails.
    RuntimeError
        With ``workers > 1``, if a worker process dies; run with
        ``workers=1`` to see the underlying error. An error raised while a
        worker builds its model is re-raised as such.

    Notes
    -----
    Results do not depend on run order or on which worker ran a condition,
    because every run starts with a reset that fully restores the starting
    state: the default baseline reset, or a custom ``reset_func`` that
    restores every state variable (compartments, counters, flags). Such a
    function must set compartment volumes and event-assigned variables
    *before* calling ``model._te.reset()`` as well as after it (or start
    with ``model.reset(mode='baseline')``): RoadRunner's reset converts
    initial concentrations to amounts at the current volumes, so resetting
    at the previous run's end volume makes the starting concentrations
    depend on that run. For the same reason a condition's ``preserve`` may
    name settings (a spike cap) but not state variables. Every run also
    starts from the caller's settings (global parameters and boundary
    species): they are written back before each condition, so a
    ``reset_func`` that changes a setting -- such as the shipped model's
    ``reset_spike_cap`` -- cannot leak into later runs.

    On macOS and Windows, ``workers > 1`` starts processes with ``spawn``:
    each worker re-imports the main module, so call :func:`run_many` /
    :func:`scan` under ``if __name__ == '__main__':`` in scripts, and define
    metrics and any ``reset_func`` at module level (not as lambdas or
    closures). Such workers also import nskinetics afresh, which takes
    several seconds each; prefer ``workers=1`` for small batches. RoadRunner
    instances are never shared between threads or processes.

    Examples
    --------
    >>> conds = [cond.with_parameters(k=k) for k in (0.2, 0.4, 0.8)]  # doctest: +SKIP
    >>> batch = nsk.run_many(km, conds,                                # doctest: +SKIP
    ...                      metrics={'P_end': lambda res: res.final['[P]']})
    >>> batch.metrics[['label', 'status', 'P_end']]                    # doctest: +SKIP
    """
    conditions, metrics, extra, workers = _validate_batch(
        model, conditions, metrics, workers, keep_results, on_error)
    records = _run_records(model, conditions, metrics, extra, workers,
                           keep_results, on_error)
    return RunBatch(metrics=_table(conditions, records, metrics),
                    results=[rec[3] for rec in records],
                    conditions=conditions)


def _grid_values(name, values):
    if isinstance(values, (str, bytes)) or isinstance(values, Mapping):
        raise TypeError(f'grid[{name!r}] must be a sequence of numbers.')
    try:
        arr = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as e:
        raise TypeError(f'grid[{name!r}] must be a sequence of numbers: '
                        f'{e}') from e
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError(f'grid[{name!r}] must be a non-empty 1-D sequence.')
    if not np.all(np.isfinite(arr)):
        raise ValueError(f'grid[{name!r}] must be finite.')
    return tuple(float(v) for v in arr)


def _placement_for(model, condition, key):
    """Where a perturbed ``key`` belongs in a condition derived from
    ``condition``: ``'initial'`` or ``'parameters'``.

    A key the base condition already sets stays in that mapping, so the new
    value replaces the base one instead of colliding with it; any other key
    is placed by :func:`_auto_placement`. Analyses that perturb a base
    condition (:func:`scan`, finite differences, calibration targets) use
    this with :meth:`Condition.with_initial` /
    :meth:`Condition.with_parameters`.
    """
    if key in condition.initial:
        return 'initial'
    if key in condition.parameters:
        return 'parameters'
    return _auto_placement(model, key)


def _auto_placement(model, key):
    """'initial' for state ids (species, compartments, rate-rule and
    event-assigned variables, any '[x]'), else 'parameters'."""
    r = model._te
    if key.startswith('[') and key.endswith(']'):
        return 'initial'
    species = set(r.getFloatingSpeciesIds()) | set(r.getBoundarySpeciesIds())
    if key in species or key in r.getCompartmentIds():
        return 'initial'
    if key in model._state_ids():
        return 'initial'
    return 'parameters'


def _format_value(value):
    return f'{value:g}'


def scan(model, condition, grid, *, metrics, mode='product', workers=1,
         keep_results=True, on_error='record', placement=None):
    """Run a condition over a grid of parameter and/or initial values.

    Parameters
    ----------
    model : KineticModel
        The model to run (never modified; see :func:`run_many`).
    condition : Condition
        The base condition; must have ``reset=True``. Its own
        ``parameters`` and ``initial`` values are kept, with the grid values
        merged over them (:meth:`Condition.with_parameters` /
        :meth:`Condition.with_initial`).
    grid : mapping of str to sequence of float
        ``{id: values}``, e.g. ``{'k': [0.2, 0.4], '[S]': [30, 50]}``, in the
        model's own units. A key the base condition already sets stays in
        that mapping (its value is replaced). Any other key is placed in
        ``Condition.initial`` if it is state -- a bracketed selection
        (``'[S]'``), a species or compartment id, or a rate-rule or
        event-assigned variable -- and in ``Condition.parameters`` otherwise
        (global parameters such as ``'k'``). ``placement`` overrides this
        for keys the base condition does not set.
    metrics : mapping or sequence
        As for :func:`run_many` (required).
    mode : {'product', 'zip'}
        ``'product'`` (default): every combination, the first key varying
        slowest (``itertools.product`` order). ``'zip'``: the i-th values of
        every key together; all sequences must have the same length.
    workers : int
        As for :func:`run_many` (see its Notes for ``workers > 1`` on macOS
        and Windows). Defaults to ``1``.
    keep_results : bool
        As for :func:`run_many`. Defaults to ``True``.
    on_error : {'record', 'raise'}
        As for :func:`run_many`. Defaults to ``'record'``.
    placement : mapping of str to {'initial', 'parameters'}, optional
        Explicit placement of some grid keys, e.g. ``{'n_max': 'initial'}``.
        It may not contradict the base condition (a key the base condition
        sets in ``initial`` cannot be placed in ``parameters``, and vice
        versa).

    Returns
    -------
    ScanResult
        A :class:`RunBatch` whose metrics table also has one column per grid
        variable; :meth:`ScanResult.pivot` and :meth:`ScanResult.contourf`
        arrange a metric over a 1- or 2-variable grid. Generated conditions
        are labelled ``'k=0.2, [S]=30'`` (prefixed by the base label, if
        any).

    Raises
    ------
    TypeError
        For arguments of the wrong type (see :func:`run_many`), or a grid
        that is not a mapping of numeric sequences.
    ValueError
        Before anything runs: for an empty grid, an unknown ``mode``,
        ``'zip'`` sequences of different lengths, repeated values of one key
        in ``'product'`` mode, a ``placement`` naming a non-grid key or an
        unknown place or contradicting the base condition, a grid key the
        model does not have or named like a table column (``label``,
        ``status``, ``error``), metric names that clash with grid keys or
        table columns, or any reason listed for :func:`run_many`.
    KineticSimulationError
        With ``on_error='raise'``, if a run fails.
    RuntimeError
        As for :func:`run_many`.

    Notes
    -----
    As for :func:`run_many`: results do not depend on run order or worker
    (given a reset that restores the volumes before ``RoadRunner.reset()``
    and a ``preserve`` of settings only; see its Notes), and
    ``workers > 1`` under ``spawn`` (macOS, Windows) needs an
    ``if __name__ == '__main__':`` guard and module-level metrics.
    """
    if not isinstance(condition, Condition):
        raise TypeError(f'condition must be a Condition; got '
                        f'{type(condition).__name__}.')
    if not isinstance(grid, Mapping):
        raise TypeError(f'grid must be a mapping {{id: sequence of values}}; '
                        f'got {type(grid).__name__}.')
    if not grid:
        raise ValueError('grid must not be empty.')
    if mode not in _SCAN_MODES:
        raise ValueError(f'mode must be one of {_SCAN_MODES}; got {mode!r}.')
    if not isinstance(model, KineticModel):
        raise TypeError(f'model must be a KineticModel; got '
                        f'{type(model).__name__}.')
    values = {}
    for key, vals in grid.items():
        if not isinstance(key, str) or not key:
            raise TypeError(f'grid keys must be non-empty strings; got '
                            f'{key!r}.')
        if key in _RESERVED_COLUMNS:
            raise ValueError(f'grid key {key!r} clashes with the batch table '
                             f'columns {list(_RESERVED_COLUMNS)}.')
        values[key] = _grid_values(key, vals)
    if mode == 'zip':
        lengths = {k: len(v) for k, v in values.items()}
        if len(set(lengths.values())) > 1:
            raise ValueError(f"mode='zip' needs equally long sequences; got "
                             f'{lengths}.')
        points = list(zip(*values.values()))
    else:
        dupes = [k for k, v in values.items() if len(set(v)) != len(v)]
        if dupes:
            raise ValueError(f'grid values repeat for {dupes}; a product grid '
                             'needs distinct values per key.')
        points = list(itertools.product(*values.values()))
    placement = dict(placement or {})
    bad = {k: p for k, p in placement.items()
           if k not in values or p not in _PLACEMENTS}
    if bad:
        raise ValueError(f'placement {bad} must map grid keys to one of '
                         f'{_PLACEMENTS}.')
    model._validate_ids(values)
    where = {}
    for k in values:
        given = placement.get(k)
        base = ('initial' if k in condition.initial else
                'parameters' if k in condition.parameters else None)
        if given is not None and base is not None and given != base:
            raise ValueError(
                f'placement={{{k!r}: {given!r}}} contradicts the base '
                f'condition, which sets grid key {k!r} in its {base}; drop '
                f'the placement or remove {k!r} from the base condition.')
        where[k] = given or _placement_for(model, condition, k)
    keys = list(values)
    conditions = []
    for point in points:
        params = {k: v for k, v in zip(keys, point)
                  if where[k] == 'parameters'}
        initial = {k: v for k, v in zip(keys, point) if where[k] == 'initial'}
        tag = ', '.join(f'{k}={_format_value(v)}' for k, v in zip(keys, point))
        label = f'{condition.label}: {tag}' if condition.label else tag
        conditions.append(condition.with_parameters(params)
                          .with_initial(initial).replace(label=label))
    conditions, metrics, extra, workers = _validate_batch(
        model, conditions, metrics, workers, keep_results, on_error,
        reserved=tuple(keys))
    records = _run_records(model, conditions, metrics, extra, workers,
                           keep_results, on_error)
    inputs = {k: [p[i] for p in points] for i, k in enumerate(keys)}
    return ScanResult(metrics=_table(conditions, records, metrics, inputs),
                      results=[rec[3] for rec in records],
                      conditions=conditions, grid=values, mode=mode,
                      placement=where)
