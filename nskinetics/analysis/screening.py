# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Global screening of kinetics-level metrics: Morris and Sobol'.

:func:`morris` (elementary effects on one-at-a-time trajectories) ranks
parameters cheaply; :func:`sobol` (variance-based first-order and total
indices) quantifies their shares of the metric's variance, interactions
included. Both sample the parameter ranges, run every sample through
:func:`nskinetics.run_many` (serially or on worker processes, each sample a
:class:`~nskinetics.Condition` derived from the base condition with
:meth:`~nskinetics.Condition.with_parameters` /
:meth:`~nskinetics.Condition.with_initial`) and evaluate metrics on the
runs. Failed runs are excluded and counted, and so are non-finite metric
values of successful runs (e.g. a yield with a zero denominator).

The estimators are pure-numpy private cores (:func:`_morris_design`,
:func:`_morris_indices`, :func:`_sobol_design`, :func:`_sobol_indices`)
that work on any function values.

References
----------
Morris, M. D. (1991). Factorial sampling plans for preliminary
computational experiments. Technometrics 33(2), 161-174.
Campolongo, F., Cariboni, J., Saltelli, A. (2007). An effective screening
design for sensitivity analysis of large models. Environmental Modelling &
Software 22, 1509-1518.
Saltelli, A., Annoni, P., Azzini, I., Campolongo, F., Ratto, M.,
Tarantola, S. (2010). Variance based sensitivity analysis of model output.
Design and estimator for the total sensitivity index. Computer Physics
Communications 181, 259-270.
Jansen, M. J. W. (1999). Analysis of variance designs for model output.
Computer Physics Communications 117(1-2), 35-43.
Sobol', I. M. (2001). Global sensitivity indices for nonlinear mathematical
models and their Monte Carlo estimates. Mathematics and Computers in
Simulation 55(1-3), 271-280.
"""
import json
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from ..engine.conditions import Condition
from ..engine.kinetic_model import KineticModel
from ..engine import execution as _ex
from ..engine._plotting import _set_font, _style_axes
from . import _fd
from ._sampling import (ParamRange, _as_ranges, _child_seed, _rng, _scale,
                        _seed_record, _seed_sequence, _unit_sample)
from ._common import (_CSV_COLUMNS, _EMPTY, _check_ids, _count,
                      _incomplete, _read_long_csv, _scatter)
from .metrics import _check_metrics_condition

__all__ = ('MorrisResult', 'SobolResult', 'morris', 'sobol')

#: ``scipy.stats.norm.ppf(0.975)``: half-width factor of the 95% intervals.
_Z95 = 1.959963984540054
#: Share of failed runs above which a ``UserWarning`` is issued.
FAIL_WARN = 0.05


# --- statistics --------------------------------------------------------------
def _check_n_boot(n_boot):
    """``n_boot`` as an int that is 0 (no bootstrap) or at least 2: one
    resample has no spread, so its half-widths would all be ``NaN``."""
    n_boot = _count('n_boot', n_boot, 0)
    if n_boot == 1:
        raise ValueError(f'n_boot must be 0 (no bootstrap) or at least 2; '
                         f'got {n_boot}.')
    return n_boot


def _bootstrap_std(draws):
    """Sample standard deviation (``ddof=1``) over the finite draws (axis
    0); ``NaN`` where fewer than two are finite."""
    draws = np.asarray(draws, dtype=float)
    ok = np.isfinite(draws)
    n = ok.sum(axis=0)
    with np.errstate(all='ignore'):
        mean = np.where(ok, draws, 0.0).sum(axis=0) / n
        sq = np.where(ok, (draws - mean) ** 2, 0.0).sum(axis=0)
        std = np.sqrt(sq / (n - 1))
    return np.where(n >= 2, std, np.nan)


# --- Morris (pure numpy) -----------------------------------------------------
def _check_levels(levels):
    levels = _count('levels', levels, 2)
    if levels % 2:
        raise ValueError(f'levels must be even (the Morris step '
                         f'p / (2 (p - 1)) then lies on the grid); got '
                         f'{levels}.')
    return levels


def _morris_design(d, r, levels=4, seed=None):
    """Morris one-at-a-time trajectories in the unit cube.

    Each trajectory starts at a random point of the ``levels``-grid
    ``{0, 1/(p-1), ..., 1}`` and moves every coordinate once, in random
    order and direction, by ``delta = p / (2 (p - 1))`` (Morris 1991).

    Parameters
    ----------
    d, r : int
        Dimensions and number of trajectories.
    levels : int
        Even number of grid levels ``p``.
    seed : None, int or SeedSequence

    Returns
    -------
    numpy.ndarray
        ``(r * (d + 1), d)``: trajectory ``t`` is rows
        ``t * (d + 1) ... t * (d + 1) + d``.
    """
    d, r = _count('d', d, 1), _count('r', r, 1)
    levels = _check_levels(levels)
    rng = _rng(seed)
    delta = levels / (2 * (levels - 1))
    out = np.empty((r, d + 1, d))
    for t in range(r):
        base = rng.integers(0, levels // 2, size=d) / (levels - 1)
        up = rng.integers(0, 2, size=d).astype(bool)
        order = rng.permutation(d)
        x = np.where(up, base, base + delta)
        out[t, 0] = x
        for s, c in enumerate(order, 1):
            x = x.copy()
            x[c] = base[c] + delta if up[c] else base[c]
            out[t, s] = x
    return out.reshape(r * (d + 1), d)


def _morris_indices(unit, y, *, r, n_boot=200, seed=None):
    """Morris statistics from trajectories and function values.

    Elementary effects are ``EE = (f(x') - f(x)) / (x'_i - x_i)`` with
    ``x`` in unit coordinates, i.e. the change of ``f`` per full parameter
    range (per unit of ``ln p`` range for a log range).

    Parameters
    ----------
    unit : array_like
        ``(r * (d + 1), d)`` trajectories (e.g. from :func:`_morris_design`);
        each step must move exactly one coordinate, each coordinate once.
    y : array_like
        ``(r * (d + 1),)`` or ``(r * (d + 1), m)`` function values;
        ``NaN`` (a failed run) drops the effects that use it.
    r : int
        Number of trajectories.
    n_boot : int
        Bootstrap resamples (of trajectories) for ``mu_star_conf``: 0 (gives
        ``NaN``) or at least 2.
    seed : None, int or SeedSequence
        For the bootstrap.

    Returns
    -------
    dict
        ``ee`` ``(r, d, m)``; and ``(d, m)`` arrays ``mu`` (mean effect),
        ``mu_star`` (mean absolute effect, Campolongo et al. 2007),
        ``sigma`` (sample standard deviation, ``ddof=1``),
        ``mu_star_conf`` (95% bootstrap half-width of ``mu_star``) and
        ``n`` (number of finite effects).

    Raises
    ------
    ValueError
        For inconsistent shapes or steps that are not one-at-a-time.
    """
    unit = np.asarray(unit, dtype=float)
    y = np.asarray(y, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    r = _count('r', r, 1)
    n_boot = _check_n_boot(n_boot)
    if unit.ndim != 2 or unit.shape[0] != r * (unit.shape[1] + 1) \
            or y.ndim != 2 or y.shape[0] != unit.shape[0]:
        raise ValueError(f'unit must be (r * (d + 1), d) and y have one row '
                         f'per point; got {unit.shape} and {y.shape} for '
                         f'r={r}.')
    d, m = unit.shape[1], y.shape[1]
    u = unit.reshape(r, d + 1, d)
    f = y.reshape(r, d + 1, m)
    du = np.diff(u, axis=1)
    moved = du != 0
    if not (np.all(moved.sum(axis=2) == 1) and np.all(moved.sum(axis=1) == 1)):
        raise ValueError('Each trajectory step must move exactly one '
                         'coordinate, and each coordinate exactly once.')
    coord = np.argmax(moved, axis=2)
    rows = np.arange(r)[:, None]
    step = du[rows, np.arange(d)[None, :], coord]
    ee = np.empty((r, d, m))
    with np.errstate(all='ignore'):
        ee[rows, coord, :] = np.diff(f, axis=1) / step[:, :, None]
    valid = np.isfinite(ee)
    n = valid.sum(axis=0)
    with np.errstate(all='ignore'):
        mu = np.where(valid, ee, 0.0).sum(axis=0) / n
        absolute = np.where(valid, np.abs(ee), 0.0)
        mu_star = absolute.sum(axis=0) / n
        sq = np.where(valid, (ee - mu) ** 2, 0.0).sum(axis=0)
        sigma = np.where(n >= 2, np.sqrt(sq / (n - 1)), np.nan)
    mu = np.where(n > 0, mu, np.nan)
    mu_star = np.where(n > 0, mu_star, np.nan)
    conf = np.full((d, m), np.nan)
    if n_boot:
        rng = _rng(seed)
        draws = np.empty((n_boot, d, m))
        for b in range(n_boot):
            idx = rng.integers(0, r, size=r)
            with np.errstate(all='ignore'):
                draws[b] = absolute[idx].sum(axis=0) / valid[idx].sum(axis=0)
        conf = _Z95 * _bootstrap_std(draws)
    return {'ee': ee, 'mu': mu, 'mu_star': mu_star, 'sigma': sigma,
            'mu_star_conf': conf, 'n': n}


# --- Sobol' (pure numpy) -----------------------------------------------------
def _check_power_of_two(n):
    n = _count('n', n, 2)
    if n & (n - 1):
        raise ValueError(f'n must be a power of 2 (balanced Sobol\' points); '
                         f'got {n}. Try {1 << (n.bit_length() - 1)} or '
                         f'{1 << n.bit_length()}.')
    return n


def _sobol_design(d, n, seed=None):
    """Saltelli (2010) ``A, B, AB_i`` design in the unit cube.

    ``A`` and ``B`` are the two halves of a scrambled Sobol' sample of
    dimension ``2d`` (``scipy.stats.qmc.Sobol``); ``AB_i`` is ``A`` with
    column ``i`` taken from ``B``.

    Parameters
    ----------
    d : int
        Dimensions.
    n : int
        Base sample size, a power of 2.
    seed : None, int or SeedSequence

    Returns
    -------
    numpy.ndarray
        ``(n * (d + 2), d)``: ``A`` (rows ``0 .. n-1``), ``B``
        (``n .. 2n-1``), then ``AB_1 ... AB_d`` (``n`` rows each); see
        :func:`_sobol_split`.
    """
    d = _count('d', d, 1)
    n = _check_power_of_two(n)
    unit = _unit_sample(n, 2 * d, method='sobol', seed=seed)
    a, b = unit[:, :d], unit[:, d:]
    ab = np.repeat(a[None], d, axis=0)
    for i in range(d):
        ab[i, :, i] = b[:, i]
    return np.concatenate([a, b, ab.reshape(d * n, d)])


def _sobol_split(y, n, d):
    """``(fA, fB, fAB)`` from values stacked like :func:`_sobol_design`:
    shapes ``(n, m)``, ``(n, m)`` and ``(d, n, m)``."""
    y = np.asarray(y, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    if y.shape[0] != n * (d + 2):
        raise ValueError(f'Expected {n * (d + 2)} values; got {y.shape[0]}.')
    return y[:n], y[n:2 * n], y[2 * n:].reshape(d, n, y.shape[1])


def _jansen(fa, fb, fab):
    """Jansen first-order and total indices (1-D ``fa``, ``fb``; ``fab``
    ``(d, n)``); ``NaN`` for zero variance or no samples."""
    if fa.size == 0:
        nan = np.full(fab.shape[0], np.nan)
        return nan, nan.copy()
    with np.errstate(all='ignore'):
        var = np.var(np.concatenate([fa, fb]))
        s1 = (var - 0.5 * np.mean((fb - fab) ** 2, axis=1)) / var
        st = 0.5 * np.mean((fa - fab) ** 2, axis=1) / var
    if not var > 0:
        s1 = np.full_like(s1, np.nan)
        st = np.full_like(st, np.nan)
    return s1, st


def _sobol_indices(fa, fb, fab, *, n_boot=200, seed=None):
    """First-order and total Sobol' indices with Jansen estimators.

    With ``V`` the variance of ``f`` over ``A`` and ``B`` together,
    ``S1_i = (V - mean((f(B) - f(AB_i))^2) / 2) / V`` and
    ``ST_i = mean((f(A) - f(AB_i))^2) / (2 V)`` (Jansen 1999; Saltelli et
    al. 2010, Table 2).

    Parameters
    ----------
    fa, fb : array_like
        ``(n,)`` or ``(n, m)``: values on ``A`` and ``B``.
    fab : array_like
        ``(d, n)`` or ``(d, n, m)``: values on ``AB_i``.
    n_boot : int
        Bootstrap resamples (of rows) for the 95% half-widths: 0 (gives
        ``NaN``) or at least 2.
    seed : None, int or SeedSequence
        For the bootstrap.

    Returns
    -------
    dict
        ``(d, m)`` arrays ``S1``, ``ST``, ``S1_conf``, ``ST_conf`` and the
        ``(m,)`` array ``n`` of rows used. A row with a non-finite value in
        ``A``, ``B`` or any ``AB_i`` is dropped (per metric).

    Raises
    ------
    ValueError
        For inconsistent shapes.
    """
    fa, fb, fab = (np.asarray(v, dtype=float) for v in (fa, fb, fab))
    if fa.ndim == 1:
        fa, fb, fab = fa[:, None], fb[:, None], fab[..., None]
    n_boot = _check_n_boot(n_boot)
    if (fa.ndim != 2 or fb.shape != fa.shape or fab.ndim != 3
            or fab.shape[1:] != fa.shape):
        raise ValueError(f'Inconsistent shapes {fa.shape}, {fb.shape}, '
                         f'{fab.shape}: expected (n, m), (n, m), (d, n, m).')
    d, m = fab.shape[0], fa.shape[1]
    out = {key: np.full((d, m), np.nan)
           for key in ('S1', 'ST', 'S1_conf', 'ST_conf')}
    used = np.zeros(m, dtype=int)
    rng = _rng(seed)
    for c in range(m):
        a, b, ab = fa[:, c], fb[:, c], fab[:, :, c]
        valid = np.isfinite(a) & np.isfinite(b) & np.all(np.isfinite(ab),
                                                         axis=0)
        a, b, ab = a[valid], b[valid], ab[:, valid]
        used[c] = a.size
        out['S1'][:, c], out['ST'][:, c] = _jansen(a, b, ab)
        if n_boot and a.size:
            s1 = np.empty((n_boot, d))
            st = np.empty((n_boot, d))
            for k in range(n_boot):
                idx = rng.integers(0, a.size, size=a.size)
                s1[k], st[k] = _jansen(a[idx], b[idx], ab[:, idx])
            out['S1_conf'][:, c] = _Z95 * _bootstrap_std(s1)
            out['ST_conf'][:, c] = _Z95 * _bootstrap_std(st)
    out['n'] = used
    return out


# --- results -----------------------------------------------------------------
def _range_to_json(r):
    return [r.lo, r.hi, r.log]


@dataclass(eq=False, repr=False, kw_only=True)
class _ScreeningResult:
    """Shared fields, views and CSV IO of :class:`MorrisResult` and
    :class:`SobolResult` (see there)."""
    parameters: list
    metric_names: list
    ranges: dict
    unit_design: np.ndarray
    metric_values: np.ndarray
    status: list
    errors: list
    seed: int | None = None
    n_boot: int = 200
    units: dict = field(default_factory=dict)
    condition: Condition | None = None

    #: The index tables (``DataFrame`` fields), in CSV order.
    _TABLES = ()

    def __post_init__(self):
        self.parameters = list(self.parameters)
        self.metric_names = list(self.metric_names)
        self.ranges = {k: self.ranges[k] for k in self.parameters}
        self.unit_design = np.asarray(self.unit_design, dtype=float)
        self.metric_values = np.asarray(self.metric_values,
                                        dtype=float).reshape(
            len(self.unit_design), len(self.metric_names))
        self.status = list(self.status)
        self.errors = list(self.errors)
        self.units = dict(self.units)

    def __repr__(self):
        return (f'{type(self).__name__}(parameters={self.parameters}, '
                f'metrics={self.metric_names}, n_runs={self.n_runs}, '
                f'n_failed={self.n_failed})')

    # --- views ---------------------------------------------------------------
    @property
    def n_runs(self):
        """Number of runs (int)."""
        return len(self.status)

    @property
    def n_failed(self):
        """Number of failed runs, excluded from the indices (int)."""
        return sum(s != 'ok' for s in self.status)

    @property
    def design(self):
        """The sampled parameter values (``pandas.DataFrame``, one row per
        run, one column per parameter, model units)."""
        return pd.DataFrame(_scale(self.unit_design, self.ranges.values()),
                            columns=list(self.parameters))

    @property
    def samples(self):
        """The design with each run's ``status``, ``error`` and metric
        values (``pandas.DataFrame``; metrics as recorded, failed rows
        included)."""
        table = self.design
        table['status'] = pd.Series(self.status, dtype=object)
        table['error'] = pd.Series(self.errors, dtype=object)
        for k, name in enumerate(self.metric_names):
            table[name] = self.metric_values[:, k]
        return table

    def _valid_metric_values(self):
        """Metric values with failed runs set to ``NaN`` (excluded)."""
        ok = np.array([s == 'ok' for s in self.status], dtype=bool)
        return np.where(ok[:, None], self.metric_values, np.nan)

    def _metric(self, metric):
        if metric is None:
            return self.metric_names[0]
        if metric not in self.metric_names:
            raise KeyError(f'{metric!r} is not a metric; metrics: '
                           f'{self.metric_names}.')
        return metric

    # --- IO ------------------------------------------------------------------
    def _meta(self):
        return {'nskinetics': type(self).__name__,
                'parameters': self.parameters,
                'metric_names': self.metric_names, 'seed': self.seed,
                'n_boot': self.n_boot, 'units': self.units,
                'condition': (None if self.condition is None
                              else self.condition.to_dict()),
                'n_runs': self.n_runs}

    def to_csv(self, path):
        """Write the result to a long-format CSV.

        The file has the columns ``quantity, key, parameter, index, value``:
        ``meta`` rows (``key``: field), ``range`` (per ``parameter``:
        ``[lo, hi, log]``), the design in unit coordinates
        (:attr:`unit_design`, quantity ``unit_design``: ``parameter``,
        ``index`` = run), the metric values of the runs
        (:attr:`metric_values`, quantity ``metric_value``: ``key`` =
        metric, ``index``; ``status`` and ``error``: ``index``) and the
        index tables (``key`` = table, e.g. ``mu_star`` or ``ST``, then
        ``index`` = metric, ``parameter``). Values are JSON, so
        :meth:`from_csv` restores the result exactly.

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
            If the condition's ``reset_kwargs`` or the units are not
            JSON-safe.
        """
        dumps = json.dumps
        rows = [('meta', key, '', '', dumps(value, allow_nan=False))
                for key, value in self._meta().items()]
        for name, r in self.ranges.items():
            rows.append(('range', '', name, '', dumps(_range_to_json(r))))
        for j, name in enumerate(self.parameters):
            rows += [('unit_design', '', name, i, dumps(float(v)))
                     for i, v in enumerate(self.unit_design[:, j])]
        for k, name in enumerate(self.metric_names):
            rows += [('metric_value', name, '', i, dumps(float(v)))
                     for i, v in enumerate(self.metric_values[:, k])]
        rows += [('status', '', '', i, dumps(s))
                 for i, s in enumerate(self.status)]
        rows += [('error', '', '', i, dumps(e))
                 for i, e in enumerate(self.errors)]
        for table in self._TABLES:
            frame = getattr(self, table)
            for metric in self.metric_names:
                for par in self.parameters:
                    rows.append(('table', table, par, metric,
                                 dumps(float(frame.at[par, metric]))))
        pd.DataFrame(rows, columns=list(_CSV_COLUMNS)).to_csv(
            path, index=False, lineterminator='\n', encoding='utf-8')

    @classmethod
    def _read(cls, path, required=()):
        """Fields shared by both results, the meta dict and the tables;
        ``required`` lists the subclass's own ``meta`` keys."""
        kind = cls.__name__
        meta, groups = _read_long_csv(path, kind, required=(
            'parameters', 'metric_names', 'seed', 'n_boot', 'units',
            'condition', 'n_runs', *required))
        params, metrics = meta['parameters'], meta['metric_names']
        n_runs, n_p, n_m = meta['n_runs'], len(params), len(metrics)
        par_idx = {p: j for j, p in enumerate(params)}
        met_idx = {m: k for k, m in enumerate(metrics)}

        def scatter(source, quantity, shape, *axes, **kwargs):
            return _scatter(path, kind, source, quantity, shape, axes,
                            **kwargs)

        cells = scatter(groups, 'range', (n_p,), ('parameter', par_idx),
                        parse=lambda v: ParamRange(v[0], v[1], log=v[2]),
                        dtype=object)
        ranges = dict(zip(params, cells))
        unit = scatter(groups, 'unit_design', (n_runs, n_p),
                       ('index', None), ('parameter', par_idx))
        metric_values = scatter(groups, 'metric_value', (n_runs, n_m),
                                ('index', None), ('key', met_idx))
        status, errors = (
            list(scatter(groups, quantity, (n_runs,), ('index', None),
                         parse=lambda v: v, dtype=object))
            for quantity in ('status', 'error'))
        g = groups.get('table', _EMPTY)
        by_table = {key: rows for key, rows in g.groupby('key', sort=False)}
        unknown = sorted(set(by_table) - set(cls._TABLES))
        if unknown:
            raise _incomplete(path, kind, f'unknown tables {unknown}')
        tables = {name: scatter(by_table, name, (n_p, n_m),
                                ('parameter', par_idx), ('index', met_idx))
                  for name in cls._TABLES}
        condition = (None if meta['condition'] is None
                     else Condition.from_dict(meta['condition']))
        fields = dict(parameters=params, metric_names=metrics, ranges=ranges,
                      unit_design=unit, metric_values=metric_values,
                      status=status,
                      errors=errors, seed=meta['seed'],
                      n_boot=meta['n_boot'], units=meta['units'],
                      condition=condition)
        frames = {name: _frame(values, params, metrics)
                  for name, values in tables.items()}
        return meta, fields, frames


def _frame(values, parameters, metrics):
    """An index table: index ``parameter``, one column per metric."""
    return pd.DataFrame(np.asarray(values, dtype=float).reshape(
        len(parameters), len(metrics)),
        index=pd.Index(list(parameters), name='parameter', dtype=object),
        columns=pd.Index(list(metrics), name='metric', dtype=object))


@dataclass(eq=False, repr=False, kw_only=True)
class MorrisResult(_ScreeningResult):
    """Outcome of :func:`morris`.

    Index tables have one row per parameter (index ``parameter``) and one
    column per metric, e.g. ``g.mu_star['P_titer']``. Effects are changes
    of the metric per full parameter range (for a log range, per full range
    of ``ln p``), in the metric's units.

    Parameters
    ----------
    parameters : list of str
        Screened ids, in design-column order.
    metric_names : list of str
    ranges : dict
        ``{id: ParamRange}``.
    unit_design : numpy.ndarray
        ``(r * (d + 1), d)`` design in unit-cube coordinates (see
        :attr:`design` for model units).
    metric_values : numpy.ndarray
        ``(n_runs, n_metrics)`` metric values of every run as recorded
        (``NaN`` where a run or metric failed).
    status, errors : list
        Per run: ``'ok'``/``'failed'`` and the error message or ``None``.
    mu, mu_star, sigma, mu_star_conf : pandas.DataFrame
        Mean elementary effect, mean absolute effect (the ranking measure),
        standard deviation of the effects (nonlinearity and interactions)
        and the 95% bootstrap half-width of ``mu_star``.
    r : int
        Number of trajectories.
    levels : int
        Grid levels ``p``.
    seed : int or None
        Reproduces the design and the bootstrap when passed back as
        ``seed=``.
    n_boot : int
    units : dict
        The model's units (not to be confused with :attr:`unit_design`).
    condition : Condition or None
        The base condition.

    Raises
    ------
    ValueError
        If ``metric_values`` cannot be shaped ``(n_runs, n_metrics)`` for
        the design.
    KeyError
        If ``ranges`` lacks a parameter.
    """
    mu: pd.DataFrame
    mu_star: pd.DataFrame
    sigma: pd.DataFrame
    mu_star_conf: pd.DataFrame
    r: int
    levels: int

    _TABLES = ('mu', 'mu_star', 'sigma', 'mu_star_conf')

    @property
    def elementary_effects(self):
        """Elementary effects (``numpy.ndarray``, ``(r, d, n_metrics)``;
        ``NaN`` where a run failed)."""
        return _morris_indices(self.unit_design, self._valid_metric_values(),
                               r=self.r, n_boot=0)['ee']

    @property
    def n_used(self):
        """Number of finite elementary effects behind each index
        (``pandas.DataFrame``, parameter x metric; at most :attr:`r`)."""
        n = _morris_indices(self.unit_design, self._valid_metric_values(),
                            r=self.r, n_boot=0)['n']
        return _frame(n, self.parameters, self.metric_names).astype('int64')

    def rank(self, metric=None):
        """Parameters ranked by ``mu_star`` for one metric.

        Parameters
        ----------
        metric : str, optional
            Defaults to the first metric.

        Returns
        -------
        pandas.DataFrame
            Columns ``parameter``, ``mu``, ``mu_star``, ``sigma``,
            ``mu_star_conf``, sorted by decreasing ``mu_star`` (``NaN``
            last).

        Raises
        ------
        KeyError
            For an unknown metric.
        """
        metric = self._metric(metric)
        table = pd.DataFrame({
            'parameter': pd.Series(self.parameters, dtype=object),
            **{name: getattr(self, name)[metric].to_numpy(dtype=float)
               for name in self._TABLES}})
        key = table['mu_star'].to_numpy()
        order = np.argsort(np.where(np.isnan(key), np.inf, -key),
                           kind='stable')
        return table.iloc[order].reset_index(drop=True)

    def plot(self, metric=None, ax=None):
        """``sigma`` against ``mu_star`` for one metric (Morris plot).

        Parameters
        ----------
        metric : str, optional
            Defaults to the first metric.
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure otherwise.

        Returns
        -------
        (matplotlib.figure.Figure, matplotlib.axes.Axes)

        Raises
        ------
        KeyError
            For an unknown metric.

        Notes
        -----
        Points far right matter most; points high above the diagonal act
        nonlinearly or through interactions. Horizontal bars are the 95%
        bootstrap intervals of ``mu_star``.
        """
        metric = self._metric(metric)
        _set_font()
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        x = self.mu_star[metric].to_numpy(dtype=float)
        y = self.sigma[metric].to_numpy(dtype=float)
        conf = self.mu_star_conf[metric].to_numpy(dtype=float)
        ax.errorbar(x, y, xerr=np.where(np.isfinite(conf), conf, 0.0),
                    fmt='o', color='steelblue', ecolor='gray', capsize=2)
        for name, xi, yi in zip(self.parameters, x, y):
            if np.isfinite(xi) and np.isfinite(yi):
                ax.annotate(name, (xi, yi), textcoords='offset points',
                            xytext=(4, 4))
        ax.set_xlabel('mu*')
        ax.set_ylabel('sigma')
        ax.set_title(metric)
        _style_axes(ax)
        return fig, ax

    def _meta(self):
        return {**super()._meta(), 'r': self.r, 'levels': self.levels}

    @classmethod
    def from_csv(cls, path):
        """Read a result written by :meth:`to_csv`.

        Parameters
        ----------
        path : str or path-like
            CSV written by :meth:`MorrisResult.to_csv`.

        Returns
        -------
        MorrisResult
            Equal to the one written (exact values, ``NaN`` in place).

        Raises
        ------
        ValueError
            If the file was not written by :meth:`MorrisResult.to_csv`, or is
            incomplete (a missing ``meta`` key, or an entry of the design,
            the metric values, the run status or the index tables missing
            or duplicated, e.g. a truncated or edited file).
        """
        meta, fields, frames = cls._read(path, required=('r', 'levels'))
        return cls(**fields, **frames, r=meta['r'], levels=meta['levels'])


@dataclass(eq=False, repr=False, kw_only=True)
class SobolResult(_ScreeningResult):
    """Outcome of :func:`sobol`.

    Index tables have one row per parameter (index ``parameter``) and one
    column per metric, e.g. ``s.ST['P_titer']``.

    Parameters
    ----------
    parameters, metric_names, ranges, status, errors
        As for :class:`MorrisResult`.
    seed, n_boot, units, condition
        As for :class:`MorrisResult`.
    unit_design : numpy.ndarray
        ``(n * (d + 2), d)`` design in unit-cube coordinates: ``A``, ``B``,
        then ``AB_1 ... AB_d`` (``n`` rows each).
    metric_values : numpy.ndarray
        ``(n_runs, n_metrics)`` metric values of every run as recorded.
    S1, ST : pandas.DataFrame
        First-order and total Sobol' indices (Jansen estimators); ``NaN``
        for a metric that does not vary over the design (zero variance).
    S1_conf, ST_conf : pandas.DataFrame
        95% bootstrap half-widths: the intervals are ``S1 +- S1_conf`` and
        ``ST +- ST_conf``.
    n : int
        Base sample size (a power of 2).

    Raises
    ------
    ValueError, KeyError
        As for :class:`MorrisResult`.
    """
    S1: pd.DataFrame
    ST: pd.DataFrame
    S1_conf: pd.DataFrame
    ST_conf: pd.DataFrame
    n: int

    _TABLES = ('S1', 'ST', 'S1_conf', 'ST_conf')

    @property
    def n_used(self):
        """Number of base rows behind each metric's indices
        (``pandas.Series`` indexed by metric; at most :attr:`n`): rows whose
        ``A``, ``B`` and every ``AB_i`` run gave a finite value."""
        fa, fb, fab = _sobol_split(self._valid_metric_values(), self.n,
                                   len(self.parameters))
        valid = (np.isfinite(fa) & np.isfinite(fb)
                 & np.all(np.isfinite(fab), axis=0))
        return pd.Series(valid.sum(axis=0).astype('int64'),
                         index=pd.Index(self.metric_names, name='metric',
                                        dtype=object), name='n_used')

    def rank(self, metric=None):
        """Parameters ranked by ``ST`` for one metric.

        Parameters
        ----------
        metric : str, optional
            Defaults to the first metric.

        Returns
        -------
        pandas.DataFrame
            Columns ``parameter``, ``S1``, ``ST``, ``S1_conf``, ``ST_conf``,
            sorted by decreasing ``ST`` (``NaN`` last).

        Raises
        ------
        KeyError
            For an unknown metric.
        """
        metric = self._metric(metric)
        table = pd.DataFrame({
            'parameter': pd.Series(self.parameters, dtype=object),
            **{name: getattr(self, name)[metric].to_numpy(dtype=float)
               for name in self._TABLES}})
        key = table['ST'].to_numpy()
        order = np.argsort(np.where(np.isnan(key), np.inf, -key),
                           kind='stable')
        return table.iloc[order].reset_index(drop=True)

    def plot(self, metric=None, ax=None):
        """Bar chart of ``S1`` and ``ST`` with 95% intervals for one metric.

        Parameters
        ----------
        metric : str, optional
            Defaults to the first metric.
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure otherwise.

        Returns
        -------
        (matplotlib.figure.Figure, matplotlib.axes.Axes)

        Raises
        ------
        KeyError
            For an unknown metric.
        """
        metric = self._metric(metric)
        _set_font()
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        x = np.arange(len(self.parameters))
        width = 0.38
        for offset, name, color in ((-width / 2, 'S1', 'steelblue'),
                                    (width / 2, 'ST', 'darkorange')):
            values = getattr(self, name)[metric].to_numpy(dtype=float)
            conf = getattr(self, f'{name}_conf')[metric].to_numpy(dtype=float)
            ax.bar(x + offset, np.where(np.isfinite(values), values, 0.0),
                   width, yerr=np.where(np.isfinite(conf), conf, 0.0),
                   label=name, color=color, capsize=2)
        ax.set_xticks(x)
        ax.set_xticklabels(self.parameters, rotation=90)
        ax.set_ylabel('Sobol\' index')
        ax.set_title(metric)
        ax.axhline(0.0, color='gray', linewidth=0.6)
        ax.legend(edgecolor='white')
        _style_axes(ax)
        return fig, ax

    def _meta(self):
        return {**super()._meta(), 'n': self.n}

    @classmethod
    def from_csv(cls, path):
        """Read a result written by :meth:`to_csv`.

        Parameters
        ----------
        path : str or path-like
            CSV written by :meth:`SobolResult.to_csv`.

        Returns
        -------
        SobolResult
            Equal to the one written (exact values, ``NaN`` in place).

        Raises
        ------
        ValueError
            If the file was not written by :meth:`SobolResult.to_csv`, or is
            incomplete (a missing ``meta`` key, or an entry of the design,
            the metric values, the run status or the index tables missing
            or duplicated, e.g. a truncated or edited file).
        """
        meta, fields, frames = cls._read(path, required=('n',))
        return cls(**fields, **frames, n=meta['n'])


# --- public API --------------------------------------------------------------
def _prepare(model, condition, bounds, metrics, workers):
    """Validate the shared arguments (before any run); return the ranges,
    normalized metrics, placements and workers."""
    if not isinstance(model, KineticModel):
        raise TypeError(f'model must be a KineticModel; got '
                        f'{type(model).__name__}.')
    if not isinstance(condition, Condition):
        raise TypeError(f'condition must be a Condition; got '
                        f'{type(condition).__name__}.')
    ranges = _as_ranges(bounds)
    if metrics is None:
        raise ValueError('metrics are required: {name: metric}.')
    _check_ids(model, 'bounds', ranges)
    _check_ids(model, 'metrics',
               _ex._metric_selections(_ex._normalize_metrics(metrics)),
               settable=False)
    _, metrics, _, workers = _ex._validate_batch(
        model, [condition], metrics, workers, False, 'record',
        reserved=tuple(ranges))
    if not metrics:
        raise ValueError('metrics must name at least one metric.')
    # e.g. At(t) outside the condition's time span would fail every run
    _check_metrics_condition(metrics, condition)
    placement = {name: _ex._placement_for(model, condition, name)
                 for name in ranges}
    # Every sample condition must be constructible (e.g. '[S]' cannot join
    # an 'S' amount in the base condition); the values are irrelevant.
    try:
        _fd.pin(model, condition, {n: r.lo for n, r in ranges.items()},
                placement)
    except ValueError as e:
        raise ValueError(
            f'bounds: {e} List each species once, in the same form (amount '
            "id or '[concentration]') as the condition's own initial "
            'values.') from e
    # user -> morris/sobol -> _prepare -> the helper -> warnings.warn
    _fd.warn_initial_assignment_inputs(model, ranges, stacklevel=4)
    return ranges, metrics, placement, workers


def _evaluate(model, condition, ranges, metrics, placement, unit, workers):
    """Run every design point; ``(metric_values, status, errors)``."""
    names = list(ranges)
    natural = _scale(unit, ranges.values())
    conditions = [_fd.pin(model, condition,
                          {n: float(v) for n, v in zip(names, row)}, placement)
                  for row in natural]
    batch = _ex.run_many(model, conditions, metrics=metrics, workers=workers,
                         keep_results=False)
    table = batch.metrics
    values = table[list(metrics)].to_numpy(dtype=float)
    status = table['status'].tolist()
    errors = table['error'].tolist()
    n_failed = sum(s != 'ok' for s in status)
    if n_failed > FAIL_WARN * len(status):
        first = next(e for s, e in zip(status, errors) if s != 'ok')
        warnings.warn(f'{n_failed} of {len(status)} runs failed '
                      f'({100 * n_failed / len(status):.1f}%) and are '
                      f'excluded from the indices; first error: {first}',
                      UserWarning, stacklevel=3)
    return values, status, errors


def _warn_nonfinite(values, ok, metric_names):
    """Warn (``UserWarning``) about metrics that are not finite on more than
    ``FAIL_WARN`` of the successful runs (their rows/effects are excluded,
    silently otherwise: such runs are not counted as failed)."""
    n_ok = int(ok.sum())
    bad = np.sum(~np.isfinite(values[ok]), axis=0)
    hit = [f'{name!r} on {int(b)} of {n_ok}'
           for name, b in zip(metric_names, bad) if b > FAIL_WARN * n_ok]
    if hit:
        warnings.warn(f"Metric values are not finite (NaN/inf) on successful "
                      f"runs: {', '.join(hit)}; those runs are excluded from "
                      'the indices of that metric (see n_used).',
                      UserWarning, stacklevel=3)


def morris(model, condition, bounds, metrics, *, r=20, levels=4, seed=None,
           workers=1, n_boot=200):
    """Morris elementary-effects screening of metrics over parameter ranges.

    ``r`` one-at-a-time trajectories through the ranges (``r * (d + 1)``
    runs for ``d`` parameters) give, per parameter and metric, the mean
    effect ``mu``, the mean absolute effect ``mu_star`` (the ranking
    measure; Campolongo et al. 2007) and the spread ``sigma``
    (nonlinearity or interactions).

    Parameters
    ----------
    model : KineticModel
        Never modified (runs use private copies; see
        :func:`nskinetics.run_many`).
    condition : Condition
        The base run (``reset=True``); each sample merges its values into it
        with :meth:`~nskinetics.Condition.with_parameters` /
        :meth:`~nskinetics.Condition.with_initial` (placement as in
        :func:`nskinetics.scan`).
    bounds : mapping
        ``{id: (lo, hi)}`` (uniform) or ``{id: ParamRange(lo, hi, log=True)}``
        (uniform in ``ln``), model units: global parameters, initial values
        (``'[S]'``), compartments...
    metrics : mapping or sequence
        ``{name: metric}`` or named metrics (see :class:`Metric`); at least
        one.
    r : int
        Number of trajectories (``>= 2``; default 20).
    levels : int
        Even number of grid levels ``p`` (default 4); the step is
        ``p / (2 (p - 1))`` of the range.
    seed : None, int or numpy.random.SeedSequence
        Seeds the design and the bootstrap; the same seed gives identical
        results. ``None`` draws fresh entropy, recorded in
        :attr:`MorrisResult.seed`.
    workers : int
        Processes for the runs (default ``1``: serial).
    n_boot : int
        Bootstrap resamples for ``mu_star_conf`` (default 200; 0 skips;
        otherwise at least 2).

    Returns
    -------
    MorrisResult

    Raises
    ------
    TypeError, ValueError
        For invalid arguments, unknown or unsettable ids, metric names that
        clash with parameter names or table columns, a metric that cannot
        be evaluated on the condition (e.g. :class:`At` with ``t`` outside
        its time span), or a condition a batch cannot run (see
        :func:`nskinetics.run_many`), before any run; the message names the
        argument (``bounds:``, ``metrics:``) where it is not obvious.

    Warns
    -----
    UserWarning
        If more than 5% of the runs fail; failed runs are excluded (an
        elementary effect needs both of its runs) and counted in
        :attr:`MorrisResult.n_failed`. Also if a metric is not finite on
        more than 5% of the successful runs (e.g. a zero denominator); those
        effects are excluded too (see :attr:`MorrisResult.n_used`). And if a
        bounds id is read by an SBML initial assignment (see Notes).

    Notes
    -----
    An elementary effect is the change of the metric per full range of the
    parameter (of ``ln p`` for a log range), in the metric's units, so
    ``mu_star`` compares parameters on the scale of their ranges. Effects
    are differences across ``2/3`` of a range (``levels=4``), so
    they capture nonlinearity that local derivatives miss, but ``mu_star``
    is a screening measure, not a variance share: use :func:`sobol` on the
    parameters that matter. Results hold for the condition given.

    Every sampled id is set in every run, so the other sampled ids stay at
    their sampled values: a compartment is varied at fixed species amounts
    unless a concentration ``'[S]'`` is also sampled, which then holds the
    concentration. An id read by an SBML initial assignment (e.g. ``S_init``
    in ``S = S_init``) is set after the reset has evaluated it, so the
    assigned initial values do not follow and the id looks inert; sample
    ``'init(S_init)'`` instead.

    Examples
    --------
    >>> bounds = {'k': ParamRange(0.05, 5, log=True), 'threshold': (1, 20)}
    >>> g = morris(km, cond, bounds, metrics, r=30, seed=0)  # doctest: +SKIP
    >>> g.rank('P_titer')                                    # doctest: +SKIP
    """
    r = _count('r', r, 2)
    levels = _check_levels(levels)
    n_boot = _check_n_boot(n_boot)
    ss = _seed_sequence(seed)
    ranges, metrics, placement, workers = _prepare(model, condition, bounds,
                                                   metrics, workers)
    unit = _morris_design(len(ranges), r, levels, seed=_child_seed(ss, 0))
    values, status, errors = _evaluate(model, condition, ranges, metrics,
                                       placement, unit, workers)
    ok = np.array([s == 'ok' for s in status], dtype=bool)
    _warn_nonfinite(values, ok, list(metrics))
    stats = _morris_indices(unit, np.where(ok[:, None], values, np.nan),
                            r=r, n_boot=n_boot, seed=_child_seed(ss, 1))
    names, metric_names = list(ranges), list(metrics)
    return MorrisResult(
        parameters=names, metric_names=metric_names, ranges=ranges,
        unit_design=unit, metric_values=values, status=status,
        errors=errors,
        **{key: _frame(stats[key], names, metric_names)
           for key in MorrisResult._TABLES},
        r=r, levels=levels, seed=_seed_record(ss), n_boot=n_boot,
        units=dict(model._units), condition=condition)


def sobol(model, condition, bounds, metrics, *, n=512, seed=None, workers=1,
          n_boot=200):
    """Variance-based (Sobol') indices of metrics over parameter ranges.

    Uses the Saltelli (2010) design -- matrices ``A`` and ``B`` from one
    scrambled Sobol' sample (``scipy.stats.qmc.Sobol``) and ``AB_i`` (``A``
    with column ``i`` from ``B``), ``n * (d + 2)`` runs for ``d``
    parameters -- and Jansen estimators for the first-order index ``S1``
    (the share of the metric's variance due to the parameter alone) and the
    total index ``ST`` (including all its interactions).

    Parameters
    ----------
    model : KineticModel
        Never modified (see :func:`morris`).
    condition : Condition
        The base run (see :func:`morris`).
    bounds : mapping
        ``{id: (lo, hi)}`` or ``{id: ParamRange}`` (see :func:`morris`).
        Parameters are sampled independently and uniformly (in ``ln`` for a
        log range).
    metrics : mapping or sequence
        As for :func:`morris`.
    n : int
        Base sample size, a power of 2 (default 512).
    seed : None, int or numpy.random.SeedSequence
        As for :func:`morris`.
    workers : int
        Processes for the runs (default ``1``: serial).
    n_boot : int
        Bootstrap resamples for the 95% half-widths (default 200; 0 skips;
        otherwise at least 2).

    Returns
    -------
    SobolResult

    Raises
    ------
    TypeError, ValueError
        For invalid arguments (``n`` not a power of 2 included), unknown
        ids, clashing names, a metric that cannot be evaluated on the
        condition or an unrunnable condition, before any run (see
        :func:`morris`).

    Warns
    -----
    UserWarning
        If more than 5% of the runs fail; a base row with any failed run
        (in ``A``, ``B`` or an ``AB_i``) is dropped for that metric and the
        failures are counted in :attr:`SobolResult.n_failed`. Also if a
        metric is not finite on more than 5% of the successful runs; such
        rows are dropped for that metric too (see
        :attr:`SobolResult.n_used`). And if a bounds id is read by an SBML
        initial assignment (see :func:`morris`), or if a metric does not
        vary over the design (see Notes).

    Notes
    -----
    ``S1`` near 0 with a large ``ST`` marks a parameter that acts through
    interactions; ``sum(S1) < 1`` measures how much of the variance comes
    from interactions. Estimates carry sampling error (see the ``_conf``
    tables): small negative values are noise. A metric that does not vary
    over the design (zero variance, e.g. a constant, or an
    :class:`EventCount` that always reaches its cap) has undefined indices:
    its ``S1``, ``ST`` and their ``_conf`` entries are ``NaN``, with a
    ``UserWarning``. Cost: ``n * (d + 2)`` runs
    (e.g. 6144 for ``n = 512``, ``d = 10``); screen with :func:`morris`
    first. Results hold for the condition and the ranges given. Sampled ids
    are set as for :func:`morris` (see its Notes on compartments and ids
    read by initial assignments).

    Examples
    --------
    >>> bounds = {'k': ParamRange(0.05, 5, log=True), 'threshold': (1, 20)}
    >>> s = sobol(km, cond, bounds, metrics, n=512, seed=0)  # doctest: +SKIP
    >>> s.ST['P_titer']                                      # doctest: +SKIP
    """
    n = _check_power_of_two(n)
    n_boot = _check_n_boot(n_boot)
    ss = _seed_sequence(seed)
    ranges, metrics, placement, workers = _prepare(model, condition, bounds,
                                                   metrics, workers)
    d = len(ranges)
    unit = _sobol_design(d, n, seed=_child_seed(ss, 0))
    values, status, errors = _evaluate(model, condition, ranges, metrics,
                                       placement, unit, workers)
    ok = np.array([s == 'ok' for s in status], dtype=bool)
    _warn_nonfinite(values, ok, list(metrics))
    fa, fb, fab = _sobol_split(np.where(ok[:, None], values, np.nan), n, d)
    stats = _sobol_indices(fa, fb, fab, n_boot=n_boot, seed=_child_seed(ss, 1))
    names, metric_names = list(ranges), list(metrics)
    # With rows to use, all-NaN indices mean a zero variance (see _jansen).
    constant = [m for k, m in enumerate(metric_names)
                if stats['n'][k] > 0 and np.isnan(stats['S1'][:, k]).all()]
    if constant:
        warnings.warn(f'Metric(s) {constant} do not vary over the design '
                      "(zero variance), so their Sobol' indices are "
                      'undefined (NaN).', UserWarning, stacklevel=2)
    return SobolResult(
        parameters=names, metric_names=metric_names, ranges=ranges,
        unit_design=unit, metric_values=values, status=status,
        errors=errors,
        **{key: _frame(stats[key], names, metric_names)
           for key in SobolResult._TABLES},
        n=n, seed=_seed_record(ss), n_boot=n_boot,
        units=dict(model._units), condition=condition)
