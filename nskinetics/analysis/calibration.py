# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Fed-batch-aware parameter estimation with uncertainty and identifiability.

Measured time courses (:class:`Dataset`) are compared with model
:class:`Observable` values in one or more :class:`Experiment` objects, each a
:class:`~nskinetics.Condition` that carries its own feeding parameters and
initial state. :class:`FitParameter` objects say what is estimated -- global
parameters, initial values (also one per experiment, ``scope='local'``) and
observable scales or offsets -- and a :class:`CalibrationProblem` ties them
together as a weighted least-squares problem.
:func:`fit_parameters` solves it from several starts (or a global stage and
a polish) and returns a :class:`FitResult` with estimates, Wald confidence
intervals, correlations, and identifiability and goodness-of-fit flags.

Every model evaluation is a full :meth:`~nskinetics.KineticModel.run` of
an experiment's condition on a private copy of the model, through the
event-aware finite-difference engine of :mod:`nskinetics.analysis._fd`,
which also gives the Jacobian (central differences in each parameter's
fitting space, with event-count flags). The caller's model is never
modified.

Noise models (``sigma_i`` of data point ``i``; residuals are
``(model_i - data_i) / sigma_i``):

- ``'additive'``: ``sigma_i = sigma`` (absolute, in the observable's units);
- ``'proportional'``: ``sigma_i = rel_sigma * max(|y_i|, relative_floor *
  max |y|)``, with ``y`` the *data* of that observable in that experiment
  (the floor is a detection limit, as in
  :func:`~nskinetics.analysis.identifiability`);
- ``'combined'``: ``sigma_i = sqrt(sigma**2 + proportional_i**2)``.

A ``sigma`` column in the data overrides the observable's noise model for
those points. Every ``sigma_i`` is computed once from the data and held
fixed, never from the model's predictions, so the objective is a plain
weighted least-squares problem (its Jacobian is the Gauss-Newton Fisher
information that the covariance uses); a prediction-based sigma would bias
the estimates toward larger predictions.

With ``absolute_sigma=True`` (the default) the sigmas are the standard
deviations of the data. With ``absolute_sigma=False`` they are weights
only: the noise level is estimated from the residuals, and the covariance
is scaled by the reduced chi-square. In this module *relative* always
refers to the proportional noise model (``rel_sigma``, ``relative_floor``),
never to ``absolute_sigma=False``.

Uncertainty of a fit, from cheapest to most faithful:

- Wald (:attr:`FitResult.ci`): ``x -+ 1.96 SE`` from the Gauss-Newton
  covariance ``(J^T J)^-1`` at the optimum (a pseudo-inverse if the
  Jacobian is rank deficient; see :attr:`FitResult.covariance`), in the
  fitting space (``ln p``
  for log parameters; with ``absolute_sigma=False`` the covariance is
  scaled by the reduced chi-square and ``1.96`` becomes the Student-t
  quantile ``t(0.975, dof)``);
- Laplace sampling (:meth:`FitResult.sample`, :meth:`FitResult.predict`):
  draws from the same local Gaussian (a Student-t with ``dof`` degrees of
  freedom when ``absolute_sigma=False``), for propagating the uncertainty
  to predictions;
- profile likelihood (:func:`profile_likelihood`): each parameter is fixed
  on a grid and the others re-optimized, so the interval follows the
  actual objective (non-quadratic, bounded or flat).

The first two share the Gaussian approximation: they understate the
uncertainty of a parameter whose profile is not quadratic. Check flagged
parameters with a profile.
"""
import dataclasses
import functools
import inspect
import json
import math
import pickle
import textwrap
import warnings
from collections.abc import Mapping
from dataclasses import KW_ONLY, dataclass, field
from types import MappingProxyType

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, NullLocator, ScalarFormatter
from scipy.optimize import least_squares
from scipy.stats import chi2 as chi2_distribution
from scipy.stats import norm
from scipy.stats import t as t_distribution

from ..engine import execution as _ex
from ..engine._plotting import _bf, _set_font, _style_axes
from ..engine.conditions import Condition, _hashable
from ..engine.kinetic_model import KineticModel
from ..exceptions import FittingError, KineticSimulationError
from ..utils._fit import hybrid_global_optimize
from . import _fd
from ._common import (_CSV_COLUMNS, _EMPTY, _check_ids, _count, _finite,
                      _incomplete, _read_long_csv, _scatter)
from ._sampling import (ParamRange, _child_seed, _gaussian, _seed_record,
                        _seed_sequence, _unit_sample)
from .identifiability import (NULL_LOADING, RANK_RTOL, REL_SE_MAX,
                              RELATIVE_FLOOR, _check_rank_rtol,
                              _check_relative_floor, _identifiable, _rel_se,
                              _relative_sigma, _svd_analysis)

__all__ = ('Dataset', 'Observable', 'Experiment', 'FitParameter',
           'CalibrationProblem', 'fit_parameters', 'FitResult',
           'profile_likelihood', 'ProfileResult')

#: Loss functions of ``scipy.optimize.least_squares`` accepted by
#: :class:`CalibrationProblem`.
LOSSES = ('linear', 'soft_l1', 'huber', 'cauchy')
#: Observable noise models.
NOISE_MODELS = ('additive', 'proportional', 'combined')
#: :class:`FitParameter` scopes.
SCOPES = ('global', 'local')
#: :func:`fit_parameters` methods.
FIT_METHODS = ('multistart', 'hybrid')
#: Start-point samplers (scrambled Sobol' or Latin hypercube; see
#: :func:`nskinetics.analysis._sampling._unit_sample`).
SAMPLERS = ('sobol', 'lhs')
#: Pairs of estimates whose correlation exceeds this (in absolute value) are
#: flagged.
CORRELATION_MAX = 0.95
#: At full rank, an equilibrated Fisher information ``D J^T J D`` (every
#: column scaled to unit norm, ``D_j = 1 / ||J_j||``, see
#: :attr:`FitResult.covariance`; a simulation column the outputs do not
#: measurably depend on is set to 0 first, which makes the Jacobian rank
#: deficient instead) whose condition number
#: (``FitResult.condition_number ** 2``) exceeds this is flagged as
#: ill-conditioned: its columns are nearly collinear. Beyond
#: ``1 / rank_rtol`` in ``cond(J D)`` (``1e6`` by default) the Jacobian
#: counts as rank deficient instead.
CONDITION_MAX = 1e8
#: Smallest residual (in units of sigma) given to every data point of an
#: evaluation that failed (integration failure or non-finite values). Within
#: one optimizer start the penalty is ``max(PENALTY_RESIDUAL, 10 * max |r|)``
#: over the residuals ``r`` of its start point, so a failed point always
#: looks worse than the start and a trust-region step into it is rejected.
PENALTY_RESIDUAL = 1e3
#: More free parameters than this trigger a ``UserWarning``.
N_FREE_WARN = 15
#: Two-sided 95% standard-normal quantile of the Wald intervals with
#: absolute sigmas (with ``absolute_sigma=False``, Student's t with ``dof``
#: degrees of freedom is used).
Z_95 = float(norm.ppf(0.975))
#: A fit with absolute sigmas whose ``chi2`` (with a robust loss: ``2 *
#: cost``, the robust objective, which is never above ``chi2``) exceeds this
#: quantile of chi-square with ``dof`` degrees of freedom is flagged
#: ``lack_of_fit`` (a chi-square goodness-of-fit test at the 0.1% level;
#: conservative for a robust loss).
LACK_OF_FIT_LEVEL = 0.999
#: Linear-space parameters are differenced with the absolute step
#: ``rel_step * max(|p|, LINEAR_STEP_FLOOR * (hi - lo))``: the floor guards
#: a parameter near 0 against integration noise, but a ``|p|`` far below
#: ``LINEAR_STEP_FLOOR * (hi - lo)`` gets a step that is not small relative
#: to ``p`` (see :class:`FitParameter`).
LINEAR_STEP_FLOOR = 1e-3
#: Options of :func:`~nskinetics.utils.hybrid_global_optimize` used by
#: ``fit_parameters(method='hybrid')`` unless overridden by
#: ``global_options``.
HYBRID_DEFAULTS = {'de_popsize': 10, 'de_maxiter': 30, 'do_basinhop': False,
                   'polish_top_k': 1, 'local_maxiter': 50}
# What global_options may set: the keyword options of hybrid_global_optimize
# except those fit_parameters sets itself; the counts among them must be >= 1.
_RESERVED_OPTIONS = ('objective', 'bounds', 'parallel', 'workers',
                     'random_state')
_GLOBAL_OPTIONS = tuple(
    name for name, p in inspect.signature(
        hybrid_global_optimize).parameters.items()
    if p.kind is inspect.Parameter.KEYWORD_ONLY
    and name not in _RESERVED_OPTIONS)
_COUNT_OPTIONS = ('de_popsize', 'de_maxiter', 'bh_niter', 'local_maxiter',
                  'polish_top_k')
#: Confidence level of profile-likelihood intervals.
PROFILE_LEVEL = 0.95
#: Profile threshold with absolute sigmas: the ``PROFILE_LEVEL`` quantile of
#: chi-square with one degree of freedom (3.84), on ``2 * (cost -
#: cost_min)`` -- the increase of the weighted sum of squared residuals for
#: the linear loss. With ``absolute_sigma=False`` the profile, divided by the
#: reduced chi-square, is compared with ``F(1, dof)``'s quantile (``t**2``)
#: instead.
PROFILE_THRESHOLD = float(chi2_distribution.ppf(PROFILE_LEVEL, 1))
#: A profile that goes this far (in the threshold's units) below the fit's
#: optimum triggers a ``UserWarning``: the fit had not found the minimum.
PROFILE_DROP_WARN = 0.01
#: Width at which :meth:`FitResult.summary` wraps its text lines.
SUMMARY_WIDTH = 88
#: Default number of output times of the model curves of
#: :meth:`FitResult.plot` (per experiment, data times added).
PLOT_POINTS = 200
_CSV_FORMAT = 'FitResult'
_PROFILE_FORMAT = 'ProfileResult'
_PROFILE_SCHEMA = {'side': object, 'x': float, 'value': float,
                   'cost': float, 'chi2': float, 'status': object,
                   'success': bool, 'nfev': np.int64,
                   'n_failed_evals': np.int64, 'message': object}
_DATA_COLUMNS = ('experiment', 'observable', 'time', 'value', 'sigma')
_RESIDUAL_SCHEMA = {'experiment': object, 'observable': object,
                    'time': float, 'data': float, 'model': float,
                    'sigma': float, 'weighted_residual': float}
_START_SCHEMA = {'start': np.int64, 'status': object, 'success': bool,
                 'best': bool, 'cost': float, 'chi2': float,
                 'nfev': np.int64, 'njev': np.int64,
                 'n_failed_evals': np.int64, 'message': object}
_FLAG_SCHEMA = {'flag': object, 'parameters': object, 'value': float,
                'message': object}


class _Text(str):
    """A ``str`` that displays as itself (e.g. :meth:`FitResult.summary` in
    a notebook shows the formatted text, not a quoted repr)."""

    def __repr__(self):
        return str(self)


def _wald_quantile(absolute_sigma, dof):
    """The two-sided ``PROFILE_LEVEL`` quantile of a Wald interval: the
    standard normal's for absolute sigmas; Student's t with ``dof`` degrees
    of freedom when the noise level is estimated from the residuals
    (``absolute_sigma=False``; ``NaN`` without degrees of freedom)."""
    if absolute_sigma:
        return Z_95
    if dof > 0:
        return float(t_distribution.ppf(0.5 + 0.5 * PROFILE_LEVEL, dof))
    return math.nan


def _rows(mask, limit=5):
    """The positions of the first ``limit`` true entries, for messages."""
    idx = np.flatnonzero(np.asarray(mask))
    out = idx[:limit].tolist()
    return f'{out}{" ..." if idx.size > limit else ""}'


def _typed_frame(columns, schema):
    """A DataFrame with exactly the columns and dtypes of ``schema``, built
    from ``{column: sequence}`` (the one constructor of every frame a
    :class:`FitResult` stores, so a CSV round trip gives equal frames)."""
    return pd.DataFrame({c: pd.Series(list(columns[c]), dtype=dtype)
                         for c, dtype in schema.items()})


def _py(value):
    """A JSON-safe Python scalar for a numpy or Python value."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, tuple):
        return list(value)
    return value


# --- data --------------------------------------------------------------------
def _labels(arg, series):
    """Labels as a list of non-empty strings (``str()`` of each value)."""
    missing = series.isna().to_numpy()
    if missing.any():
        raise ValueError(f'Dataset: {arg} labels are missing in rows '
                         f'{_rows(missing)}.')
    out = [str(v) for v in series]
    empty = np.array([not v for v in out], dtype=bool)
    if empty.any():
        raise ValueError(f'Dataset: {arg} labels are empty in rows '
                         f'{_rows(empty)}.')
    return out


def _numbers(arg, series):
    """A column as a float array (``ValueError`` for non-numeric entries)."""
    if series.dtype == bool:
        raise TypeError(f'Dataset: {arg} must be numeric; got booleans.')
    try:
        return pd.to_numeric(series, errors='raise').to_numpy(dtype=float)
    except (TypeError, ValueError) as e:
        raise ValueError(f'Dataset: {arg} must be numeric ({e}).') from e


class Dataset:
    """Measured time courses in long format.

    One row per measurement, with the columns ``experiment`` (label),
    ``observable`` (name of an :class:`Observable`), ``time`` (in the
    model's time units), ``value`` (in the observable's units) and
    ``sigma`` (the measurement's standard deviation in the same units;
    ``NaN`` where not given, i.e. the observable's noise model applies).

    Parameters
    ----------
    frame : pandas.DataFrame
        With columns ``experiment``, ``observable``, ``time``, ``value`` and
        optionally ``sigma``; other columns are ignored. Labels are
        converted with ``str()``.

    Raises
    ------
    TypeError
        If ``frame`` is not a DataFrame.
    ValueError
        If a column is missing, there are no rows, a label is missing or
        empty, a time or value is not finite, or a given sigma is not
        finite and positive.

    Notes
    -----
    ``dataset[label]`` is the :class:`Dataset` of one experiment
    (``KeyError`` if there is none), ``label in dataset`` tests for one,
    and ``len(dataset)`` is the number of rows. Times are checked against
    each experiment's simulated span (``t0`` to ``t_end`` of its condition)
    by :class:`Experiment`.
    """

    def __init__(self, frame):
        if not isinstance(frame, pd.DataFrame):
            raise TypeError(f'Dataset needs a pandas DataFrame; got '
                            f'{type(frame).__name__} (see Dataset.from_frame, '
                            'from_long_csv and from_wide).')
        missing = [c for c in _DATA_COLUMNS[:4] if c not in frame.columns]
        if missing:
            raise ValueError(f'Dataset: missing column(s) {missing}; got '
                             f'{list(frame.columns)}.')
        if not len(frame):
            raise ValueError('Dataset: no rows.')
        experiment = _labels('experiment', frame['experiment'])
        observable = _labels('observable', frame['observable'])
        time = _numbers('time', frame['time'])
        value = _numbers('value', frame['value'])
        sigma = (_numbers('sigma', frame['sigma']) if 'sigma' in frame.columns
                 else np.full(len(frame), np.nan))
        for arg, arr in (('time', time), ('value', value)):
            bad = ~np.isfinite(arr)
            if bad.any():
                raise ValueError(f'Dataset: {arg} is not finite in rows '
                                 f'{_rows(bad)}.')
        with np.errstate(invalid='ignore'):
            bad = ~np.isnan(sigma) & ~(np.isfinite(sigma) & (sigma > 0))
        if bad.any():
            raise ValueError(f'Dataset: sigma must be finite and > 0 (or '
                             f'missing) in rows {_rows(bad)}.')
        self._frame = pd.DataFrame({
            'experiment': pd.Series(experiment, dtype=object),
            'observable': pd.Series(observable, dtype=object),
            'time': time, 'value': value, 'sigma': sigma})

    # --- constructors -------------------------------------------------------
    @classmethod
    def from_frame(cls, df, *, experiment='experiment',
                   observable='observable', time='time', value='value',
                   sigma=None):
        """Build a dataset from a long-format DataFrame.

        Parameters
        ----------
        df : pandas.DataFrame
            One row per measurement.
        experiment, observable, time, value : str
            Names of the columns holding the experiment label, the
            observable name, the time (model time units) and the measured
            value.
        sigma : str, optional
            Name of a column of standard deviations (``NaN`` cells: the
            observable's noise model applies). None (default): the column
            named ``sigma`` if ``df`` has one, as ``Dataset(df)`` does; drop
            that column to ignore it.

        Returns
        -------
        Dataset

        Raises
        ------
        TypeError, ValueError
            As for :class:`Dataset`, or if a named column is missing.
        """
        if not isinstance(df, pd.DataFrame):
            raise TypeError(f'df must be a pandas DataFrame; got '
                            f'{type(df).__name__}.')
        if sigma is None and 'sigma' in df.columns:
            sigma = 'sigma'
        cols = {'experiment': experiment, 'observable': observable,
                'time': time, 'value': value}
        if sigma is not None:
            cols['sigma'] = sigma
        missing = [v for v in cols.values() if v not in df.columns]
        if missing:
            raise ValueError(f'Dataset: column(s) {missing} not found; got '
                             f'{list(df.columns)}.')
        return cls(pd.DataFrame({k: df[v].to_numpy() for k, v in
                                 cols.items()}))

    @classmethod
    def from_long_csv(cls, path, *, experiment='experiment', time='time',
                      observable='observable', value='value', sigma=None):
        """Read a long-format CSV file (one row per measurement).

        Parameters
        ----------
        path : str or path-like
            The CSV file.
        experiment, time, observable, value : str
            Names of the columns holding the experiment label, the time
            (model time units), the observable name and the value. The
            defaults are the canonical names (``'experiment'``, ``'time'``,
            ``'observable'``, ``'value'``), as in :meth:`from_frame`.
        sigma : str, optional
            Name of a column of standard deviations (empty cells: the
            observable's noise model applies). None (default): the column
            named ``sigma`` if the file has one (see :meth:`from_frame`).

        Returns
        -------
        Dataset

        Raises
        ------
        ValueError
            As for :meth:`from_frame`.
        """
        df = pd.read_csv(path, dtype={experiment: str, observable: str})
        return cls.from_frame(df, experiment=experiment,
                              observable=observable, time=time, value=value,
                              sigma=sigma)

    @classmethod
    def from_wide(cls, df, *, experiment, time, columns, sigma=None):
        """Build a dataset from a wide DataFrame (one column per observable).

        Parameters
        ----------
        df : pandas.DataFrame
            One row per sampling time and experiment.
        experiment, time : str
            Names of the experiment-label and time columns (time in model
            time units). A single run without a label column needs one
            added, e.g. ``Dataset.from_wide(df.assign(run='B1'),
            experiment='run', time='t_h', columns=['S', 'P'])``; any label
            works, since :class:`Experiment` relabels a one-experiment
            dataset to its own name.
        columns : sequence of str or mapping
            The measured columns: names (each is also the observable name)
            or ``{column: observable name}``.
        sigma : mapping, optional
            ``{measured column: sigma column}`` for columns with
            per-measurement standard deviations.

        Returns
        -------
        Dataset
            Rows grouped by observable (in the order of ``columns``), each
            in the row order of ``df``. Empty (``NaN``) cells are not
            measurements and are dropped.

        Raises
        ------
        TypeError, ValueError
            As for :class:`Dataset`, if a named column is missing, or if
            ``columns`` is empty.
        """
        if not isinstance(df, pd.DataFrame):
            raise TypeError(f'df must be a pandas DataFrame; got '
                            f'{type(df).__name__}.')
        if isinstance(columns, str):
            columns = [columns]
        mapping = (dict(columns) if isinstance(columns, Mapping)
                   else {c: c for c in columns})
        if not mapping:
            raise ValueError('from_wide needs at least one measured column.')
        sigma = dict(sigma or {})
        unknown = [c for c in sigma if c not in mapping]
        if unknown:
            raise ValueError(f'sigma names {unknown}, which are not among the '
                             f'measured columns {list(mapping)}.')
        needed = [experiment, time, *mapping, *sigma.values()]
        missing = [c for c in needed if c not in df.columns]
        if missing:
            raise ValueError(f'Dataset: column(s) {missing} not found; got '
                             f'{list(df.columns)}.')
        parts = []
        for col, name in mapping.items():
            values = df[col]
            keep = values.notna().to_numpy()
            part = pd.DataFrame({
                'experiment': df[experiment].to_numpy()[keep],
                'observable': [name] * int(keep.sum()),
                'time': df[time].to_numpy()[keep],
                'value': values.to_numpy()[keep],
                'sigma': (df[sigma[col]].to_numpy()[keep] if col in sigma
                          else np.full(int(keep.sum()), np.nan))})
            parts.append(part)
        return cls(pd.concat(parts, ignore_index=True))

    # --- access -------------------------------------------------------------
    @property
    def frame(self):
        """The rows (``pandas.DataFrame``, a copy) with columns
        ``experiment, observable, time, value, sigma``."""
        return self._frame.copy()

    @property
    def experiments(self):
        """Experiment labels in order of first appearance (list of str)."""
        return list(dict.fromkeys(self._frame['experiment']))

    @property
    def observables(self):
        """Observable names in order of first appearance (list of str)."""
        return list(dict.fromkeys(self._frame['observable']))

    def __len__(self):
        return len(self._frame)

    def __contains__(self, experiment):
        return str(experiment) in set(self._frame['experiment'])

    def __getitem__(self, experiment):
        """The rows of one experiment, as a :class:`Dataset` (the label is
        compared after ``str()``, as labels are stored).

        Raises
        ------
        KeyError
            If there is no such experiment.
        """
        experiment = str(experiment)
        mask = (self._frame['experiment'] == experiment).to_numpy()
        if not mask.any():
            raise KeyError(f'No experiment {experiment!r} in the dataset; '
                           f'have {self.experiments}.')
        return Dataset(self._frame[mask].reset_index(drop=True))

    def __eq__(self, other):
        if not isinstance(other, Dataset):
            return NotImplemented
        return self._frame.equals(other._frame)

    __hash__ = None

    def __repr__(self):
        return (f'Dataset({len(self)} rows, experiments={self.experiments}, '
                f'observables={self.observables})')


# --- observables -------------------------------------------------------------
def _scale_value(arg, value):
    """A scale/offset: a non-empty str (a fit-parameter name) or a finite
    float."""
    if isinstance(value, str):
        if not value:
            raise ValueError(f'{arg} must be a number or a non-empty '
                             'fit-parameter name.')
        return value
    return _finite(arg, value)


def _positive(arg, value):
    value = _finite(arg, value)
    if value <= 0:
        raise ValueError(f'{arg} must be > 0; got {value!r}.')
    return value


@dataclass(frozen=True)
class Observable:
    """A measured quantity: a model selection, transformed, scaled and
    offset, with a noise model.

    The model value is ``scale * transform(x) + offset`` with ``x`` the
    recorded column ``selection`` of a run (on the run's output grid).

    Parameters
    ----------
    name : str
        Name used in the data (:class:`Dataset` ``observable`` column).
    selection : str
        RoadRunner selection to record, e.g. ``'[S]'``.
    scale, offset : float or str
        Numbers, or the name of a :class:`FitParameter` (whose ``target`` is
        that name) to estimate them. Defaults ``1.0`` and ``0.0``.
    noise : {'additive', 'proportional', 'combined'}
        Noise model, giving the standard deviation ``sigma_i`` of data
        point ``i`` (residuals are ``(model_i - data_i) / sigma_i``):

        - ``'additive'``: ``sigma_i = sigma``;
        - ``'proportional'``: ``sigma_i = rel_sigma * max(|y_i|,
          relative_floor * max|y|)``, with ``y`` the data of this
          observable in that experiment and ``relative_floor`` from
          :class:`CalibrationProblem` (default ``1e-3``): the bounded
          relative form of ``sigma=('relative', rel)`` in
          :func:`~nskinetics.analysis.identifiability` (same
          ``relative_floor``), computed from the data rather than from
          the nominal model output;
        - ``'combined'``: ``sigma_i = sqrt(sigma**2 +
          proportional_i**2)``.

        A value in the :class:`Dataset` ``sigma`` column overrides the
        noise model for that point. All sigmas are computed once from the
        data and held fixed. See :ref:`calibration-noise`. Defaults to
        ``'additive'``.
    sigma : float, optional
        Absolute standard deviation (observable units); required by
        ``'combined'``; for ``'additive'``, data points without a ``sigma``
        column value need it.
    rel_sigma : float, optional
        Relative standard deviation; required by ``'proportional'`` and
        ``'combined'``.
    transform : callable, optional
        ``transform(x) -> array`` of the same shape, applied to the recorded
        column (no strings are evaluated). It should be vectorized. It acts
        on this one selection: a signal that combines several species
        (e.g. total sugars, or the OD of several biomass pools) needs an
        assignment rule in the model, whose id is then the selection. With
        ``workers > 1`` under the ``spawn`` start method (macOS, Windows)
        it must be picklable, i.e. a module-level function such as
        ``np.log10``, not a lambda or a nested function. It must not
        raise for model values outside its domain: it should return
        ``NaN`` or ``inf`` there instead (as numpy ufuncs such as
        ``np.log`` do), which is penalized like a failed run and counted
        in ``n_failed_evals``. An exception it raises propagates out of
        :meth:`CalibrationProblem.residuals`, :func:`fit_parameters`,
        :func:`profile_likelihood` and :meth:`FitResult.plot`.

    Raises
    ------
    TypeError, ValueError
        For an invalid field: an empty name or selection, a zero scale, an
        unknown noise model, a missing or non-positive sigma or rel_sigma
        the noise model needs, a sigma or rel_sigma it does not use, the
        same fit-parameter name for scale and offset, or a non-callable
        transform.
    """
    name: str
    selection: str
    _: KW_ONLY
    scale: float | str = 1.0
    offset: float | str = 0.0
    noise: str = 'additive'
    sigma: float | None = None
    rel_sigma: float | None = None
    transform: object = None

    def __post_init__(self):
        put = functools.partial(object.__setattr__, self)
        for arg in ('name', 'selection'):
            v = getattr(self, arg)
            if not isinstance(v, str) or not v:
                raise TypeError(f'Observable.{arg} must be a non-empty '
                                f'string; got {v!r}.')
        put('scale', _scale_value(f'scale of {self.name!r}', self.scale))
        if not isinstance(self.scale, str) and self.scale == 0:
            raise ValueError(f'Observable {self.name!r}: scale must be '
                             'nonzero (a zero scale makes the observable a '
                             'constant).')
        put('offset', _scale_value(f'offset of {self.name!r}', self.offset))
        if isinstance(self.scale, str) and self.scale == self.offset:
            raise ValueError(f'Observable {self.name!r}: scale and offset '
                             f'name the same fit parameter {self.scale!r}.')
        if self.noise not in NOISE_MODELS:
            raise ValueError(f'Observable {self.name!r}: noise must be one of '
                             f'{NOISE_MODELS}; got {self.noise!r}.')
        needs_abs = self.noise == 'combined'
        needs_rel = self.noise in ('proportional', 'combined')
        if self.sigma is not None:
            if self.noise == 'proportional':
                raise ValueError(f'Observable {self.name!r}: sigma is not '
                                 "used by noise='proportional'; give "
                                 "rel_sigma, or use noise='combined'.")
            put('sigma', _positive(f'sigma of {self.name!r}', self.sigma))
        elif needs_abs:
            raise ValueError(f'Observable {self.name!r}: noise={self.noise!r} '
                             'needs sigma.')
        if self.rel_sigma is not None:
            if not needs_rel:
                raise ValueError(f'Observable {self.name!r}: rel_sigma is not '
                                 "used by noise='additive'; use "
                                 "noise='proportional' or 'combined'.")
            put('rel_sigma', _positive(f'rel_sigma of {self.name!r}',
                                       self.rel_sigma))
        elif needs_rel:
            raise ValueError(f'Observable {self.name!r}: noise={self.noise!r} '
                             'needs rel_sigma.')
        if self.transform is not None and not callable(self.transform):
            raise TypeError(f'Observable {self.name!r}: transform must be '
                            'callable.')

    @property
    def parameters(self):
        """Names of the fit parameters this observable uses for ``scale``
        and ``offset`` (tuple of str; empty if both are numbers)."""
        return tuple(v for v in (self.scale, self.offset)
                     if isinstance(v, str))

    def _apply(self, x, scale, offset):
        """``scale * transform(x) + offset`` for a recorded column ``x``."""
        x = np.asarray(x, dtype=float)
        with np.errstate(all='ignore'):
            y = x if self.transform is None else np.asarray(
                self.transform(x), dtype=float)
            if y.shape != x.shape:
                raise ValueError(f'Observable {self.name!r}: transform '
                                 f'returned shape {y.shape} for input shape '
                                 f'{x.shape}.')
            return scale * y + offset

    def values(self, result, parameters=None):
        """Model values of this observable on a run's output grid.

        Parameters
        ----------
        result : SimulationResult
            A run that recorded ``selection``.
        parameters : mapping, optional
            ``{name: value}`` for a ``scale`` / ``offset`` given by name.

        Returns
        -------
        numpy.ndarray
            ``scale * transform(result[selection]) + offset``, one value per
            output time (observable units).

        Raises
        ------
        KeyError
            If ``selection`` was not recorded.
        ValueError
            If a named scale or offset has no value in ``parameters``, or
            ``transform`` changes the shape.
        """
        parameters = parameters or {}
        missing = [p for p in self.parameters if p not in parameters]
        if missing:
            raise ValueError(f'Observable {self.name!r} needs values for '
                             f'{missing}.')
        scale, offset = (parameters[v] if isinstance(v, str) else v
                         for v in (self.scale, self.offset))
        return self._apply(result[self.selection], scale, offset)

    def noise_sigma(self, data, *, relative_floor=RELATIVE_FLOOR):
        """Standard deviations of data points under this noise model.

        Parameters
        ----------
        data : array_like
            Measured values of this observable in one experiment (all of
            them: the proportional floor is relative to their largest
            magnitude).
        relative_floor : float
            Floor of the proportional term, as a fraction of
            ``max |data|``.

        Returns
        -------
        numpy.ndarray
            One sigma per point (observable units); ``NaN`` everywhere for
            ``'additive'`` noise without ``sigma`` (the data must then
            supply it).
        """
        y = np.asarray(data, dtype=float)
        if self.noise == 'additive':
            value = math.nan if self.sigma is None else self.sigma
            return np.full(y.shape, value)
        with np.errstate(all='ignore'):
            peak = float(np.abs(y).max(initial=0.0))
            prop = _relative_sigma(y, self.rel_sigma, peak, relative_floor)
            if self.noise == 'proportional':
                return prop
            return np.sqrt(self.sigma ** 2 + prop ** 2)


# --- experiments and parameters ----------------------------------------------
class Experiment:
    """One experiment: a run condition and the data measured in it.

    Parameters
    ----------
    name : str
        Unique name (local parameters are called ``'target@name'``).
    condition : Condition
        How the experiment is run (initial state, feeding parameters,
        ``t0``, ``t_end`` or ``times``, reset options). It must reset
        (``reset=True``). Calibration runs it with its output grid replaced
        by :attr:`times` and its selections by the observables' selections
        (``Condition.replace``); every other field is kept.
    data : Dataset or pandas.DataFrame
        The measurements: a :class:`Dataset` of one experiment (relabelled
        to ``name``; e.g. ``data['FB1']``), or of several, of which the rows
        labelled ``name`` are used; or a DataFrame with columns
        ``observable, time, value`` (and optionally ``sigma``).
    observables : sequence of str or Observable, optional
        This experiment's complete set of observables, which must define
        every observable in its data. Each entry is the name of a problem
        observable, or an :class:`Observable` used instead of the
        same-named problem observable (or as an extra one under a new
        name). Problem observables that are not listed are not used in
        this experiment. For example, to give ``S`` its own sigma in an
        experiment that also measures ``P``: ``observables=['P',
        Observable('S', '[S]', sigma=0.5)]``. A single name or
        :class:`Observable` is accepted. ``None`` (default): the problem's
        observables.

    Attributes
    ----------
    name : str
        As given.
    condition : Condition
        As given (the output grid and selections are replaced only in the
        runs of a :class:`CalibrationProblem`).
    data : Dataset
        This experiment's rows only, relabelled to ``name``.
    observables : tuple or None
        As given (a single name or :class:`Observable` as a 1-tuple);
        ``None`` means the problem's observables.
    times : numpy.ndarray
        Output times of every run of the experiment: the sorted union of
        ``t0`` and the data times (model time units).

    Raises
    ------
    TypeError
        For arguments of the wrong type.
    ValueError
        For an empty name, a condition with ``reset=False``, no data, a
        multi-experiment dataset without rows named ``name``, data times
        outside ``[t0, t_end]`` of the condition (``times[-1]`` for a
        condition with ``times``), data only at ``t0``, or duplicate
        observable names.
    """

    def __init__(self, name, condition, data, observables=None):
        if not isinstance(name, str) or not name:
            raise TypeError(f'Experiment name must be a non-empty string; got '
                            f'{name!r}.')
        if not isinstance(condition, Condition):
            raise TypeError(f'Experiment {name!r}: condition must be a '
                            f'Condition; got {type(condition).__name__}.')
        if not condition.reset:
            raise ValueError(f'Experiment {name!r}: the condition has '
                             'reset=False; every evaluation must start from '
                             'the baseline (set the starting state through '
                             'Condition.initial).')
        if isinstance(data, pd.DataFrame):
            frame = data.copy()
            if 'experiment' not in frame.columns:
                frame.insert(0, 'experiment', name)
            data = Dataset(frame)
        if not isinstance(data, Dataset):
            raise TypeError(f'Experiment {name!r}: data must be a Dataset or '
                            f'a DataFrame; got {type(data).__name__}.')
        if len(data.experiments) > 1:
            if name not in data:
                raise ValueError(f'Experiment {name!r}: the dataset holds '
                                 f'experiments {data.experiments} but none '
                                 f'named {name!r}; pass data[label].')
            data = data[name]
        frame = data.frame
        frame['experiment'] = pd.Series([name] * len(frame), dtype=object)
        data = Dataset(frame)
        t0 = condition.t0
        t_end = (condition.t_end if condition.times is None
                 else condition.times[-1])
        t = frame['time'].to_numpy(dtype=float)
        outside = (t < t0) | (t > t_end)
        if outside.any():
            bad = sorted(set(t[outside].tolist()))[:5]
            raise ValueError(
                f'Experiment {name!r}: data times {bad} lie outside the '
                f'simulated span [{t0:g}, {t_end:g}] of its condition (times '
                'are in model time units).')
        times = np.unique(np.concatenate([[t0], t]))
        if times.size < 2:
            raise ValueError(f'Experiment {name!r}: all data are at t0 = '
                             f'{t0:g}; at least one later time is needed.')
        if observables is not None:
            if isinstance(observables, (str, Observable)):
                observables = [observables]
            observables = tuple(observables)
            names = []
            for o in observables:
                if isinstance(o, Observable):
                    names.append(o.name)
                elif isinstance(o, str) and o:
                    names.append(o)
                else:
                    raise TypeError(f'Experiment {name!r}: observables must '
                                    'be names or Observable objects; got '
                                    f'{o!r}.')
            dup = sorted({n for n in names if names.count(n) > 1})
            if dup:
                raise ValueError(f'Experiment {name!r}: observables {dup} are '
                                 'listed more than once.')
        self.name = name
        self.condition = condition
        self.data = data
        self.observables = observables
        times.setflags(write=False)
        self.times = times

    def __repr__(self):
        return (f'Experiment({self.name!r}, {len(self.data)} data points, '
                f'observables={self.data.observables})')


@dataclass(frozen=True)
class FitParameter:
    """A quantity to estimate, with its search bounds.

    Parameters
    ----------
    target : str
        A global parameter id (``'k'``), an initial-value selection
        (``'[S]'``, a compartment, ``'init(k)'``), or the name an
        :class:`Observable` uses for its ``scale`` or ``offset``.
    bounds : (float, float)
        ``(lo, hi)`` in the target's own (model or observable) units;
        ``lo < hi``, both finite, and ``lo > 0`` when ``log=True``.
    log : bool
        Estimate ``ln(value)`` (default ``True``): bounds, starts, steps and
        the covariance are in log space and the confidence interval is
        asymmetric in natural units. With ``log=False`` the
        finite-difference step is ``rel_step * max(|p|, 1e-3 (hi - lo))``
        (see :class:`CalibrationProblem`), so it depends on the width of
        the bounds: when ``|p|`` is far below ``1e-3 (hi - lo)`` (more
        than about ``1e4``-``1e5`` times smaller than ``hi - lo``, e.g. a
        rate constant of 0.5 on ``(0, 1e6)``) the central difference is
        not small relative to ``p`` and can leave the bounds, and the
        Jacobian, the estimate and the SE are inaccurate without any flag.
        Use ``log=True`` for a positive quantity, or bounds within a few
        decades of the expected size.
    scope : {'global', 'local'}
        ``'global'`` (default): one value shared by every experiment.
        ``'local'``: one free value per experiment, named
        ``'target@experiment'`` (e.g. the initial substrate of each run).
    initial : float or mapping, optional
        Start value of the first optimizer start (natural units, within the
        bounds); for ``scope='local'`` also ``{experiment: value}`` (an
        experiment not listed, and a ``None`` value, use the default below;
        names of experiments not in a problem are ignored, so one parameter
        can serve several problems). ``None`` (default): the value the run
        starts from, i.e. the experiment's condition entry if it sets the
        target, else the model's current value (a global parameter: in the
        first experiment it applies to); the midpoint of the bounds in the
        fitting space (geometric for ``log=True``) for observable
        scales/offsets, for values outside the bounds, or if that run
        fails.

    Raises
    ------
    TypeError, ValueError
        For an empty target, invalid bounds, a non-positive lower bound with
        ``log=True``, an unknown scope, or an initial value outside the
        bounds (or a mapping for a global parameter).
    """
    target: str
    bounds: tuple
    _: KW_ONLY
    log: bool = True
    scope: str = 'global'
    initial: object = None

    def __post_init__(self):
        put = functools.partial(object.__setattr__, self)
        if not isinstance(self.target, str) or not self.target:
            raise TypeError(f'FitParameter.target must be a non-empty string; '
                            f'got {self.target!r}.')
        if not isinstance(self.log, (bool, np.bool_)):
            raise TypeError(f'FitParameter {self.target!r}: log must be a '
                            'bool.')
        put('log', bool(self.log))
        if self.scope not in SCOPES:
            raise ValueError(f'FitParameter {self.target!r}: scope must be '
                             f'one of {SCOPES}; got {self.scope!r}.')
        b = self.bounds
        if isinstance(b, (str, bytes, Mapping)) or not hasattr(b, '__len__') \
                or len(b) != 2:
            raise TypeError(f'FitParameter {self.target!r}: bounds must be a '
                            f'(lo, hi) pair; got {b!r}.')
        try:
            rng = ParamRange(b[0], b[1], log=self.log)
        except (TypeError, ValueError) as e:
            hint = ''
            if self.log and isinstance(e, ValueError):
                try:
                    ParamRange(b[0], b[1], log=False)
                except (TypeError, ValueError):
                    pass
                else:
                    # valid on a linear scale: only the log default fails
                    hint = (' Pass log=False for a quantity that can be 0 or '
                            'negative, such as an offset.')
            raise type(e)(f'FitParameter {self.target!r}: {e}{hint}') from e
        put('bounds', rng.bounds)
        init = self.initial
        if init is None:
            return
        if isinstance(init, Mapping):
            if self.scope != 'local':
                raise TypeError(f'FitParameter {self.target!r}: a mapping of '
                                'initial values needs scope=\'local\'.')
            values = {}
            for k, v in init.items():
                if not isinstance(k, str):
                    raise TypeError(f'FitParameter {self.target!r}: '
                                    'initial keys must be experiment names; '
                                    f'got {k!r}.')
                # None: that experiment uses the default start (see
                # initial_for)
                values[k] = None if v is None else self._inside(v)
            put('initial', MappingProxyType(values))
        else:
            put('initial', self._inside(init))

    def __hash__(self):
        # the generated hash fails on a mapping ``initial`` (a read-only
        # MappingProxyType); hash it as a frozenset, as Condition does
        return hash(tuple(_hashable(getattr(self, f.name))
                          for f in dataclasses.fields(self)))

    def __getstate__(self):
        # MappingProxyType cannot be pickled (nor deep-copied): ship a
        # plain dict and re-validate/re-freeze on the way back in
        state = {f.name: getattr(self, f.name)
                 for f in dataclasses.fields(self)}
        if isinstance(state['initial'], Mapping):
            state['initial'] = dict(state['initial'])
        return state

    def __setstate__(self, state):
        self.__init__(**state)

    def _inside(self, value):
        value = _finite(f'initial value of {self.target!r}', value)
        lo, hi = self.bounds
        if not lo <= value <= hi:
            raise ValueError(f'FitParameter {self.target!r}: initial value '
                             f'{value!r} is outside the bounds {self.bounds}.')
        return value

    @property
    def range(self):
        """The bounds as a :class:`~nskinetics.analysis.ParamRange`
        (``log`` as given)."""
        return ParamRange(*self.bounds, log=self.log)

    def initial_for(self, experiment=None):
        """The start value of the first optimizer start.

        Parameters
        ----------
        experiment : str, optional
            Experiment name (used by a ``{experiment: value}`` mapping of a
            local parameter; ignored otherwise).

        Returns
        -------
        float or None
            The value in natural units, or ``None`` for the default (the
            value the run starts from, else the midpoint of the bounds in
            the fitting space; see ``initial``), which a
            :class:`CalibrationProblem` resolves.
        """
        if isinstance(self.initial, Mapping):
            return self.initial.get(experiment)
        return self.initial


# --- the problem -------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class _Block:
    """The data of one observable in one experiment: absolute positions in
    the residual vector, positions on the experiment's output grid, values,
    sigmas, and the scale and offset (numbers, or free-parameter names)."""
    observable: Observable
    rows: np.ndarray
    grid: np.ndarray
    scale: float | str
    offset: float | str


def _linear_scale(value, lo, hi):
    """The size of a linear-space parameter, ``max(|p|, LINEAR_STEP_FLOOR *
    (hi - lo))`` (natural units; positive for finite bounds ``lo < hi``):
    the basis of its finite-difference step, independent of the
    parameter's unit (the rank test of :class:`FitResult` equilibrates
    the Jacobian's columns instead)."""
    with np.errstate(all='ignore'):
        return np.maximum(np.abs(value),
                          LINEAR_STEP_FLOOR * (np.asarray(hi, dtype=float)
                                               - np.asarray(lo, dtype=float)))


@dataclass(frozen=True, eq=False)
class _Spec:
    """Everything an evaluation needs except the model (picklable: sent to
    process workers through the initializer).

    ``targets[j]`` is ``None`` for an observable scale/offset (a virtual
    finite-difference parameter); ``applies[j]`` is ``None`` (every
    experiment) or ``(i,)`` for a local parameter.
    """
    names: tuple
    targets: tuple
    log: np.ndarray
    applies: tuple
    lo: np.ndarray
    hi: np.ndarray
    conditions: tuple
    blocks: tuple
    data: np.ndarray
    sigma: np.ndarray
    counters: tuple
    rel_step: float
    max_refinements: int
    loss: str
    f_scale: float

    @property
    def n_data(self):
        return self.data.size

    def natural(self, x):
        """Natural values of the transformed point ``x``."""
        x = np.asarray(x, dtype=float)
        with np.errstate(all='ignore'):
            return np.where(self.log, np.exp(x), x)

    def fd_parameters(self, x):
        """The :class:`_fd.FDParameter` objects of point ``x``: the step is
        ``rel_step`` in log space, and ``rel_step * max(|p|,
        LINEAR_STEP_FLOOR * (hi - lo))`` (:func:`_linear_scale`) in linear
        space."""
        out = []
        for j, value in enumerate(self.natural(x)):
            log = bool(self.log[j])
            if log:
                step = self.rel_step
            else:
                step = self.rel_step * float(
                    _linear_scale(value, self.lo[j], self.hi[j]))
            out.append(_fd.FDParameter(self.names[j], self.targets[j],
                                       float(value), log, step,
                                       self.applies[j]))
        return out

    def predictions(self, results, values):
        """Model values at every data point (residual order)."""
        out = np.empty(self.n_data)
        for result, blocks in zip(results, self.blocks):
            for b in blocks:
                scale, offset = (values[v] if isinstance(v, str) else v
                                 for v in (b.scale, b.offset))
                y = b.observable._apply(result[b.observable.selection],
                                        scale, offset)
                out[b.rows] = y[b.grid]
        return out

    def cost(self, r):
        """``0.5 * sum(rho(r**2))`` with ``least_squares``'s ``loss`` and
        ``f_scale`` (``rho`` scaled as
        ``f_scale**2 * rho((r/f_scale)**2)``)."""
        with np.errstate(all='ignore'):
            z = (np.asarray(r, dtype=float) / self.f_scale) ** 2
            if self.loss == 'linear':
                rho = z
            elif self.loss == 'soft_l1':
                rho = 2.0 * (np.sqrt(1.0 + z) - 1.0)
            elif self.loss == 'huber':
                rho = np.where(z <= 1.0, z, 2.0 * np.sqrt(z) - 1.0)
            else:
                rho = np.log1p(z)
            return float(0.5 * self.f_scale ** 2 * np.sum(rho))


def _residual_vector(spec, results, values):
    """``(model - data) / sigma`` at every data point (the ``fn`` of the
    finite-difference engine)."""
    with np.errstate(all='ignore'):
        return (spec.predictions(results, values) - spec.data) / spec.sigma


class _Evaluator:
    """Residuals and their Jacobian on one (private) model, with failures
    turned into a constant penalty vector and counted.

    The last point and the last Jacobian are cached, so ``least_squares``'s
    ``fun(x)`` and ``jac(x)`` at the same ``x`` share the nominal runs. The
    cached residual array (also the nominal ``f`` of the Jacobian) is
    read-only and never handed out: :meth:`residuals` returns a copy, as
    ``least_squares`` may rescale the array it gets in place (SciPy 1.15
    does for robust losses).
    """

    def __init__(self, model, spec, penalty=PENALTY_RESIDUAL):
        self.model = model
        self.spec = spec
        self.penalty = float(penalty)
        self.fn = functools.partial(_residual_vector, spec)
        self.n_failed = 0
        self.error = None
        self._point = None      # (key, FDPoint | None, r, ok)
        self._jac = None        # (key, FDJacobian)

    def evaluate(self, x):
        """``(r, ok)`` at transformed ``x``: the weighted residuals, or the
        penalty vector (``ok=False``) if a run failed or a residual is not
        finite."""
        x = np.asarray(x, dtype=float)
        key = x.tobytes()
        if self._point is not None and self._point[0] == key:
            return self._point[2], self._point[3]
        point = _fd.evaluate(self.model, self.spec.conditions,
                             self.spec.fd_parameters(x), self.fn,
                             event_counters=self.spec.counters)
        ok = point.ok and bool(np.all(np.isfinite(point.f)))
        if ok:
            r = point.f
        else:
            self.n_failed += 1
            self.error = point.error or 'non-finite model values'
            r = np.full(self.spec.n_data, self.penalty)
        r.flags.writeable = False   # shared by the cache and the Jacobian
        self._point = (key, point if ok else None, r, ok)
        return r, ok

    def residuals(self, x):
        """A copy of the residual vector at ``x`` (``least_squares``'s
        ``fun``; the caller may modify it)."""
        return self.evaluate(x)[0].copy()

    def fd_jacobian(self, x):
        """The :class:`_fd.FDJacobian` at ``x`` (``None`` where the point
        fails), reusing the nominal runs of :meth:`evaluate`."""
        x = np.asarray(x, dtype=float)
        key = x.tobytes()
        if self._jac is not None and self._jac[0] == key:
            return self._jac[1]
        _, ok = self.evaluate(x)
        if not ok:
            return None
        jac = _fd.fd_jacobian(self.model, self.spec.conditions,
                              self.spec.fd_parameters(x), self.fn,
                              method='central',
                              event_counters=self.spec.counters,
                              max_refinements=self.spec.max_refinements,
                              nominal=self._point[1])
        self._jac = (key, jac)
        return jac

    def jacobian(self, x):
        """``dr/dx`` at ``x`` for the optimizer: non-finite entries (a column
        whose differences failed) are 0, and a failed point gives zeros."""
        jac = self.fd_jacobian(x)
        if jac is None:
            return np.zeros((self.spec.n_data, len(self.spec.names)))
        return np.where(np.isfinite(jac.jacobian), jac.jacobian, 0.0)

    def predictions(self):
        """Model values at the data points of the last successful point."""
        point = self._point[1]
        return self.spec.predictions(point.results, point.values)


class CalibrationProblem:
    """A weighted least-squares calibration problem.

    Parameters
    ----------
    model : KineticModel
        The model; never modified (evaluations run on a private copy taken
        here, so later edits of ``model`` do not affect the problem).
    experiments : sequence of Experiment
        Uniquely named experiments.
    observables : sequence of Observable
        Uniquely named observables shared by the experiments (an
        experiment may instead list its own complete set: problem
        observables by name and its own :class:`Observable` objects; see
        :class:`Experiment`).
    parameters : sequence of FitParameter
        What to estimate (at least one; each target once).
    loss : {'linear', 'soft_l1', 'huber', 'cauchy'}
        Loss of ``scipy.optimize.least_squares`` (``'linear'``: plain
        least squares; the others down-weight outliers and need absolute
        sigmas, ``absolute_sigma=True``).
    f_scale : float
        Residual size (in units of sigma) at which a robust loss starts to
        down-weight. Defaults to ``1.0``.
    absolute_sigma : bool
        ``True`` (default): the sigmas are absolute and the covariance is
        ``(J^T J)^-1`` (see :attr:`FitResult.covariance`). ``False``: they
        are weights only (the noise level is estimated from the residuals)
        and the covariance is scaled by the reduced chi-square; only with
        ``loss='linear'``, whose fit does not depend on the scale of the
        weights (a robust loss's ``f_scale`` would refer to that arbitrary
        scale).
    relative_floor : float
        Floor of proportional noise, as a fraction of each observable's
        largest measured magnitude per experiment (default ``1e-3``, the
        ``RELATIVE_FLOOR`` of
        :func:`~nskinetics.analysis.identifiability`), in ``(0, 1]``.
    rel_step : float
        Finite-difference step in log space (relative), in ``(0, 1]``;
        linear-space parameters use ``rel_step * max(|p|, 1e-3 (hi - lo))``.
        Defaults to ``1e-3``. The floor guards a ``log=False`` parameter
        near 0 against integration noise, but ties its step to the width
        of the bounds: for ``|p|`` more than about ``1e4``-``1e5`` times
        smaller than ``hi - lo`` the step is not small relative to ``p``,
        and the Jacobian, the estimate and the SE are inaccurate without
        any flag (see :class:`FitParameter` ``log``). A smaller
        ``rel_step`` (like a looser ``rtol``, see ``tolerances``) narrows
        the margin of the rank test and column filter of
        :attr:`FitResult.covariance` over the integration noise, which is
        thin at the defaults on fed-batch and other event models, so a
        parameter the outputs do not depend on can show as
        ``high_rel_se`` or ``wide_interval`` instead of
        ``rank_deficient`` (see :ref:`calibration-rank`).
    tolerances : (float, float), optional
        ``(atol, rtol)`` of the private copy's integrator. ``None``
        (default) tightens them to at most ``(1e-12, 1e-10)`` (and the
        absolute tolerance further to outputs that start at 0), as
        :func:`~nskinetics.analysis.local_sensitivity` does. A looser
        ``rtol`` raises the integration noise of every column, as a
        smaller ``rel_step`` does, so a parameter the outputs do not
        depend on can escape the zeroing in :attr:`FitResult.covariance`
        and show as ``high_rel_se`` or ``wide_interval``. A tighter one
        (e.g. ``(1e-14, 1e-12)``) widens that margin, which is thin on
        fed-batch and other event models and when no dominant rate
        constant is fitted (see :ref:`calibration-rank`).
    event_counters : sequence of str, optional
        Event counters whose firings the Jacobian watches; ``None``
        (default) detects them (e.g. a FeedSpike's ``count_var``).
    max_refinements : int
        Step refinements when a perturbation changes an event count
        (default 2).

    Attributes
    ----------
    model, experiments, observables, parameters, loss, f_scale,
    absolute_sigma, relative_floor, rel_step, tolerances, max_refinements
        As given (sequences as tuples).
    event_counters : tuple of str
        The event counters watched (as given, or detected).
    names : list of str
        Free-parameter names: a global parameter's target, ``'target@exp'``
        for each experiment of a local one (a local observable scale only
        for the experiments that use it).
    n_free, n_data : int
        Numbers of free parameters and data points.

    Raises
    ------
    TypeError
        For arguments of the wrong type.
    ValueError
        Before any simulation: no parameters or experiments, duplicate
        experiment, observable or target names, unknown targets (model ids
        are checked like :meth:`~nskinetics.KineticModel.run` checks a
        condition), observables whose selection the model cannot record,
        data of an observable an experiment does not define, a named scale
        or offset without a :class:`FitParameter`, a :class:`FitParameter`
        for the scale or offset of an observable that has no data in any
        experiment (it would affect no residual), a data point without a
        usable sigma, more free parameters than data points, experiment
        conditions the targets cannot be merged into, or a robust loss with
        ``absolute_sigma=False``.

    Warns
    -----
    UserWarning
        With more than 15 free parameters (fit an identifiable subset;
        see :func:`~nskinetics.analysis.identifiability`), for a global
        target that an experiment's condition sets (its one fitted value
        replaces the per-experiment values), for targets read by SBML
        initial assignments, and for event counters that cannot be
        watched.

    Notes
    -----
    Residuals are ``(model_i - data_i) / sigma_i`` with fixed, data-based
    sigmas (see :class:`Observable` ``noise`` and
    :ref:`calibration-noise`). Each evaluation runs every
    experiment through :meth:`~nskinetics.KineticModel.run` from its
    baseline, with the fitted values merged into the experiment's own
    condition (:meth:`~nskinetics.Condition.with_parameters` /
    :meth:`~nskinetics.Condition.with_initial`; local values only into
    their experiment's), so per-experiment settings never leak between
    experiments. A global target replaces
    what the conditions set for it (with a warning); fit values that
    differ between experiments with ``scope='local'``.

    Data times must be in model time units.
    """

    def __init__(self, model, experiments, observables, parameters, *,
                 loss='linear', f_scale=1.0, absolute_sigma=True,
                 relative_floor=RELATIVE_FLOOR, rel_step=1e-3,
                 tolerances=None, event_counters=None, max_refinements=2):
        if not isinstance(model, KineticModel):
            raise TypeError(f'model must be a KineticModel; got '
                            f'{type(model).__name__}.')
        experiments = self._sequence('experiments', experiments, Experiment)
        observables = self._sequence('observables', observables, Observable,
                                     allow_empty=True)
        parameters = self._sequence('parameters', parameters, FitParameter)
        if loss not in LOSSES:
            raise ValueError(f'loss must be one of {LOSSES}; got {loss!r}.')
        f_scale = _positive('f_scale', f_scale)
        if not isinstance(absolute_sigma, (bool, np.bool_)):
            raise TypeError(f'absolute_sigma must be a bool; got '
                            f'{absolute_sigma!r}.')
        if loss != 'linear' and not absolute_sigma:
            raise ValueError(
                f'loss={loss!r} needs absolute sigmas (absolute_sigma=True): '
                'a robust loss down-weights residuals beyond f_scale standard '
                'deviations; with weights-only sigmas (absolute_sigma=False) '
                'f_scale would refer to their arbitrary scale, so the '
                'estimate and the profile would change when every sigma is '
                'rescaled. Give absolute sigmas, or use loss=\'linear\' with '
                'absolute_sigma=False.')
        relative_floor = _check_relative_floor(relative_floor)
        rel_step = _fd.check_step(rel_step)
        tolerances = _fd.check_tolerances(tolerances)
        max_refinements = _count('max_refinements', max_refinements, 0)
        counters = (None if event_counters is None
                    else _fd.check_counters(event_counters))
        self._check_names('experiment', [e.name for e in experiments])
        self._check_names('observable', [o.name for o in observables])
        self._check_names('target', [p.target for p in parameters])

        # Observables of each experiment, in the order its data name them.
        shared = {o.name: o for o in observables}
        used = [self._resolve(e, shared) for e in experiments]
        obs_params = {}
        for obs in (o for u in used for o in u.values()):
            for name in obs.parameters:
                obs_params.setdefault(name, []).append(obs.name)
        for obs in {id(o): o for u in used for o in u.values()}.values():
            _check_ids(model, f'observable {obs.name!r} selection',
                       [obs.selection], settable=False)
        by_target = {p.target: p for p in parameters}
        unfitted = sorted(n for n in obs_params if n not in by_target)
        if unfitted:
            raise ValueError(
                f'Observable scale/offset name(s) {unfitted} have no '
                'FitParameter; add FitParameter(name, bounds) for each, or '
                'give the scale/offset as a number.')
        keys = model._model_keys()
        clash = sorted(n for n in obs_params if n in keys)
        if clash:
            raise ValueError(f'Observable scale/offset name(s) {clash} are '
                             'also model ids; rename them.')
        # A scale/offset of an observable no experiment has data for is in
        # no residual (and not a model id): say so, rather than "unknown id"
        defined = {}
        own = (o for e in experiments for o in (e.observables or ())
               if isinstance(o, Observable))
        for obs in (*observables, *own):
            for name in obs.parameters:
                defined.setdefault(name, set()).add(obs.name)
        orphan = sorted(p.target for p in parameters
                        if p.target in defined and p.target not in obs_params
                        and p.target not in keys)
        if orphan:
            owners = sorted({o for n in orphan for o in defined[n]})
            raise ValueError(
                f'FitParameter target(s) {orphan} are the scale/offset of '
                f'observable(s) {owners} that have no data in any experiment, '
                'so they affect no residual; remove them or add data.')
        model_targets = [p.target for p in parameters
                         if p.target not in obs_params]
        _check_ids(model, 'parameters', model_targets)
        if counters:
            _check_ids(model, 'event_counters', counters, settable=False)

        # Free parameters: names, transformed bounds, start values.
        names, targets, applies, log, kinds, exps, target_of = (
            [], [], [], [], [], [], [])
        lo, hi, x_init, nat_lo, nat_hi = [], [], [], [], []
        for p in parameters:
            virtual = p.target in obs_params
            if p.scope == 'global':
                scopes = [(None, None)]
            else:
                scopes = [(i, e) for i, e in enumerate(experiments)
                          if not virtual or any(
                              p.target in o.parameters
                              for o in used[i].values())]
            rng = p.range
            a, b = rng.transformed_bounds
            for i, e in scopes:
                names.append(p.target if e is None else f'{p.target}@{e.name}')
                targets.append(None if virtual else p.target)
                applies.append(None if i is None else (i,))
                kinds.append('observable' if virtual else 'model')
                target_of.append(p.target)
                exps.append(None if e is None else e.name)
                log.append(p.log)
                lo.append(a)
                hi.append(b)
                nat_lo.append(p.bounds[0])     # exact, not exp(ln(lo))
                nat_hi.append(p.bounds[1])
                start = p.initial_for(None if e is None else e.name)
                # None: resolved on first use of x0 (see _start)
                x_init.append(None if start is None
                              else rng.to_transformed(start))
        self._check_names('free parameter', names)

        # Run conditions, and every target merged into them (validates the
        # merge, the ids, the selections and the reset semantics; no runs).
        conditions = []
        for i, (e, obs) in enumerate(zip(experiments, used)):
            sels = tuple(dict.fromkeys(o.selection for o in obs.values()))
            conditions.append(e.condition.replace(
                times=tuple(float(t) for t in e.times), selections=sels,
                label=e.name))
        pinned = []
        for i, c in enumerate(conditions):
            mine = {targets[j]: 1.0 for j in range(len(names))
                    if targets[j] is not None
                    and (applies[j] is None or applies[j] == (i,))}
            try:
                pinned.append(_fd.pin(model, c, mine))
            except ValueError as e:
                raise ValueError(
                    f'Experiment {experiments[i].name!r}: the targets cannot '
                    f'be merged into its condition ({e}); give each species '
                    "in the same form (amount id or '[concentration]') as "
                    "the condition's own initial values.") from e
        _ex._validate_batch(model, pinned, None, 1, True, 'record')

        # Data, sigmas and the residual layout.
        blocks, rows_data, rows_sigma, frame_rows = [], [], [], []
        offset = 0
        for i, (e, obs) in enumerate(zip(experiments, used)):
            frame = e.data.frame
            y = frame['value'].to_numpy(dtype=float)
            s = frame['sigma'].to_numpy(dtype=float)
            t = frame['time'].to_numpy(dtype=float)
            grid = np.searchsorted(e.times, t)
            sigma = np.empty(len(frame))
            mine = []
            for name, o in obs.items():
                rows = np.flatnonzero((frame['observable'] == name).to_numpy())
                if not rows.size:
                    continue
                given = s[rows]
                model_sigma = o.noise_sigma(y[rows],
                                            relative_floor=relative_floor)
                sig = np.where(np.isfinite(given), given, model_sigma)
                with np.errstate(invalid='ignore'):
                    bad = ~(np.isfinite(sig) & (sig > 0))
                if bad.any():
                    hint = ("give Observable(sigma=...) or a sigma column"
                            + (" (e.g. sigma=1.0 with absolute_sigma=False "
                               "to estimate the noise level from the "
                               "residuals)" if loss == 'linear' else
                               ' (absolute standard deviations: a robust '
                               'loss needs them)')
                            if o.noise == 'additive' else
                            'its data are all 0; give a sigma column or '
                            "use noise='combined'")
                    raise ValueError(
                        f'Experiment {e.name!r}, observable {name!r}: no '
                        f'usable sigma for data rows {_rows(bad)}; {hint}.')
                sigma[rows] = sig
                scale, off = (self._free_name(v, e.name, by_target)
                              for v in (o.scale, o.offset))
                mine.append(_Block(observable=o, rows=rows + offset,
                                   grid=grid[rows], scale=scale, offset=off))
            blocks.append(tuple(mine))
            rows_data.append(y)
            rows_sigma.append(sigma)
            frame_rows.append(pd.DataFrame({
                'experiment': pd.Series(list(frame['experiment']),
                                        dtype=object),
                'observable': pd.Series(list(frame['observable']),
                                        dtype=object),
                'time': t, 'data': y, 'sigma': sigma}))
            offset += len(frame)
        n_data, n_free = offset, len(names)
        if n_free > n_data:
            raise ValueError(f'{n_free} free parameters {names} exceed the '
                             f'{n_data} data points; fix some parameters or '
                             'add data.')

        # Warnings, then the private model (a copy; no runs here).
        if n_free > N_FREE_WARN:
            warnings.warn(
                f'{n_free} free parameters: typical data identify only a few. '
                'Run nskinetics.analysis.local_sensitivity and '
                'identifiability first and fit an identifiable subset '
                '(IdentifiabilityResult.best_subset()).', UserWarning,
                stacklevel=2)
        for p in parameters:
            if p.scope != 'global' or p.target in obs_params:
                continue
            setting = [e.name for e in experiments
                       if p.target in e.condition.parameters
                       or p.target in e.condition.initial]
            if setting:
                warnings.warn(
                    f'Global FitParameter {p.target!r}: its single fitted '
                    'value replaces the value set by the condition(s) of '
                    f'experiment(s) {setting}; use scope=\'local\' to fit one '
                    'value per experiment, or remove it from those '
                    'conditions.', UserWarning, stacklevel=2)
        _fd.warn_initial_assignment_inputs(model, model_targets, stacklevel=3)
        private = _fd.private_model(model, tolerances=tolerances,
                                    rel_step=rel_step, stacklevel=3)
        if counters is None:
            counters = private._fd_counters
            _fd.warn_unwatched_counters(private._fd_unwatched, stacklevel=3)

        self.model = model
        self.experiments = experiments
        self.observables = observables
        self.parameters = parameters
        self.loss = loss
        self.f_scale = f_scale
        self.absolute_sigma = bool(absolute_sigma)
        self.relative_floor = relative_floor
        self.rel_step = rel_step
        self.tolerances = tolerances
        self.max_refinements = max_refinements
        self.event_counters = tuple(counters)
        self.names = names
        self.n_free = n_free
        self.n_data = n_data
        self._kinds = kinds
        self._target_of = target_of
        self._experiment_of = exps
        lo, hi = np.array(lo, dtype=float), np.array(hi, dtype=float)
        self._natural_lo = np.array(nat_lo, dtype=float)
        self._natural_hi = np.array(nat_hi, dtype=float)
        self._spec = _Spec(
            names=tuple(names), targets=tuple(targets),
            log=np.array(log, dtype=bool), applies=tuple(applies), lo=lo,
            hi=hi, conditions=tuple(conditions), blocks=tuple(blocks),
            data=np.concatenate(rows_data), sigma=np.concatenate(rows_sigma),
            counters=tuple(counters), rel_step=rel_step,
            max_refinements=max_refinements, loss=loss, f_scale=f_scale)
        self._x_initial = tuple(x_init)
        self._x_start = None
        self._data = pd.concat(frame_rows, ignore_index=True)
        self._private = private
        self._prepared = tolerances is not None
        self._evaluator = _Evaluator(private, self._spec)
        self._units = dict(model._units)

    # --- construction helpers ------------------------------------------------
    @staticmethod
    def _sequence(arg, items, kind, allow_empty=False):
        if isinstance(items, (kind, str, bytes, Mapping)) or not hasattr(
                items, '__iter__'):
            raise TypeError(f'{arg} must be a sequence of {kind.__name__} '
                            'objects.')
        items = tuple(items)
        for i, item in enumerate(items):
            if not isinstance(item, kind):
                raise TypeError(f'{arg}[{i}] is a {type(item).__name__}, not '
                                f'a {kind.__name__}.')
        if not items and not allow_empty:
            raise ValueError(f'{arg} must not be empty.')
        return items

    @staticmethod
    def _check_names(what, names):
        dup = sorted({n for n in names if names.count(n) > 1})
        if dup:
            raise ValueError(f'Duplicate {what} name(s) {dup}.')

    @staticmethod
    def _resolve(experiment, shared):
        """``{name: Observable}`` of the observables an experiment's data
        name, in data order."""
        if experiment.observables is None:
            available = dict(shared)
        else:
            available = {}
            for o in experiment.observables:
                if isinstance(o, Observable):
                    available[o.name] = o
                elif o in shared:
                    available[o] = shared[o]
                else:
                    raise ValueError(f'Experiment {experiment.name!r} names '
                                     f'observable {o!r}, which is not among '
                                     'the problem observables '
                                     f'{list(shared)}.')
        missing = [o for o in experiment.data.observables
                   if o not in available]
        if missing:
            raise ValueError(f'Experiment {experiment.name!r} has data for '
                             f'observable(s) {missing} that it does not '
                             f'define; defined: {list(available)} (an '
                             "experiment's observables list is its complete "
                             'set; name the shared observables it also '
                             'uses).')
        return {o: available[o] for o in experiment.data.observables}

    @staticmethod
    def _free_name(value, experiment, by_target):
        """A scale/offset as a number or the free-parameter name holding it
        in ``experiment``."""
        if not isinstance(value, str):
            return value
        return (value if by_target[value].scope == 'global'
                else f'{value}@{experiment}')

    # --- public API ----------------------------------------------------------
    @property
    def log(self):
        """Whether each free parameter is fitted in log space
        (``numpy.ndarray`` of bool)."""
        return self._spec.log.copy()

    @property
    def bounds(self):
        """``(lo, hi)``: the bounds in the fitting space (``ln`` for log
        parameters), as two float arrays."""
        return self._spec.lo.copy(), self._spec.hi.copy()

    @property
    def x0(self):
        """The first start (fitting space, ``numpy.ndarray``).

        Each parameter's ``initial`` value; by default the value the run
        starts from (the experiment condition's entry if it sets the
        target, else the model's value; see :class:`FitParameter`), or the
        midpoint of the bounds in the fitting space. The point is clipped
        into the bounds (a value on a bound is kept). The first access may
        run each experiment once (to read model values at its start time).
        """
        return _clip(self._start(), self._spec.lo, self._spec.hi)

    def _start(self):
        """The first start before it is clipped into the bounds (cached).

        A default (``None``) start value of a model target is read from
        the experiment it applies to (a global one: the first experiment)
        with :func:`_fd.nominal_values` -- the condition's own entry, or
        the first row of one run on the private model. It is used if it is
        finite and within the bounds; otherwise, if that run fails, and for
        observable scales/offsets, the start is the midpoint of the bounds
        in the fitting space.
        """
        if self._x_start is not None:
            return self._x_start
        spec = self._spec
        x = list(self._x_initial)
        pending = {}    # experiment -> free parameters read from its run
        for j, value in enumerate(x):
            if value is None and spec.targets[j] is not None:
                i = 0 if spec.applies[j] is None else spec.applies[j][0]
                pending.setdefault(i, []).append(j)
        if pending:
            private = self._ready()
            for i, js in pending.items():
                targets = list(dict.fromkeys(spec.targets[j] for j in js))
                try:
                    values = _fd.nominal_values(private, spec.conditions[i],
                                                targets)
                except KineticSimulationError:
                    continue        # the midpoint; starts never raise
                for j in js:
                    v = values[spec.targets[j]]
                    if (math.isfinite(v) and self._natural_lo[j] <= v
                            <= self._natural_hi[j]
                            and (v > 0 or not spec.log[j])):
                        x[j] = math.log(v) if spec.log[j] else v
        self._x_start = np.array(
            [0.5 * (spec.lo[j] + spec.hi[j]) if v is None else v
             for j, v in enumerate(x)], dtype=float)
        return self._x_start

    @property
    def data(self):
        """The data points in residual order (``pandas.DataFrame`` with
        columns ``experiment, observable, time, data, sigma``; sigma as
        used, in observable units)."""
        return self._data.copy()

    @property
    def n_failed_evals(self):
        """Failed evaluations of :meth:`residuals` so far (int)."""
        return self._evaluator.n_failed

    def to_natural(self, x):
        """Natural values of a point in the fitting space.

        Parameters
        ----------
        x : array_like or mapping
            The point, as for :meth:`residuals`.

        Returns
        -------
        dict
            ``{name: value}`` in natural (model or observable) units.

        Raises
        ------
        ValueError
            For an invalid point.

        Warns
        -----
        UserWarning
            As :meth:`residuals`.
        """
        return dict(zip(self.names, map(float, self._spec.natural(
            self._as_x(x)))))

    def _as_x(self, theta, stacklevel=3, arg='theta'):
        """``theta`` (an array in the fitting space, or ``{name: natural
        value}``) as a transformed float array; an array outside the bounds
        warns (``stacklevel`` is that of the warning; ``arg`` names
        ``theta`` in messages)."""
        if isinstance(theta, Mapping):
            unknown = [k for k in theta if k not in self.names]
            missing = [n for n in self.names if n not in theta]
            if unknown or missing:
                raise ValueError(f'{arg} must give every free parameter '
                                 f'{self.names}; unknown {unknown}, missing '
                                 f'{missing}.')
            x = []
            for j, name in enumerate(self.names):
                v = _finite(f'{arg}[{name!r}]', theta[name])
                if self._spec.log[j]:
                    if v <= 0:
                        raise ValueError(f'{arg}[{name!r}] = {v!r} must be '
                                         '> 0 (fitted in log space).')
                    v = math.log(v)
                x.append(v)
            return np.array(x, dtype=float)
        x = np.asarray(theta, dtype=float)
        if x.shape != (self.n_free,):
            raise ValueError(f'{arg} must have shape ({self.n_free},); got '
                             f'{x.shape}.')
        if not np.all(np.isfinite(x)):
            raise ValueError(f'{arg} must be finite.')
        spec = self._spec
        outside = [n for n, v, a, b in zip(self.names, x, spec.lo, spec.hi)
                   if not a <= v <= b]
        if outside:
            # the optimizer and profiles never leave the bounds: most likely
            # natural values given as an array
            warnings.warn(
                f'{arg} lies outside the bounds for {outside}: an array is '
                'read in the fitting space (ln p for log parameters; see '
                'CalibrationProblem.bounds). Pass {name: natural value} to '
                'give natural values.', UserWarning, stacklevel=stacklevel)
        return x

    def _ready(self):
        """Prepare the private model once before the first evaluation:
        without explicit tolerances, tighten its absolute tolerance to the
        outputs that start at 0 (one run per experiment)."""
        if not self._prepared:
            sels = [o.observable.selection for blocks in self._spec.blocks
                    for o in blocks]
            _fd.scale_atol_to_outputs(self._private, self._spec.conditions,
                                      sels)
            self._prepared = True
        return self._private

    def residuals(self, theta):
        """Weighted residuals ``(model - data) / sigma`` at a point.

        Parameters
        ----------
        theta : array_like or mapping
            The point: an array in the fitting space (``ln p`` for log
            parameters, in the order of :attr:`names`), or ``{name:
            natural value}`` for every free parameter.

        Returns
        -------
        numpy.ndarray
            ``(n_data,)``, in the order of :attr:`data`. If a run fails (or a
            model value is not finite), every entry is the penalty
            ``PENALTY_RESIDUAL`` (``1e3``) and :attr:`n_failed_evals` is
            incremented; nothing is raised.

        Raises
        ------
        ValueError
            For an invalid ``theta``.

        Warns
        -----
        UserWarning
            If an array ``theta`` lies outside :attr:`bounds` (fitting
            space): natural values must be given as a mapping.
        """
        return self._residuals_at(self._as_x(theta))

    def _residuals_at(self, x):
        self._ready()
        return self._evaluator.evaluate(x)[0].copy()

    def cost(self, theta):
        """The objective ``0.5 * sum(rho(r**2))`` at a point.

        Parameters
        ----------
        theta : array_like or mapping
            The point, as for :meth:`residuals`.

        Returns
        -------
        float
            The objective with the problem's ``loss`` and ``f_scale`` (for
            the linear loss half the weighted sum of squared residuals).
            A failed evaluation gives the cost of the penalty residuals.

        Raises
        ------
        ValueError
            For an invalid ``theta``.

        Warns
        -----
        UserWarning
            As :meth:`residuals`.
        """
        return self._spec.cost(self._residuals_at(self._as_x(theta)))

    def __repr__(self):
        return (f'CalibrationProblem({len(self.experiments)} experiment(s), '
                f'{self.n_data} data points, free parameters {self.names}, '
                f'loss={self.loss!r})')


# --- fitting -----------------------------------------------------------------
def _clip(x, lo, hi):
    """``x`` clipped into the bounds ``[lo, hi]`` (fitting space)."""
    return np.clip(np.asarray(x, dtype=float), lo, hi)


def _run_start(evaluator, x0, options):
    """One ``least_squares`` run from ``x0`` (fitting space); returns a
    plain-data record (also from a process worker)."""
    spec = evaluator.spec
    # least_squares works on u = x - c. SciPy's trf sets its first trust
    # radius from the scaled norm of the start, so a start near 0 in the
    # fitting space (ln p near 0, i.e. p near 1, as a log-midpoint of
    # bounds symmetric about 1 gives, or a linear start at 0) would stall
    # there after two evaluations and report convergence. It also stops
    # when norm(step) < xtol * (xtol + norm(u)), so the shift must stay
    # O(1): it must not grow with the bounds or the start, or one wide
    # range would loosen the stopping test of every parameter. Coordinates
    # with |x0| < 1 therefore start at +-1; all others are left unshifted
    # (plain SciPy). The shift is constant, so the Jacobian is unchanged.
    # x0 is rounded through u so the optimizer's first call hits the
    # evaluator's cached point; the record keeps the start as given.
    start = np.asarray(x0, dtype=float)
    c = np.where(np.abs(start) < 1.0,
                 start - np.where(start < 0, -1.0, 1.0), 0.0)
    u0 = start - c
    x0 = u0 + c
    record = {'x0': start, 'x': np.full(x0.size, np.nan), 'ok': False,
              'success': False, 'cost': math.nan, 'chi2': math.nan,
              'nfev': 1, 'njev': 0, 'active_mask': np.zeros(x0.size, int),
              'message': '', 'r': None, 'model': None, 'jacobian': None,
              'fd_flags': None}
    r0, ok0 = evaluator.evaluate(x0)
    if not ok0:
        record['message'] = f'The start point failed: {evaluator.error}'
        record['n_failed'] = evaluator.n_failed
        return record
    with np.errstate(all='ignore'):
        evaluator.penalty = max(PENALTY_RESIDUAL,
                                10.0 * float(np.max(np.abs(r0))))
    res = least_squares(lambda u: evaluator.residuals(u + c), u0,
                        jac=lambda u: evaluator.jacobian(u + c),
                        bounds=(spec.lo - c, spec.hi - c), method='trf',
                        loss=spec.loss, f_scale=spec.f_scale, x_scale='jac',
                        max_nfev=options.get('max_nfev'))
    # u + c can overshoot a bound by an ulp
    x = np.clip(np.asarray(res.x, dtype=float) + c, spec.lo, spec.hi)
    record.update(x=x, nfev=int(res.nfev), njev=int(res.njev or 0),
                  message=str(res.message),
                  active_mask=np.asarray(res.active_mask, dtype=int))
    r, ok = evaluator.evaluate(x)
    if not ok:
        record['message'] = (f'{res.message} The final point failed: '
                             f'{evaluator.error}')
        record['n_failed'] = evaluator.n_failed
        return record
    jac = evaluator.fd_jacobian(x)
    record.update(ok=True, success=bool(res.status > 0),
                  cost=spec.cost(r), chi2=float(np.sum(r * r)), r=r.copy(),
                  model=evaluator.predictions(), jacobian=jac.jacobian.copy(),
                  fd_flags=jac.flags, n_failed=evaluator.n_failed)
    return record


def _start_task(model, x0, *, spec, options):
    """One start on a worker's model (the function of the process pool)."""
    return _run_start(_Evaluator(model, spec), x0, options)


def _problem_observables(problem):
    """The problem's shared observables and the experiments' own ones."""
    return list(problem.observables) + [
        o for e in problem.experiments for o in (e.observables or ())
        if isinstance(o, Observable)]


def _check_transforms_picklable(problem):
    """``spawn_check`` of the process pools of :func:`fit_parameters` and
    :func:`profile_likelihood`: name an observable whose transform cannot
    be pickled, ahead of the generic error that would print the whole
    problem description."""
    for obs in _problem_observables(problem):
        if obs.transform is None:
            continue
        try:
            pickle.dumps(obs.transform)
        except (pickle.PicklingError, TypeError, AttributeError) as e:
            raise ValueError(
                f'Observable {obs.name!r}: transform {obs.transform!r} is '
                f'not picklable ({e}); with workers > 1 under the spawn '
                'start method (macOS, Windows) it must be a module-level '
                'function such as np.log10 (not a lambda or a nested '
                'function), and scripts must call fit_parameters/'
                "profile_likelihood under `if __name__ == '__main__':`."
            ) from e


def _check_fit_args(problem, method, n_starts, sampler, workers, max_nfev,
                    global_options, initial):
    if not isinstance(problem, CalibrationProblem):
        raise TypeError(f'problem must be a CalibrationProblem; got '
                        f'{type(problem).__name__}.')
    if method not in FIT_METHODS:
        raise ValueError(f'method must be one of {FIT_METHODS}; got '
                         f'{method!r}.')
    n_starts = _count('n_starts', n_starts, 1)
    if sampler not in SAMPLERS:
        raise ValueError(f'sampler must be one of {SAMPLERS}; got '
                         f'{sampler!r}.')
    workers = _ex._check_workers(workers)
    if max_nfev is not None:
        max_nfev = _count('max_nfev', max_nfev, 1)
    if method == 'hybrid':
        if workers > 1:
            raise ValueError("workers > 1 parallelizes the starts of "
                             "method='multistart'; method='hybrid' runs "
                             'serially (parallel=\'none\').')
        if global_options is not None and not isinstance(global_options,
                                                         Mapping):
            raise TypeError('global_options must be a mapping.')
        options = dict(global_options or {})
        reserved = [k for k in options if k in _RESERVED_OPTIONS]
        if reserved:
            raise ValueError(f'global_options cannot set {reserved}.')
        unknown = sorted(str(k) for k in options if k not in _GLOBAL_OPTIONS)
        if unknown:
            raise ValueError(f'Unknown global_options key(s) {unknown}; '
                             f'allowed: {list(_GLOBAL_OPTIONS)}.')
        for key in _COUNT_OPTIONS:
            if key in options:
                _count(f'global_options[{key!r}]', options[key], 1)
        if initial is not None:
            raise ValueError("initial sets start 0 of method='multistart'; "
                             "method='hybrid' starts from its global stage.")
    elif global_options is not None:
        raise ValueError("global_options is only used by method='hybrid'.")
    x_init = None
    if initial is not None:
        if isinstance(initial, FitResult):
            initial = initial.estimates.to_dict()
        # user -> fit_parameters -> _check_fit_args -> _as_x
        x_init = problem._as_x(initial, stacklevel=4, arg='initial')
    return n_starts, workers, max_nfev, x_init


def fit_parameters(problem, *, method='multistart', n_starts=10,
                   sampler='sobol', seed=None, workers=1, max_nfev=None,
                   global_options=None, initial=None):
    """Estimate the free parameters of a calibration problem.

    Parameters
    ----------
    problem : CalibrationProblem
        The problem to solve (not modified; its private model is used).
    method : {'multistart', 'hybrid'}
        ``'multistart'`` (default): ``scipy.optimize.least_squares``
        (``method='trf'`` with the bounds, the problem's ``loss`` and
        ``f_scale``, ``x_scale='jac'`` and the event-aware central
        finite-difference Jacobian) from ``n_starts`` points in the fitting
        space. ``'hybrid'``: a global stage,
        :func:`~nskinetics.utils.hybrid_global_optimize` (serial,
        ``parallel='none'``) on the cost ``0.5 * sum(rho(r**2))``, then one
        ``least_squares`` polish from its best point; preferable when event
        timing makes the objective non-smooth (see Notes).
    n_starts : int
        Number of starts (``'multistart'``): start 0 is ``initial`` if
        given, else :attr:`CalibrationProblem.x0` (each parameter's
        ``initial`` value; by default the value the experiment's run starts
        from, i.e. what the model and conditions already hold, else the
        midpoint of its bounds in the fitting space), the others are
        quasi-random points of ``sampler`` over the bounds. Every start is
        clipped into the bounds (a start on a bound is kept). Ignored by
        ``'hybrid'``.
    sampler : {'sobol', 'lhs'}
        Scrambled Sobol' or Latin-hypercube start points. Ignored by
        ``'hybrid'``.
    seed : None, int or numpy.random.SeedSequence
        Seed of the start points and of the global stage; the same seed
        gives the same fit. ``None`` draws fresh entropy (recorded in
        :attr:`FitResult.seed`).
    workers : int
        Processes over which the starts run (``'multistart'`` only; each
        builds its own model once). ``1`` (default) runs serially in this
        process. Under the ``spawn`` start method (macOS, Windows) every
        observable ``transform`` must be picklable (a module-level
        function, not a lambda), and a script must call this under
        ``if __name__ == '__main__':``.
    max_nfev : int, optional
        Maximum residual evaluations per ``least_squares`` run (``None``:
        SciPy's default, ``100 * n_free``).
    global_options : mapping, optional
        Options of :func:`~nskinetics.utils.hybrid_global_optimize` for
        ``'hybrid'``, over the defaults ``HYBRID_DEFAULTS``
        (``{'de_popsize': 10, 'de_maxiter': 30, 'do_basinhop': False,
        'polish_top_k': 1, 'local_maxiter': 50}``), e.g.
        ``{'de_maxiter': 100}``. Only its keyword options other than
        ``parallel``, ``workers`` and ``random_state`` are accepted (not
        ``scipy.optimize.minimize`` arguments; the ``least_squares`` polish
        follows anyway); any other key, or a count below 1, raises before
        any run.
    initial : FitResult, mapping or array_like, optional
        Start 0 of ``'multistart'`` in place of
        :attr:`CalibrationProblem.x0`, e.g. a previous fit's estimates to
        refit with another loss or more data: a :class:`FitResult` (its
        natural-unit ``estimates`` by name), ``{name: natural value}``, or
        an array in the fitting space (as for
        :meth:`CalibrationProblem.residuals`). It must name every free
        parameter of ``problem`` exactly. The starts table records it
        (clipped into the bounds) in ``'initial:<name>'``. Not used by
        ``'hybrid'`` (``ValueError``).

    Returns
    -------
    FitResult
        The best start (lowest cost among the starts whose final point
        evaluates).

    Raises
    ------
    TypeError, ValueError
        For invalid arguments, before any run.
    FittingError
        If every start fails (its start or final point cannot be
        evaluated).

    Warns
    -----
    UserWarning
        If an array ``initial`` lies outside the bounds of the fitting
        space (natural values must be given as a mapping).

    Notes
    -----
    A failed evaluation inside the optimizer (an integration failure or
    non-finite model values) never raises: every residual becomes a
    constant penalty (at least ``PENALTY_RESIDUAL``, ``1e3``, and ten
    times the largest residual of the start point), the trust region
    rejects the step, and the failure is counted in
    :attr:`FitResult.n_failed_evals`.

    When the parameters move events (e.g. feed spikes of a measured
    species) across sample times, the objective is only piecewise smooth,
    and local starts can stop on a wrong piece while reporting success.
    Prefer ``method='hybrid'``, more starts, or a robust-loss fit followed
    by a linear-loss refit started from it (``fit_parameters(linear_problem,
    initial=robust_fit)``, with ``linear_problem`` the same problem built
    with ``loss='linear'``, the default), and check the ``lack_of_fit``
    flag and how many starts reach the best cost; see
    :ref:`calibration-event-timing`.

    Examples
    --------
    >>> data = Dataset.from_long_csv(                     # doctest: +SKIP
    ...     'hplc.csv', experiment='run', time='t', observable='analyte',
    ...     value='g_L')
    >>> obs = [Observable('S', '[S]', sigma=0.5)]          # doctest: +SKIP
    >>> exps = [Experiment('FB1', cond1, data['FB1']),     # doctest: +SKIP
    ...         Experiment('FB2', cond2, data['FB2'])]
    >>> params = [FitParameter('k', (0.01, 10)),           # doctest: +SKIP
    ...           FitParameter('[S]', (20, 80), scope='local')]
    >>> problem = CalibrationProblem(km, exps, obs, params)  # doctest: +SKIP
    >>> fit = fit_parameters(problem, n_starts=10, seed=0)  # doctest: +SKIP
    >>> print(fit.summary())                               # doctest: +SKIP
    >>> profiles = profile_likelihood(fit, ['k'])          # doctest: +SKIP
    >>> bands = fit.predict('FB1', seed=1)                 # doctest: +SKIP
    """
    n_starts, workers, max_nfev, x_init = _check_fit_args(
        problem, method, n_starts, sampler, workers, max_nfev,
        global_options, initial)
    ss = _seed_sequence(seed)
    spec = problem._spec
    lo, hi = spec.lo, spec.hi
    options = {'max_nfev': max_nfev}
    private = problem._ready()
    global_message = None
    if method == 'multistart':
        starts = [problem.x0 if x_init is None
                  else _clip(x_init, lo, hi)]
        if n_starts > 1:
            unit = _unit_sample(n_starts - 1, problem.n_free, method=sampler,
                                seed=_child_seed(ss, 0))
            starts += [_clip(lo + u * (hi - lo), lo, hi) for u in unit]
        if workers == 1:
            records = [_run_start(_Evaluator(private, spec), x0, options)
                       for x0 in starts]
        else:
            task = functools.partial(_start_task, spec=spec, options=options)
            records = _ex._map_with_model(
                private, task, starts, workers=workers, chunksize=1,
                spawn_check=functools.partial(_check_transforms_picklable,
                                              problem))
    else:
        evaluator = _Evaluator(private, spec)
        r0, ok0 = evaluator.evaluate(problem.x0)
        if ok0:
            with np.errstate(all='ignore'):
                evaluator.penalty = max(PENALTY_RESIDUAL,
                                        10.0 * float(np.max(np.abs(r0))))
        opts = {**HYBRID_DEFAULTS, **(global_options or {})}
        state = int(_child_seed(ss, 1).generate_state(1)[0])
        stage = hybrid_global_optimize(
            lambda x: spec.cost(evaluator.residuals(x)),
            list(zip(lo.tolist(), hi.tolist())), parallel='none',
            random_state=state, **opts)
        global_message = (f'Global stage ({stage.nfev} objective '
                          f'evaluations): {stage.message}')
        records = [_run_start(evaluator, _clip(stage.x, lo, hi), options)]
        # the global stage's evaluations count toward the cost of the fit
        records[0]['nfev'] += int(stage.nfev)
    ok = [i for i, rec in enumerate(records) if rec['ok']]
    if not ok:
        detail = '; '.join(f'start {i}: {rec["message"]}'
                           for i, rec in enumerate(records[:5]))
        raise FittingError(f'All {len(records)} start(s) failed ({detail}). '
                           'Check that the bounds exclude regions where the '
                           'model cannot be integrated.')
    best = min(ok, key=lambda i: (records[i]['cost'], i))
    return FitResult._from_records(problem, records, best, method=method,
                                   seed=_seed_record(ss),
                                   global_message=global_message)


# --- the result --------------------------------------------------------------
@dataclass(eq=False, repr=False)
class FitResult:
    """Outcome of :func:`fit_parameters`.

    Built from the fields below (what :meth:`to_csv` stores); every
    attribute listed under Attributes is derived from them, so
    :meth:`from_csv` restores the result exactly.

    Parameters
    ----------
    names : list of str
        Free-parameter names.
    targets : list of str
        What each sets: a model id or initial-value selection, or an
        observable scale/offset name.
    experiments : list
        The experiment of each local parameter (``None`` for a global one).
    kinds : list of str
        ``'model'`` or ``'observable'`` (a scale/offset).
    log : numpy.ndarray
        Whether each parameter is fitted in log space.
    lower, upper : numpy.ndarray
        Bounds in natural units (exactly the :class:`FitParameter`
        bounds).
    x : numpy.ndarray
        Estimates in the fitting space (``ln p`` for log parameters).
    cost : float
        Optimizer objective ``0.5 * sum(rho(r**2))`` at the optimum.
    chi2 : float
        Plain weighted sum of squared residuals ``sum(r**2)`` at the
        optimum, whatever the loss.
    n_data : int
        Number of data points.
    jacobian : numpy.ndarray
        ``(n_data, n_free)`` weighted-residual Jacobian ``dr/dx`` at the
        optimum, in the fitting space (central finite differences; ``NaN``
        in a column whose differences failed).
    active_mask : numpy.ndarray
        Per parameter, ``-1`` / ``1`` if the estimate is at its lower /
        upper bound (``least_squares``'s ``active_mask``), else ``0``.
    fd_flags : pandas.DataFrame
        The finite-difference flags of the Jacobian (index ``parameter``;
        columns ``value, log, step, n_refinements, event_count_changed,
        reliable, scheme, error`` as in
        :attr:`~nskinetics.analysis.LocalSensitivityResult.flags`).
    residuals : pandas.DataFrame
        One row per data point: ``experiment, observable, time, data,
        model, sigma, weighted_residual`` (observable units; time in model
        units).
    starts : pandas.DataFrame
        One row per start: ``start, status ('ok'|'failed'), success
        (optimizer converged), best, cost, chi2, nfev, njev,
        n_failed_evals, message``, then ``'initial:<name>'`` and
        ``'estimate:<name>'`` (natural units) per parameter. For
        ``method='hybrid'`` (one row, the polish), ``nfev`` includes the
        global stage's objective evaluations.
    n_failed_evals : int
        Failed model evaluations over all starts (penalized, not raised).
    success : bool
        Whether the optimizer of the best start reported convergence.
    message : str
        Its message.
    method : str
        ``'multistart'`` or ``'hybrid'``, as given to
        :func:`fit_parameters`.
    loss : str
        The problem's loss (``'linear'`` or a robust loss, see
        :class:`CalibrationProblem`).
    f_scale : float
        The robust loss's scale, in units of sigma.
    absolute_sigma : bool
        The problem's setting: whether the sigmas are the standard
        deviations of the data (``True``) or weights only, with the noise
        level estimated from the residuals (``False``).
    seed : int or None
        The seed that reproduces the starts.
    units : dict
        The model's units (e.g. ``{'time': 'h', 'conc': 'g/L'}``).
    rank_rtol : float
        Relative threshold of both the rank test and the integration-noise
        column filter (see ``covariance``); default ``1e-6``, as in
        :func:`~nskinetics.analysis.identifiability`; in ``(0, 1)``. The
        default tolerances keep that noise only just below it (measured up
        to about ``9.8e-7``; see :ref:`calibration-rank`). It is not a
        :class:`CalibrationProblem` option: to judge a fit at another
        threshold (e.g. above the integration noise of looser
        ``tolerances=``), use ``dataclasses.replace(fit, rank_rtol=...)``,
        which recomputes every derived attribute and flag from the stored
        Jacobian without refitting. A higher value also counts a fit with
        ``condition_number > 1 / rank_rtol`` as ``rank_deficient``
        (infinite SE; from ``1e-4`` on, ``ill_conditioned`` can no longer
        fire), so raise it only as far as the noise requires.
    problem : CalibrationProblem or None
        The problem (not stored; ``None`` after :meth:`from_csv`).

    Attributes
    ----------
    estimates : pandas.Series
        Natural-unit estimates.
    n_free, dof : int
        Free parameters and ``n_data - n_free``.
    reduced_chi2 : float
        ``chi2 / dof`` (``NaN`` for ``dof = 0``).
    aic, aicc, bic : float
        With ``n = n_data`` and ``K`` estimated quantities: for absolute
        sigmas ``K = n_free``, ``AIC = chi2 + 2K`` and ``BIC = chi2 + K ln
        n`` (``-2 ln L`` up to a constant that is the same for every model
        fitted to the same data and sigmas); with ``absolute_sigma=False``
        (weights-only sigmas, Gaussian noise of unknown level)
        ``n ln(chi2/n)`` replaces ``chi2`` and the estimated noise level is
        counted, ``K = n_free + 1`` (Hurvich and Tsai, 1989).
        ``AICc = AIC + 2K(K+1)/(n-K-1)`` (``NaN`` for ``n - K - 1 <= 0``).
        With a robust loss they use the plain ``chi2`` at the robust
        optimum (outliers included), not the robust cost, so compare models
        by them only between linear-loss fits.
    covariance : pandas.DataFrame
        Covariance of ``x`` (fitting space): ``D (D J^T J D)^+ D``, the SVD
        pseudo-inverse on the numerical range of the equilibrated
        Jacobian ``J D`` (rank test at ``rank_rtol``), multiplied by
        ``reduced_chi2`` when ``absolute_sigma=False``. ``D`` is the
        diagonal column equilibration ``D_j = 1 / ||J_j||`` (1 for a zero
        column), for log and linear parameters alike, a near-optimal
        diagonal scaling (van der Sluis, 1969): every column of ``J D``
        has unit norm, so the rank test and ``condition_number`` measure
        only how nearly collinear the columns are, independent of the
        units, origin and bounds of the parameters and of the noise level.
        How precisely each parameter is determined is judged separately
        (``rel_se``, the ``high_rel_se`` and ``wide_interval`` flags). At
        full rank ``D`` cancels: the covariance is ``(J^T J)^-1``. All
        ``NaN`` if ``J`` is not finite. Before equilibrating, a simulation
        column (``kinds`` ``'model'``) is set to 0 if the change of the
        weighted residuals over its nominal finite-difference step,
        ``||J_j|| h_j`` (fitting space; the step before any event-count
        refinement: ``rel_step`` in ``ln p``, ``rel_step * max(|p|, 1e-3
        (hi - lo))`` for a linear parameter; a refined column's derivative
        times that step), is at most ``rank_rtol`` times the largest such
        change of a simulation column whose difference is a derivative
        (reliable, or one-sided after an ``event_count_changed``
        refinement; not a central difference across an event jump): the
        outputs do not measurably depend on that parameter, and its column
        is integration noise, which equilibration would lift to unit norm
        and which, following the dynamics, would look correlated with the
        other columns. It is then
        an exact null direction (infinite ``se``, ``rank_deficient``,
        correlation 0 with the others) and the other parameters are
        unaffected, as in :func:`~nskinetics.analysis.identifiability`,
        whose rank test on the raw ``ln p`` columns rejects such a column
        too (for a ``log=False`` parameter the step is smallest near 0, so a
        weak one there measures a smaller change). Observable scales and
        offsets involve no integration and are never set to 0; nor do they
        set the reference (an offset's change over its step grows with its
        distance from 0). The comparison is among simulation columns only:
        if none of them is informative (for example, only observable scales
        or offsets are fitted besides an uninformative simulation
        parameter), the largest is its own reference and is kept. The
        margin over the integration noise is thin (see ``rank_rtol``), in
        particular on fed-batch and other event models, when only weakly
        influential parameters are fitted, or with a looser ``rtol`` or a
        smaller ``rel_step``: a parameter without effect can then be kept,
        show as ``high_rel_se`` or ``wide_interval`` and widen the other
        intervals; fix it and refit, tighten ``tolerances`` (e.g.
        ``(1e-14, 1e-12)``, at little extra cost), or re-judge with
        ``dataclasses.replace(fit, rank_rtol=...)`` (see
        :ref:`calibration-rank`).
    se : pandas.Series
        Standard errors in the fitting space (of ``ln p`` for log
        parameters); ``inf`` for a parameter on an unidentifiable direction.
    se_natural, rel_se : pandas.Series
        First-order natural-unit SE (``p * se`` for log parameters) and
        ``se_natural / |p|``.
    ci : pandas.DataFrame
        95% Wald interval in natural units (columns ``lower, upper``):
        ``x -+ q se`` mapped back, so asymmetric for log parameters, with
        ``q = 1.96`` for absolute sigmas and the Student-t quantile
        ``t(0.975, dof)`` with ``absolute_sigma=False`` (weights-only
        sigmas: the noise level is estimated from the residuals).
    correlation : pandas.DataFrame
        Correlation of the estimates: that of ``covariance``, except for
        parameters with an infinite ``se`` (on an unidentifiable
        direction), which ``covariance`` leaves out. Among those it is the
        limit of a vanishing ridge, the correlation of
        ``(J^T J + eps I)^-1`` as ``eps -> 0``, i.e. of the null-space
        projector: e.g. ``-1`` for ``k1`` and ``k2`` when the data determine
        only ``k1 * k2``. Between one of them and a parameter off that
        direction it is ``0``. (The ``correlation`` of
        :class:`~nskinetics.analysis.IdentifiabilityResult` is ``NaN``
        there instead.) The ridge limit is taken in the equilibrated
        coordinates of ``covariance`` (``D^-1 x``).
    numerical_rank : int or None
        Rank of the equilibrated Jacobian ``J D`` (see ``covariance``;
        ``None`` if ``J`` is not finite).
    condition_number : float
        ``s_max / s_min`` of the equilibrated Jacobian ``J D`` (that of
        ``D J^T J D`` is its square): 1 for orthogonal columns, large for
        nearly collinear ones, whatever the units, origin, bounds and
        precision of the parameters; ``inf`` with a zero column (one the
        outputs do not measurably depend on, see ``covariance``).
    singular_values : numpy.ndarray
        Singular values of the equilibrated weighted Jacobian ``J D``, in
        descending order; ``NaN`` if ``J`` is not finite.
    flags : pandas.DataFrame
        Warnings, one row each: ``flag, parameters (tuple), value,
        message``. Flags, in this order: ``not_converged``,
        ``failed_evaluations`` (value: the count, ``n_failed_evals``),
        ``lack_of_fit`` (absolute sigmas only: ``chi2`` above the 99.9%
        quantile of chi-square with ``dof`` degrees of freedom, value: the
        reduced chi-square; with a robust loss ``2 * cost`` is tested
        instead, value: ``2 * cost / dof``), ``no_covariance`` (a
        non-finite Jacobian, or ``absolute_sigma=False`` with
        ``dof = 0``), ``rank_deficient`` (rank < n_free, i.e.
        ``condition_number > 1 / rank_rtol``, ``1e6`` by default; value:
        the numerical rank) or else ``ill_conditioned`` (condition number
        of the equilibrated ``D J^T J D``, ``condition_number**2``, > 1e8:
        nearly collinear columns; value: that condition number),
        ``high_rel_se`` (``rel_se > 1``, only for a parameter whose bounds
        exclude 0; value: ``rel_se``), ``wide_interval`` (for a
        zero-capable ``log=False`` parameter, such as an offset or an
        initial value with lower bound 0, whose ``rel_se`` depends on
        where the origin is: a finite ``se`` whose 95% Wald interval is
        wider than the bounds, ``2 q se > upper - lower``, so the data do
        not determine it within the allowed range; value: ``2 q se /
        (upper - lower)``), ``high_correlation`` (pairs
        with ``|corr| > 0.95``; value: the signed correlation),
        ``at_bound`` (value: the estimate, natural units),
        ``event_count_changed`` and ``unreliable_derivative`` (from
        ``fd_flags``). ``value`` is ``NaN`` for ``not_converged``,
        ``no_covariance``, ``event_count_changed`` and
        ``unreliable_derivative``.

    Raises
    ------
    ValueError
        If the fields have inconsistent shapes, or ``rank_rtol`` is not
        finite or not in ``(0, 1)``.
    TypeError
        If ``rank_rtol`` is not a real number (a bool included).

    Notes
    -----
    With a robust loss the covariance is still the Gauss-Newton
    approximation (``covariance``) of the *plain* weighted least-squares
    objective at the optimum (the loss only changed where the optimum
    is), so it ignores the down-weighting of outliers. ``lack_of_fit`` is
    then tested on ``2 * cost = sum(rho(r**2))`` rather than on ``chi2``
    (which includes the outliers on purpose): ``rho(z) <= z`` for every
    robust loss, so the test is conservative, although a few gross
    outliers can still trigger it. Wald
    intervals assume the model is close to linear in ``x`` over the
    interval, and the covariance assumes a good fit (with absolute sigmas,
    ``chi2`` close to ``dof``; see ``lack_of_fit``); check flagged
    parameters with a profile likelihood.
    """
    names: list
    targets: list
    experiments: list
    kinds: list
    log: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    x: np.ndarray
    cost: float
    chi2: float
    n_data: int
    jacobian: np.ndarray
    active_mask: np.ndarray
    fd_flags: pd.DataFrame
    residuals: pd.DataFrame
    starts: pd.DataFrame
    n_failed_evals: int
    success: bool
    message: str
    method: str = 'multistart'
    loss: str = 'linear'
    f_scale: float = 1.0
    absolute_sigma: bool = True
    seed: int | None = None
    units: dict = field(default_factory=dict)
    rank_rtol: float = RANK_RTOL
    problem: object = field(default=None, compare=False)

    def __post_init__(self):
        self.names = [str(n) for n in self.names]
        n = len(self.names)
        self.targets = list(self.targets)
        self.experiments = list(self.experiments)
        self.kinds = list(self.kinds)
        self.log = np.ascontiguousarray(self.log, dtype=bool)
        self.lower = np.ascontiguousarray(self.lower, dtype=float)
        self.upper = np.ascontiguousarray(self.upper, dtype=float)
        self.x = np.ascontiguousarray(self.x, dtype=float)
        self.jacobian = np.ascontiguousarray(self.jacobian, dtype=float)
        self.active_mask = np.ascontiguousarray(self.active_mask, dtype=int)
        self.cost, self.chi2 = float(self.cost), float(self.chi2)
        self.n_data, self.n_failed_evals = int(self.n_data), \
            int(self.n_failed_evals)
        self.success = bool(self.success)
        self.f_scale = float(self.f_scale)
        self.rank_rtol = _check_rank_rtol(self.rank_rtol)
        self.absolute_sigma = bool(self.absolute_sigma)
        if self.jacobian.shape != (self.n_data, n):
            raise ValueError(f'jacobian must have shape ({self.n_data}, {n}); '
                             f'got {self.jacobian.shape}.')
        for arg in ('targets', 'experiments', 'kinds'):
            if len(getattr(self, arg)) != n:
                raise ValueError(f'{arg} must have one entry per parameter.')
        for arg in ('log', 'lower', 'upper', 'x', 'active_mask'):
            if getattr(self, arg).shape != (n,):
                raise ValueError(f'{arg} must have shape ({n},).')
        self._derive()

    # --- derived quantities --------------------------------------------------
    def _derive(self):
        names, n, k = self.names, self.n_data, len(self.names)
        index = pd.Index(names, name='parameter', dtype=object)
        with np.errstate(all='ignore'):
            natural = np.where(self.log, np.exp(self.x), self.x)
        # exp(ln(bound)) may miss the bound in the last bit: an estimate at
        # a bound equals it exactly
        natural = np.clip(natural, self.lower, self.upper)
        self.estimates = pd.Series(natural, index=index, name='estimate')
        self.n_free = k
        self.dof = n - k
        with np.errstate(all='ignore'):
            self.reduced_chi2 = (self.chi2 / self.dof if self.dof > 0
                                 else math.nan)
            base = (self.chi2 if self.absolute_sigma
                    else n * math.log(self.chi2 / n) if self.chi2 > 0
                    else -math.inf)
        # absolute_sigma=False: the noise level is estimated too (chi2 / n)
        n_est = k + (0 if self.absolute_sigma else 1)
        self.aic = base + 2 * n_est
        self.bic = base + n_est * math.log(n)
        self.aicc = (self.aic + 2 * n_est * (n_est + 1) / (n - n_est - 1)
                     if n - n_est - 1 > 0 else math.nan)
        scale = 1.0 if self.absolute_sigma else self.reduced_chi2
        space = ['log' if lg else 'linear' for lg in self.log]
        # The rank test and the pseudo-inverse work on the equilibrated
        # J D, D_j = 1/||J_j|| for every column (log and linear alike; 1
        # for a zero column, which is rank deficient anyway): column
        # equilibration, the near-optimal diagonal scaling (van der Sluis).
        # A relative singular-value threshold on the raw J would depend on
        # the units of linear parameters (dr/dp scales with 1/unit of p),
        # on where their origin is, and on how precisely each parameter is
        # determined; on J D the rank and the condition number measure only
        # how nearly collinear the columns are. Precision is judged
        # separately (rel_se, wide_interval). The results are mapped back
        # to x; at full rank D cancels. Equilibration would also lift a
        # column of pure integration noise (a parameter the outputs do not
        # depend on) to unit norm, where its noise, which follows the
        # dynamics, looks like a genuine direction correlated with the
        # others; such columns are zeroed first (_null_columns).
        sc = np.ones(k)
        if np.all(np.isfinite(self.jacobian)):
            self._unseen = self._null_columns()
            jac = np.where(self._unseen, 0.0, self.jacobian)
            with np.errstate(all='ignore'):
                inv_norm = 1.0 / np.linalg.norm(jac, axis=0)
            sc = np.where(np.isfinite(inv_norm) & (inv_norm > 0), inv_norm,
                          1.0)
            svd = _svd_analysis(jac * sc, rank_rtol=self.rank_rtol)
            with np.errstate(all='ignore'):
                # scale is NaN without degrees of freedom
                # (absolute_sigma=False)
                cov = svd.covariance * np.outer(sc, sc) * scale
                se = svd.se * sc * math.sqrt(scale)
            self.numerical_rank = svd.numerical_rank
            self.condition_number = svd.condition_number
            self.singular_values = svd.singular_values
            null_loading = svd.null_loading
            # the singular values are sorted, so the kept ones are a prefix
            # and the rows of vt beyond the rank span the null space
            basis = svd.vt[svd.numerical_rank:].T
            projector = basis @ basis.T
        else:
            self._unseen = np.zeros(k, dtype=bool)
            cov = np.full((k, k), math.nan)
            se = np.full(k, math.nan)
            self.numerical_rank = None
            self.condition_number = math.nan
            self.singular_values = np.full(k, math.nan)
            null_loading = np.full(k, math.nan)
            projector = np.zeros((k, k))
        self.covariance = pd.DataFrame(cov, index=index, columns=list(names))
        self.se = pd.Series(se, index=index, name='se')
        se_nat, rel = _rel_se(natural, space, se)
        self.se_natural = pd.Series(se_nat, index=index, name='se_natural')
        self.rel_se = pd.Series(rel, index=index, name='rel_se')
        q = _wald_quantile(self.absolute_sigma, self.dof)
        with np.errstate(all='ignore'):
            lo_t, hi_t = self.x - q * se, self.x + q * se
            lo_n = np.where(self.log, np.exp(lo_t), lo_t)
            hi_n = np.where(self.log, np.exp(hi_t), hi_t)
            d = np.sqrt(np.diag(cov))
            corr = cov / np.outer(d, d)
        corr = np.where(np.isfinite(corr), corr, math.nan)
        np.fill_diagonal(corr, np.where(d > 0, 1.0, math.nan))
        # Parameters on an unidentifiable direction (infinite se): the
        # pseudo-inverse drops that direction, so its correlation would be
        # that of the projection onto the identifiable subspace (e.g. +1 for
        # k1, k2 when only k1*k2 is determined). Use the limit of a vanishing
        # ridge, corr((D J^T J D + eps I)^-1) as eps -> 0, instead: the
        # correlation of the null-space projector among them (-1 there), and
        # 0 against parameters off that direction. (The diagonal D does not
        # change the correlation of a covariance; the ridge limit is that of
        # the equilibrated coordinates.)
        null = np.isinf(se)
        if null.any():
            with np.errstate(all='ignore'):
                corr_null = np.clip(
                    projector / np.outer(null_loading, null_loading), -1, 1)
            both = np.outer(null, null)
            corr = np.where(both, corr_null, corr)
            corr = np.where(np.logical_xor.outer(null, null), 0.0, corr)
            corr[np.diag_indices(k)] = np.where(null, 1.0, np.diag(corr))
        self.ci = pd.DataFrame({'lower': lo_n, 'upper': hi_n}, index=index)
        self.correlation = pd.DataFrame(corr, index=index,
                                        columns=list(names))
        self._null_loading = null_loading
        # rel_se is scale-free only for a parameter whose range excludes 0
        # (log, or linear with bounds of one sign); for a zero-capable
        # linear one (an offset, an initial value with lower bound 0) it
        # depends on where the origin is, so its precision is judged
        # against the only scale the user gave, the width of its bounds:
        # the 95% Wald interval must fit in it (the Wald counterpart of a
        # profile that crosses the threshold within the bounds)
        self._relative_ok = self.log | (self.lower > 0) | (self.upper < 0)
        with np.errstate(all='ignore'):
            self._width_ratio = 2 * q * se / (self.upper - self.lower)
        self._identifiable = np.where(
            self._relative_ok, _identifiable(natural, space, se),
            np.isfinite(se) & (self._width_ratio <= 1))
        self.flags = self._flag_table()

    def _null_columns(self):
        """Boolean mask of the simulation (``kinds == 'model'``) columns that
        the outputs do not measurably depend on.

        ``d_j = ||J_j|| h_j``, with ``h_j`` the nominal step in the
        fitting space (the step before any event-count refinement:
        ``rel_step`` for a log parameter, ``rel_step * max(|p|,
        LINEAR_STEP_FLOOR * (hi - lo))`` for a linear one; ``rel_step`` is
        common to every column and cancels, so it is left out), is the
        change of the weighted residuals (in units of sigma) over that
        step, so it does not depend on the parameter's unit (a linear
        parameter's step is smallest near 0, where a weak one measures a
        smaller change). Every simulation column carries integration noise
        of the same size (about ``_fd.FD_NOISE * rtol * ||y / sigma||``,
        whatever the parameter), so a column whose ``d_j`` is at most
        ``rank_rtol`` times the largest one of a simulation column whose
        difference is a derivative (reliable, or one-sided after an
        ``event_count_changed`` refinement; not a central difference
        across an event jump) is taken as that noise: a null direction.
        The reference is about ``s_max * ||y / sigma|| * h``, with
        ``s_max = ||d(y/sigma)/d ln p|| / ||y/sigma||`` the normalized
        sensitivity of the most influential fitted simulation parameter,
        so the noise column sits at about ``FD_NOISE * rtol / (s_max h)``
        of it. With equal relative steps this is the rank criterion of
        :func:`~nskinetics.analysis.identifiability` on the raw ``d r / d
        ln p`` columns, whose ``3 rtol / h`` estimate assumes ``s_max``
        of about 1. The nominal step, not ``fd_flags['step']``, is used
        because an event-count refinement divides a column's step (by 10
        each time) but not the integration noise of the others; a
        one-sided column's step is the largest at which its side agreed,
        not the nominal one divided by ``10**n_refinements``. Observable
        scales and offsets involve no integration and are never null here;
        nor do they set the reference, since an offset's change over its
        step grows with its distance from 0. A column without a known step
        (no ``fd_flags`` row) is left as it is.
        """
        k = len(self.names)
        model = np.array([kind == 'model' for kind in self.kinds],
                         dtype=bool)
        flags = self.fd_flags.reindex(self.names)
        step = flags['step'].to_numpy(dtype=float)
        reliable = self.fd_flags['reliable'].reindex(
            self.names, fill_value=False).to_numpy(dtype=bool)
        one_sided = flags['scheme'].isin(('forward', 'backward')).to_numpy(
            dtype=bool)
        # over the nominal step without its common factor rel_step (not
        # fd_flags' step, which an event-count refinement divides)
        scale = np.where(self.log, 1.0,
                         _linear_scale(self.x, self.lower, self.upper))
        with np.errstate(all='ignore'):
            d = np.linalg.norm(self.jacobian, axis=0) * scale
        known = model & np.isfinite(step) & np.isfinite(d)
        # A central difference across an event jump (event_count_changed,
        # scheme 'central') is inflated by the jump; it must not set the
        # reference unless nothing else can. A one-sided difference kept at
        # the largest step whose firings agree with the nominal run is a
        # proper derivative (no jump in it) and may: when a sample sits at
        # a spike, every parameter that moves the spike is flagged, and
        # the reliable columns left may all be noise.
        reference = known & (reliable | one_sided)
        if not reference.any():
            reference = known
        if not reference.any():
            return np.zeros(k, dtype=bool)
        return known & (d <= self.rank_rtol * float(np.max(d[reference])))

    def _flag_table(self):
        names = self.names
        rows = []

        def add(flag, params, value, message):
            rows.append({'flag': flag, 'parameters': tuple(params),
                         'value': float(value), 'message': message})

        if not self.success:
            add('not_converged', (), math.nan,
                f'The optimizer did not report convergence: {self.message}')
        if self.n_failed_evals:
            add('failed_evaluations', (), self.n_failed_evals,
                f'{self.n_failed_evals} model evaluation(s) failed during the '
                'fit and were given a penalty residual; check that the bounds '
                'exclude regions where the model cannot be integrated.')
        # A chi-square goodness-of-fit test. A robust loss tests 2 * cost =
        # sum(rho(r**2)) instead of chi2, which deliberately includes the
        # down-weighted outliers: rho(z) <= z for every robust loss, so
        # 2 * cost <= chi2 and the test is conservative (a false alarm rate
        # of at most 1 - LACK_OF_FIT_LEVEL on data without outliers).
        robust = self.loss != 'linear'
        stat = 2.0 * self.cost if robust else self.chi2
        if self.absolute_sigma and self.dof > 0 and np.isfinite(stat):
            limit = float(chi2_distribution.ppf(LACK_OF_FIT_LEVEL, self.dof))
            if robust and stat > limit:
                add('lack_of_fit', (), stat / self.dof,
                    f'2*cost = {stat:.4g} (the robust objective, never above '
                    f'chi2 = {self.chi2:.4g}) on {self.dof} degrees of '
                    f'freedom exceeds the {LACK_OF_FIT_LEVEL:.1%} chi-square '
                    f'quantile {limit:.4g}: even with outliers down-weighted, '
                    'the model misses the data by more than the given '
                    'standard deviations, and the Wald uncertainties assume '
                    'it does not. Causes: a wrong local minimum or piece '
                    '(compare the costs in the starts table and look at '
                    "plot(); try more starts or method='hybrid'), gross "
                    "outliers (check the residuals; refit with loss='linear' "
                    'once they are dealt with), a model that cannot describe '
                    'the data, or standard deviations smaller than the real '
                    'noise.')
            elif not robust and stat > limit:
                add('lack_of_fit', (), self.reduced_chi2,
                    f'chi2 = {self.chi2:.4g} on {self.dof} degrees of freedom '
                    f'(reduced chi2 {self.reduced_chi2:.3g}) exceeds the '
                    f'{LACK_OF_FIT_LEVEL:.1%} chi-square quantile '
                    f'{limit:.4g}: the model misses the data by more than the '
                    'given standard deviations, and the Wald uncertainties '
                    'assume it does not. Causes: a wrong local minimum '
                    '(compare the costs in the starts table and look at '
                    "plot(); try more starts or method='hybrid'), a model "
                    'that cannot describe the data, or standard deviations '
                    'smaller than the real noise or outliers (if the sigmas '
                    'are only weights, use absolute_sigma=False).')
        if self.numerical_rank is None:
            add('no_covariance', (), math.nan,
                'The residual Jacobian at the optimum is not finite (a finite '
                'difference failed), so no covariance or standard errors are '
                'available.')
        elif not self.absolute_sigma and not self.dof > 0:
            add('no_covariance', (), math.nan,
                'absolute_sigma=False scales the covariance by the reduced '
                'chi-square, which needs more data points than free '
                f'parameters (dof = {self.dof}), so no covariance, standard '
                'errors, intervals or samples are available; give absolute '
                'sigmas or add data.')
        if self.numerical_rank is not None:
            with np.errstate(all='ignore'):
                fim_condition = float(np.float64(self.condition_number) ** 2)
            if self.numerical_rank < self.n_free:
                null = [n for n, v in zip(names, self._null_loading)
                        if v > NULL_LOADING]
                unseen = [n for n, u in zip(names, self._unseen) if u]
                unseen = (f' The outputs do not measurably depend on '
                          f'{unseen} (a finite difference within the '
                          'integration noise).' if unseen else '')
                add('rank_deficient', null, self.numerical_rank,
                    f'The residual Jacobian has rank {self.numerical_rank} < '
                    f'{self.n_free} free parameters: {null} move along a '
                    'direction the data do not determine (infinite SE).'
                    f'{unseen} Fix or reparameterize them, or choose an '
                    'identifiable subset with '
                    'nskinetics.analysis.identifiability first.')
            elif fim_condition > CONDITION_MAX:
                add('ill_conditioned', (), fim_condition,
                    'The equilibrated Fisher information D J^T J D (D_j = '
                    '1/||J_j||) is ill-conditioned '
                    f'(condition number {fim_condition:.3g} > '
                    f'{CONDITION_MAX:g}; cond(J D) = '
                    f'{self.condition_number:.3g}): the columns of the '
                    'Jacobian are nearly collinear, so the estimates are '
                    'strongly coupled, and the Wald uncertainty of the weakly '
                    'determined combination rests on a linearization over a '
                    'wide range; check it with a profile likelihood.')
        for name, rel, se, ok in zip(names, self.rel_se, self.se,
                                     self._relative_ok):
            if ok and np.isfinite(se) and np.isfinite(rel) \
                    and rel > REL_SE_MAX:
                add('high_rel_se', (name,), rel,
                    f'{name}: relative standard error {rel:.3g} > '
                    f'{REL_SE_MAX:g}; it is practically not identifiable '
                    'from these data.')
        # zero-capable linear parameters (no relative scale); an infinite
        # se is already reported by rank_deficient
        q = _wald_quantile(self.absolute_sigma, self.dof)
        for j, name in enumerate(names):
            ratio, se = self._width_ratio[j], self.se.iloc[j]
            if not self._relative_ok[j] and np.isfinite(se) and ratio > 1:
                add('wide_interval', (name,), ratio,
                    f'{name}: its 95% Wald interval (width {2 * q * se:.3g}) '
                    f'is wider than its bounds ({self.lower[j]:g}, '
                    f'{self.upper[j]:g}); the data do not determine it '
                    'within the allowed range, and its weak column can '
                    'widen the other intervals (check it with a profile '
                    'likelihood, or fix it and refit).')
        corr = self.correlation.to_numpy()
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                c = corr[i, j]
                if np.isfinite(c) and abs(c) > CORRELATION_MAX:
                    add('high_correlation', (names[i], names[j]), c,
                        f'{names[i]} and {names[j]} are strongly correlated '
                        f'(corr = {c:+.4f}); the data determine mainly a '
                        'combination of them.')
        for j, name in enumerate(names):
            side = self.active_mask[j]
            if side:
                bound = self.lower[j] if side < 0 else self.upper[j]
                add('at_bound', (name,), self.estimates.iloc[j],
                    f'{name} = {self.estimates.iloc[j]:.6g} is at its '
                    f'{"lower" if side < 0 else "upper"} bound {bound:.6g}; '
                    'the optimum may lie beyond it and the Wald interval does '
                    'not apply.')
        flags = self.fd_flags
        for name in names:
            if name not in flags.index:
                continue
            row = flags.loc[name]
            if bool(row['event_count_changed']):
                if row['scheme'] in ('forward', 'backward'):
                    where, detail = '', (
                        'its derivative is one-sided and its Wald '
                        'uncertainty approximate (check it with a profile '
                        'likelihood).')
                elif isinstance(row['error'], str):
                    # one side changed the firings and the other side's
                    # run failed: no difference at all (the column is NaN)
                    where, detail = ' on one side of the estimate', (
                        'the run on the other side failed '
                        f'({row["error"]}), so it has no finite difference '
                        'and no Wald uncertainty (use a profile '
                        'likelihood).')
                else:
                    # neither side kept the nominal firings: the central
                    # difference spans the jump
                    where, detail = ' on either side of the estimate', (
                        'its finite difference spans the jump (not a '
                        'derivative) and its Wald uncertainty is not '
                        'meaningful (use a profile likelihood).')
                add('event_count_changed', (name,), math.nan,
                    f'{name}: a perturbation{where} changes how many times '
                    '(or between which output times) an event fires, so the '
                    f'objective is not smooth there; {detail}')
            elif not bool(row['reliable']):
                add('unreliable_derivative', (name,), math.nan,
                    f'{name}: the finite difference is unreliable '
                    f'({row["error"]}).')
        columns = {c: [row[c] for row in rows] for c in _FLAG_SCHEMA}
        return _typed_frame(columns, _FLAG_SCHEMA)

    # --- construction --------------------------------------------------------
    @classmethod
    def _from_records(cls, problem, records, best, *, method, seed,
                      global_message=None):
        spec = problem._spec
        rec = records[best]
        names = list(spec.names)
        data = problem._data
        r, model = rec['r'], rec['model']
        residuals = _typed_frame({
            'experiment': data['experiment'], 'observable': data['observable'],
            'time': data['time'], 'data': data['data'], 'model': model,
            'sigma': data['sigma'], 'weighted_residual': r},
            _RESIDUAL_SCHEMA)
        columns = {c: [] for c in _START_SCHEMA}
        extra = {}
        for i, s in enumerate(records):
            columns['start'].append(i)
            columns['status'].append('ok' if s['ok'] else 'failed')
            columns['success'].append(bool(s['success']))
            columns['best'].append(i == best)
            for c in ('cost', 'chi2'):
                columns[c].append(float(s[c]))
            for c in ('nfev', 'njev'):
                columns[c].append(int(s[c]))
            columns['n_failed_evals'].append(int(s.get('n_failed', 0)))
            columns['message'].append(s['message'])
            x0 = spec.natural(s['x0'])
            x = spec.natural(s['x']) if s['ok'] else np.full(len(names),
                                                             math.nan)
            for j, name in enumerate(names):
                extra.setdefault(f'initial:{name}', []).append(float(x0[j]))
            for j, name in enumerate(names):
                extra.setdefault(f'estimate:{name}', []).append(float(x[j]))
        schema = {**_START_SCHEMA, **{c: float for c in extra}}
        starts = _typed_frame({**columns, **extra}, schema)
        message = rec['message']
        if global_message:
            message = f'{global_message}; polish: {message}'
        return cls(
            names=names, targets=list(problem._target_of),
            experiments=list(problem._experiment_of),
            kinds=list(problem._kinds), log=spec.log.copy(),
            lower=problem._natural_lo.copy(),
            upper=problem._natural_hi.copy(), x=rec['x'], cost=rec['cost'],
            chi2=rec['chi2'],
            n_data=spec.n_data, jacobian=rec['jacobian'],
            active_mask=rec['active_mask'], fd_flags=rec['fd_flags'],
            residuals=residuals, starts=starts,
            n_failed_evals=int(sum(s.get('n_failed', 0) for s in records)),
            success=rec['success'], message=message, method=method,
            loss=spec.loss, f_scale=spec.f_scale,
            absolute_sigma=problem.absolute_sigma, seed=seed,
            units=dict(problem._units), problem=problem)

    # --- views ---------------------------------------------------------------
    @property
    def warnings(self):
        """The flag messages (list of str), most serious first; the
        matching flag codes are ``flags['flag']`` (:meth:`summary` prints
        each as ``[flag] message``)."""
        return list(self.flags['message'])

    @property
    def table(self):
        """One row per parameter (``pandas.DataFrame`` indexed by name):
        ``target, experiment, estimate, se_natural, rel_se, ci_lower,
        ci_upper, lower, upper, log, se, identifiable, at_bound, reliable,
        event_count_changed``.

        ``ci_lower, ci_upper`` are the 95% Wald interval (the ``lower,
        upper`` columns of ``ci``), while ``lower, upper`` are the
        :class:`FitParameter` bounds (the result's ``lower`` and ``upper``
        fields). ``log`` says whether the parameter was fitted in log space
        and ``se`` is in that fitting space. ``identifiable`` is the
        practical verdict of :func:`~nskinetics.analysis.identifiability`:
        a finite ``se`` (no loading on a null direction) and ``rel_se <=
        REL_SE_MAX`` (1); the ``rel_se`` check applies only to parameters
        whose bounds exclude 0 (log parameters, and linear ones with bounds
        of one sign). For a zero-capable linear parameter, such as an
        offset, ``rel_se`` depends on where the origin is, so it is
        identifiable when its ``se`` is finite and its 95% Wald interval
        is no wider than its bounds (``2 q se <= upper - lower``, else the
        ``wide_interval`` flag), the Wald counterpart of a profile that
        crosses the threshold within the bounds. ``reliable`` and
        ``event_count_changed`` come from ``fd_flags``."""
        fd = self.fd_flags.reindex(self.names)
        return pd.DataFrame({
            'target': pd.Series(self.targets, dtype=object).to_numpy(),
            'experiment': pd.Series(self.experiments, dtype=object).to_numpy(),
            'estimate': self.estimates.to_numpy(),
            'se_natural': self.se_natural.to_numpy(),
            'rel_se': self.rel_se.to_numpy(),
            'ci_lower': self.ci['lower'].to_numpy(),
            'ci_upper': self.ci['upper'].to_numpy(),
            'lower': self.lower, 'upper': self.upper, 'log': self.log,
            'se': self.se.to_numpy(), 'identifiable': self._identifiable,
            'at_bound': self.active_mask != 0,
            'reliable': fd['reliable'].to_numpy(dtype=bool),
            'event_count_changed': fd['event_count_changed'].to_numpy(
                dtype=bool)},
            index=self.estimates.index)

    def summary(self):
        """A text report: warnings first, then the parameter table, then
        the fit statistics and the model's units.

        Each warning is a row of ``flags``, its message prefixed with the
        flag in brackets (e.g. ``[lack_of_fit]``). Text lines are wrapped
        at ``SUMMARY_WIDTH`` (88) characters, except the rows of the
        parameter table. For a multistart fit with more
        than one evaluated start, the statistic ``starts at the best cost``
        counts the evaluated starts whose cost is within ``2 (cost_i -
        cost) / s2 < 1`` of the best (``s2`` as in
        :func:`profile_likelihood`: 1 for absolute sigmas, else the reduced
        chi-square), i.e. ``delta chi2 < 1`` for the linear loss. Only one
        start at the best cost suggests that more starts (or
        ``method='hybrid'``) may find a lower one.

        Returns
        -------
        str
            The report (it displays as plain text in a notebook; ``print``
            it in a script).
        """
        def fill(text, first='', rest=''):
            return textwrap.fill(text, width=SUMMARY_WIDTH,
                                 initial_indent=first, subsequent_indent=rest,
                                 break_long_words=False,
                                 break_on_hyphens=False) or first.rstrip()

        ok = self.starts['status'] == 'ok'
        n_ok = int(ok.sum())
        done = ('global stage + least-squares polish'
                if self.method == 'hybrid' else
                f'{n_ok}/{len(self.starts)} start(s) evaluated')
        lines = [fill(f'Calibration fit ({self.method}, loss={self.loss!r})'),
                 fill(f'{self.n_free} free parameter(s), {self.n_data} data '
                      f'points, {done}', rest='  '), '']
        if len(self.flags):
            lines.append(f'Warnings ({len(self.flags)}):')
            lines += [fill(f'[{flag}] {msg}', '  - ', '    ')
                      for flag, msg in zip(self.flags['flag'],
                                           self.flags['message'])]
        else:
            lines.append('Warnings: none')
        lines += ['', 'Parameters (natural units; SE first order; 95% Wald '
                  'CI):']
        tbl = pd.DataFrame({
            'estimate': self.estimates, 'SE': self.se_natural,
            'rel SE': self.rel_se, 'CI lower': self.ci['lower'],
            'CI upper': self.ci['upper'], 'lower bound': self.lower,
            'upper bound': self.upper,
            'space': ['log' if lg else 'linear' for lg in self.log]},
            index=self.estimates.index)
        # a renamed copy: the index object is shared with the result's
        # frames, so setting its name would change theirs
        tbl = tbl.rename_axis(None)
        text = tbl.to_string(float_format=lambda v: f'{v:.6g}')
        lines += ['  ' + line for line in text.splitlines()]
        lines += ['', 'Fit statistics:']
        stats = [('chi2 (weighted SSR)', self.chi2),
                 ('cost (0.5 sum rho)', self.cost),
                 ('n_data', self.n_data), ('n_free', self.n_free),
                 ('dof', self.dof), ('reduced chi2', self.reduced_chi2),
                 ('AIC', self.aic), ('AICc', self.aicc), ('BIC', self.bic),
                 ('rank / cond(J D)', f'{self.numerical_rank} / '
                                      f'{self.condition_number:.3g}'),
                 ('sigma', 'absolute' if self.absolute_sigma else
                  'weights only (absolute_sigma=False; covariance scaled '
                  'by reduced chi2)')]
        if self.method == 'multistart' and n_ok > 1:
            s2 = (self.reduced_chi2 if not self.absolute_sigma
                  and np.isfinite(self.reduced_chi2)
                  and self.reduced_chi2 > 0 else 1.0)
            costs = self.starts.loc[ok, 'cost'].to_numpy(dtype=float)
            with np.errstate(all='ignore'):
                at_best = int(np.sum(2.0 * (costs - self.cost) / s2 < 1.0))
            stats.append(('starts at the best cost', f'{at_best} of {n_ok}'))
        stats += [('failed evaluations', self.n_failed_evals),
                  ('success', self.success), ('message', self.message)]
        width = max(len(k) for k, _ in stats)
        for key, value in stats:
            if isinstance(value, float):
                value = f'{value:.6g}'
            first = f'  {key:<{width}}  '
            lines.append(fill(str(value), first, ' ' * len(first)))
        if self.units:
            units = ', '.join(f'{k} [{v}]' for k, v in self.units.items())
            lines += ['', fill(f'Model units: {units}; data times are in '
                               'model time units.', rest='  ')]
        return _Text('\n'.join(lines))

    def __repr__(self):
        estimates = ', '.join(f'{n}: {v:.4g}'
                              for n, v in self.estimates.items())
        return (f'FitResult(method={self.method!r}, names={self.names}, '
                f'estimates={{{estimates}}}, '
                f'chi2={self.chi2:.6g}, n_data={self.n_data}, '
                f'success={self.success}, flags={list(self.flags["flag"])})')

    # --- applying ------------------------------------------------------------
    def apply(self, model):
        """Write the global model-parameter estimates into a model.

        Global parameters are written as the model's current values (as
        :meth:`~nskinetics.KineticModel.set_value` does), which persist
        through baseline resets and are carried by
        :meth:`~nskinetics.KineticModel.copy`. Initial-value targets (species,
        ``'[x]'`` selections, compartments and other state that a baseline
        reset restores) are written as ``init(...)`` values (``'[S]'`` ->
        ``'init([S])'``), so later runs start from them. A parameter that an
        SBML initial assignment reads is written as a plain value, so the
        assigned initial values follow it at the next reset, whereas the fit
        held them fixed (:class:`CalibrationProblem` warns about such
        targets).

        Parameters
        ----------
        model : KineticModel
            The model to update (e.g. the one that was fitted).

        Returns
        -------
        dict
            ``{name: estimate}`` of what was *not* applied: local
            (per-experiment) estimates and observable scales/offsets. Use
            them in the experiments' conditions.

        Raises
        ------
        TypeError
            If ``model`` is not a :class:`~nskinetics.KineticModel`.
        ValueError
            If a target is not a settable id of ``model`` (e.g. ``init()`` of
            a species in the form its SBML does not declare); raised before
            anything is written.
        """
        if not isinstance(model, KineticModel):
            raise TypeError(f'model must be a KineticModel; got '
                            f'{type(model).__name__}.')
        writes, rest = [], {}
        for j, name in enumerate(self.names):
            value = float(self.estimates.iloc[j])
            if self.kinds[j] != 'model' or self.experiments[j] is not None:
                rest[name] = value
                continue
            target = self.targets[j]
            if (not target.startswith('init(')
                    and _ex._auto_placement(model, target) == 'initial'):
                target = f'init({target})'
            writes.append((target, value))
        _check_ids(model, 'FitResult.apply', [t for t, _ in writes])
        model._write_values(writes)
        return rest

    # --- uncertainty propagation ---------------------------------------------
    def sample(self, n, seed=None):
        """Draw parameter sets from the Laplace (Gaussian) approximation.

        The draws are centred on the estimates in the fitting space (``ln
        p`` for log parameters), with :attr:`covariance`, clipped to the
        bounds there, then mapped to natural units (so log-normal or
        log-t for log parameters). With absolute sigmas they are
        multivariate normal. With ``absolute_sigma=False`` they are
        multivariate Student-t with ``dof`` degrees of freedom: the
        posterior of a linear model with a flat prior and a noise level
        estimated from the residuals, consistent with the Student-t
        :attr:`ci` (heavier tails than a normal at small ``dof``).

        Parameters
        ----------
        n : int
            Number of draws (>= 1).
        seed : None, int or numpy.random.SeedSequence
            The same seed gives the same draws; ``None`` draws fresh entropy.

        Returns
        -------
        pandas.DataFrame
            ``n`` rows (index ``sample``), one column per free parameter
            (natural units, within the bounds).

        Raises
        ------
        TypeError, ValueError
            If ``n`` is invalid or the covariance is not finite (see the
            ``no_covariance`` flag; with ``absolute_sigma=False`` also when
            there are no degrees of freedom).

        Warns
        -----
        UserWarning
            If some parameters have infinite SE (rank-deficient Jacobian):
            the pseudo-inverse covariance gives no spread along the
            directions the data do not determine, so their draws vary far
            less than their actual uncertainty.

        Notes
        -----
        This is a *local* approximation: it assumes the objective is
        quadratic in the fitting space around the optimum. It understates
        the uncertainty of a parameter whose profile is skewed, flat or cut
        by a bound (draws beyond a bound are moved onto it, not
        redistributed); check such parameters with
        :func:`profile_likelihood`.
        """
        n = _count('n', n, 1)
        cov = self.covariance.to_numpy()
        if not np.all(np.isfinite(cov)):
            raise ValueError(
                'The covariance is not finite (see FitResult.flags), so '
                'there is no Gaussian approximation to sample from.')
        loose = [name for name, se in zip(self.names, self.se)
                 if not np.isfinite(se)]
        if loose:
            warnings.warn(
                f'{loose} have infinite standard errors (the data do not '
                'determine them): the samples do not vary along the '
                'unidentifiable direction(s) and understate their '
                'uncertainty; use profile_likelihood.', UserWarning,
                stacklevel=2)
        with np.errstate(all='ignore'):
            lo = np.where(self.log, np.log(self.lower), self.lower)
            hi = np.where(self.log, np.log(self.upper), self.upper)
        ss = _seed_sequence(seed)
        if self.absolute_sigma:
            draws = _gaussian(self.x, cov, n, seed=ss, lo=lo, hi=hi)
        else:
            # an estimated noise level: multivariate Student-t with dof
            # degrees of freedom (dof > 0, else cov is NaN), scaled before
            # clipping
            z = _gaussian(np.zeros_like(self.x), cov, n, seed=ss)
            w = np.random.default_rng(_child_seed(ss, 1)).chisquare(self.dof,
                                                                    n)
            with np.errstate(all='ignore'):
                draws = np.clip(self.x + z / np.sqrt(w / self.dof)[:, None],
                                lo, hi)
        with np.errstate(all='ignore'):
            natural = np.where(self.log, np.exp(draws), draws)
        # a draw clipped onto a bound is that bound exactly (exp(ln(bound))
        # may differ from it in the last bit)
        natural = np.where(draws <= lo, self.lower,
                           np.where(draws >= hi, self.upper, natural))
        natural = np.clip(natural, self.lower, self.upper)
        return pd.DataFrame(natural, columns=list(self.names),
                            index=pd.RangeIndex(n, name='sample'))

    def _model_for(self, model):
        """The model to run: ``model``, or the problem's private copy."""
        if model is not None:
            if not isinstance(model, KineticModel):
                raise TypeError(f'model must be a KineticModel; got '
                                f'{type(model).__name__}.')
            return model
        if self.problem is None:
            raise ValueError('This FitResult has no problem (it was read from '
                             'a file); pass model=.')
        return self.problem._ready()

    def _prediction_base(self, condition):
        """``(experiment, base condition, skipped)`` of a prediction:
        ``experiment`` is the experiment named by ``condition`` (else
        ``None``); ``skipped`` is ``(matches, names)`` of the local
        estimates an experiment's own Condition object leaves out."""
        skipped = None
        if isinstance(condition, str):
            if self.problem is None:
                raise ValueError('An experiment name needs the fitted problem '
                                 '(this FitResult was read from a file); pass '
                                 'its Condition.')
            names = [e.name for e in self.problem.experiments]
            if condition not in names:
                raise ValueError(f'Unknown experiment {condition!r}; the '
                                 f'problem has {names}.')
            return (condition,
                    self.problem.experiments[names.index(condition)].condition,
                    skipped)
        if not isinstance(condition, Condition):
            raise TypeError(f'condition must be a Condition or an experiment '
                            f'name; got {type(condition).__name__}.')
        # an experiment's own Condition object (or an equal one): its local
        # estimates are not applied, which is easy to miss
        matches = ([] if self.problem is None else
                   [e.name for e in self.problem.experiments
                    if e.condition == condition])
        names = [self.names[j] for j in range(self.n_free)
                 if self.kinds[j] == 'model' and self.experiments[j] in matches]
        if names:
            skipped = (matches, names)
        return None, condition, skipped

    def _prediction_runs(self, model, base, experiment, skipped, n_samples,
                         seed):
        """The runs of a prediction: ``(conditions, use, seed sequence)``,
        the nominal run (labelled ``'nominal'``, at the estimates) first,
        then one per draw of :meth:`sample` (``'sample <i>'``), each with
        the estimates ``use`` (indices of the global, and ``experiment``'s
        local, model parameters) pinned into ``base``."""
        use = [j for j in range(self.n_free) if self.kinds[j] == 'model'
               and self.experiments[j] in (None, experiment)]
        if not use:
            raise ValueError(
                'No fitted model parameter applies to this condition (local '
                'estimates apply only to their experiment: pass its name).')
        if experiment is None and skipped:
            matches, names = skipped
            warnings.warn(
                f'This condition equals the condition of experiment(s) '
                f'{matches}; their local estimates {names} are not applied '
                "(the condition's own values are used). Pass an experiment's "
                f'name instead (e.g. {matches[0]!r}) to apply them.',
                UserWarning, stacklevel=3)
        ss = _seed_sequence(seed)
        draws = self.sample(n_samples, seed=ss)
        targets = [self.targets[j] for j in use]
        place = {t: _ex._placement_for(model, base, t) for t in targets}
        values = [self.estimates.to_numpy()[use]] + [
            row for row in draws.to_numpy()[:, use]]
        conditions = [
            _fd.pin(model, base, dict(zip(targets, map(float, v))),
                    place).replace(label=f'sample {i - 1}' if i
                                   else 'nominal')
            for i, v in enumerate(values)]
        # user -> predict / predict_metrics -> here -> the helper
        _warn_lack_of_fit(self, 'predicted bands and percentiles',
                          stacklevel=4)
        return conditions, use, ss

    def predict(self, condition, n_samples=200, q=(5, 50, 95), seed=None,
                workers=1, outputs=None, *, model=None):
        """Prediction bands of model outputs under parameter uncertainty.

        ``n_samples`` parameter sets are drawn with :meth:`sample` and run
        with :func:`~nskinetics.run_many`; the bands are percentiles over
        the successful runs at every output time. The nominal prediction
        (at the estimates) is run too.

        Parameters
        ----------
        condition : Condition or str
            The run to predict, on its own output grid. A
            :class:`~nskinetics.Condition` gets each draw's *global*
            estimates of model parameters and initial values; these replace
            the condition's own entries for the same targets. Local
            (per-experiment) estimates are not applied; the condition's own
            values are used for those targets. Pass an experiment's name to
            run its condition with the global *and* its local estimates
            (passing that experiment's Condition object instead leaves out
            its local estimates, with a warning; see Warns).
        n_samples : int
            Number of draws (default 200).
        q : sequence of float
            Percentiles in ``[0, 100]`` (default ``(5, 50, 95)``).
        seed : None, int or numpy.random.SeedSequence
            Seed of the draws (as for :meth:`sample`: the same seed gives
            the same parameter sets).
        workers : int
            Processes of :func:`~nskinetics.run_many` (default 1, serial).
        outputs : sequence of str, optional
            Model selections to predict (e.g. ``['[S]', '[P]']``), not
            observable names (see Notes). ``None``: the condition's
            selections, or the model's default ones (floating species
            concentrations).
        model : KineticModel, optional
            The model to run. Default: the fitted problem's private copy
            (with its integrator tolerances); required after
            :meth:`from_csv`.

        Returns
        -------
        pandas.DataFrame
            Tidy table, one row per output and time: ``time`` (model time
            units), ``output``, ``nominal``, then one column per percentile
            named ``'q<value>'`` (``q5, q50, q95``). ``attrs`` holds
            ``n_samples``, ``n_failed`` (failed sampled runs),
            ``quantiles``, ``parameters`` (the names applied), ``seed``
            (the int that reproduces the draws, if any) and ``units``.

        Raises
        ------
        TypeError, ValueError
            For invalid arguments, an unknown experiment name, a condition
            no fitted model parameter applies to, or a covariance that cannot
            be sampled (see :meth:`sample`); before any run.
        KineticSimulationError
            If every run (the nominal one and every sample) fails.

        Warns
        -----
        UserWarning
            If sampled runs (or the nominal run) fail; the percentiles use
            the successful runs (``NaN`` if none). If an output is the name
            of one of the fitted problem's observables that records another
            selection (e.g. the output ``'S'`` with the observable
            ``Observable('S', '[S]')``), before any run. If ``condition``
            is a Condition equal to the condition of an experiment that has
            local estimates of model parameters, which are then not applied
            (pass the experiment's name to apply them), before any run. If
            the fit is flagged ``lack_of_fit`` (:attr:`flags`), before any
            run: the bands assume that the model matches the data, are
            local to this optimum and may be far too narrow.

        Notes
        -----
        The bands inherit the Laplace approximation of :meth:`sample`
        (normal draws, or Student-t ones with ``absolute_sigma=False``) and
        only the parameter uncertainty (no measurement noise is added).
        For scalar metrics of a run (titer, yield, productivity, spike
        count), use :meth:`predict_metrics`: a percentile of a ratio is not
        the ratio of the percentiles.

        The bands are of model selections in model units: observables
        (their scale, offset and transform) are not applied, and fitted
        observable scales/offsets are not drawn (``attrs['parameters']``
        lists what was applied). ``outputs`` are therefore model
        selections: an observable's name is not resolved as the
        observable. If the name is also a model id, it is read as that
        model quantity (e.g. ``'S'`` is the amount of species ``S``, while
        an observable ``'S'`` may record the concentration ``'[S]'``; a
        ``UserWarning`` says so); otherwise the run validation raises
        ``ValueError``. For bands of a measured signal, predict the
        observable's selection and apply the observable to the predicted
        columns, drawing its scale/offset from :meth:`sample` if it was
        fitted.
        """
        model = self._model_for(model)
        n_samples = _count('n_samples', n_samples, 1)
        q = _quantiles(q)
        workers = _ex._check_workers(workers)
        experiment, base, skipped = self._prediction_base(condition)
        if outputs is not None:
            if isinstance(outputs, str):
                outputs = [outputs]
            outputs = list(outputs)
            bad = [o for o in outputs if not isinstance(o, str) or not o]
            if bad or not outputs:
                raise ValueError(f'outputs must be non-empty selection '
                                 f'strings; got {outputs!r}.')
            if len(set(outputs)) != len(outputs):
                raise ValueError(f'Duplicate outputs in {outputs}.')
            if self.problem is not None:
                _warn_observable_outputs(self.problem, outputs)
            base = base.replace(selections=tuple(outputs))
        sels = [s for s in model._run_selections(base) if s != 'time']
        if not sels:
            raise ValueError('No output to predict; pass outputs=.')
        conditions, use, ss = self._prediction_runs(model, base, experiment,
                                                    skipped, n_samples, seed)
        batch = _ex.run_many(model, conditions, workers=workers)
        ok = batch.ok
        n_failed = int((~ok[1:]).sum())
        if n_failed:
            warnings.warn(f'{n_failed} of {n_samples} sampled runs failed; '
                          'the percentiles use the others.', UserWarning,
                          stacklevel=2)
        if not ok[0]:
            warnings.warn('The nominal run failed: '
                          f'{batch.metrics["error"].iloc[0]}', UserWarning,
                          stacklevel=2)
        if not ok.any():
            raise KineticSimulationError(
                f'Every run of the prediction failed (nominal: '
                f'{batch.metrics["error"].iloc[0]}).')
        ref = next(r for r, good in zip(batch.results, ok) if good)
        grid, n_t = ref.time, len(ref)
        nominal = (batch.results[0].df[sels].to_numpy(dtype=float) if ok[0]
                   else np.full((n_t, len(sels)), math.nan))
        stack = [r.df[sels].to_numpy(dtype=float)
                 for r, good in zip(batch.results[1:], ok[1:]) if good]
        with np.errstate(all='ignore'):
            bands = (np.percentile(np.stack(stack), q, axis=0) if stack
                     else np.full((len(q), n_t, len(sels)), math.nan))
        data = {'time': np.tile(grid, len(sels)),
                'output': pd.Series(np.repeat(sels, n_t), dtype=object),
                'nominal': nominal.T.ravel()}
        for i, value in enumerate(q):
            data[_quantile_name(value)] = bands[i].T.ravel()
        out = pd.DataFrame(data)
        out.attrs.update(n_samples=n_samples, n_failed=n_failed,
                         quantiles=tuple(q),
                         parameters=[self.names[j] for j in use],
                         seed=_seed_record(ss), units=dict(self.units))
        return out

    def predict_metrics(self, condition, metrics, n_samples=200,
                        q=(5, 50, 95), seed=None, workers=1, *, model=None):
        """Percentiles of run metrics under parameter uncertainty.

        The counterpart of :meth:`predict` for scalar metrics of a run
        (titer, yield, productivity, spike count...): the same draws (the
        same ``seed`` gives the same parameter sets as :meth:`predict` and
        :meth:`sample`), run with :func:`~nskinetics.run_many`, each metric
        evaluated on every run, and its percentiles taken over the
        successful runs. A percentile of a ratio such as a yield is not the
        ratio of the percentiles of its parts, so derived metrics need this
        rather than the bands of :meth:`predict`.

        Parameters
        ----------
        condition : Condition or str
            The run to predict, as for :meth:`predict`: a
            :class:`~nskinetics.Condition` gets the *global* estimates; an
            experiment's name also gets its local estimates.
        metrics : mapping or sequence
            ``{name: callable}`` or a sequence of named metrics, as for
            :func:`~nskinetics.run_many` (e.g. ``[Final('[P]'),
            Yield('[P]', '[S]', fed=...), EventCount('n')]`` from
            :mod:`nskinetics.analysis`); at least one. Each is in its own
            units (model time and concentration units for the
            :mod:`nskinetics.analysis` metrics).
        n_samples : int
            Number of draws (default 200).
        q : sequence of float
            Percentiles in ``[0, 100]`` (default ``(5, 50, 95)``).
        seed : None, int or numpy.random.SeedSequence
            Seed of the draws (as for :meth:`sample`).
        workers : int
            Processes of :func:`~nskinetics.run_many` (default 1, serial).
        model : KineticModel, optional
            The model to run. Default: the fitted problem's private copy;
            required after :meth:`from_csv`.

        Returns
        -------
        pandas.DataFrame
            One row per metric (index ``metric``): ``nominal`` (at the
            estimates; ``NaN`` if that run failed), one column per
            percentile named ``'q<value>'`` (``q5, q50, q95``) over the
            successful sampled runs with a finite value of that metric
            (``NaN`` if none), and ``n_ok``, the number of such runs.
            ``attrs`` holds ``n_samples``, ``n_failed`` (failed sampled
            runs), ``quantiles``, ``parameters`` (the names applied),
            ``seed`` (the int that reproduces the draws, if any) and
            ``units``.

        Raises
        ------
        TypeError, ValueError
            For invalid arguments (including the metrics, as for
            :func:`~nskinetics.run_many`), an unknown experiment name, a
            condition no fitted model parameter applies to, or a covariance
            that cannot be sampled (see :meth:`sample`); before any run.
        KineticSimulationError
            If every run (the nominal one and every sample) fails.

        Warns
        -----
        UserWarning
            As for :meth:`predict`: if sampled runs (or the nominal run)
            fail, if ``condition`` is an experiment's Condition object
            whose local estimates are then not applied, and if the fit is
            flagged ``lack_of_fit`` (the percentiles may then be far too
            narrow).

        Notes
        -----
        The percentiles inherit the Laplace approximation of :meth:`sample`
        and only the parameter uncertainty. Fitted observable scales and
        offsets are not drawn: metrics read model selections.
        """
        model = self._model_for(model)
        n_samples = _count('n_samples', n_samples, 1)
        q = _quantiles(q)
        workers = _ex._check_workers(workers)
        names = list(_ex._normalize_metrics(metrics))
        if not names:
            raise ValueError('metrics must name at least one metric.')
        experiment, base, skipped = self._prediction_base(condition)
        conditions, use, ss = self._prediction_runs(model, base, experiment,
                                                    skipped, n_samples, seed)
        batch = _ex.run_many(model, conditions, metrics=metrics,
                             workers=workers, keep_results=False)
        ok = batch.ok
        n_failed = int((~ok[1:]).sum())
        if n_failed:
            warnings.warn(f'{n_failed} of {n_samples} sampled runs failed; '
                          'the percentiles use the others.', UserWarning,
                          stacklevel=2)
        if not ok[0]:
            warnings.warn('The nominal run failed: '
                          f'{batch.metrics["error"].iloc[0]}', UserWarning,
                          stacklevel=2)
        if not ok.any():
            raise KineticSimulationError(
                f'Every run of the prediction failed (nominal: '
                f'{batch.metrics["error"].iloc[0]}).')
        values = batch.metrics[names].to_numpy(dtype=float)
        nominal = values[0] if ok[0] else np.full(len(names), math.nan)
        samples = values[1:][ok[1:]]
        data = {'nominal': nominal}
        bands, n_ok = np.full((len(q), len(names)), math.nan), []
        for k in range(len(names)):
            finite = samples[:, k][np.isfinite(samples[:, k])]
            n_ok.append(finite.size)
            if finite.size:
                bands[:, k] = np.percentile(finite, q)
        for i, value in enumerate(q):
            data[_quantile_name(value)] = bands[i]
        data['n_ok'] = np.asarray(n_ok, dtype=np.int64)
        out = pd.DataFrame(data, index=pd.Index(names, name='metric',
                                                dtype=object))
        out.attrs.update(n_samples=n_samples, n_failed=n_failed,
                         quantiles=tuple(q),
                         parameters=[self.names[j] for j in use],
                         seed=_seed_record(ss), units=dict(self.units))
        return out

    # --- plotting ------------------------------------------------------------
    def _curves(self, n_points):
        """``{experiment: (time, {observable: values}, event times)}`` of the
        model at the estimates on a dense grid (``None`` without a problem
        or if a run fails)."""
        problem = self.problem
        if problem is None:
            return None
        spec, private = problem._spec, problem._ready()
        dense = []
        for c, e in zip(spec.conditions, problem.experiments):
            end = (e.condition.t_end if e.condition.times is None
                   else e.condition.times[-1])
            grid = np.union1d(np.linspace(c.times[0], end, n_points),
                              c.times)
            dense.append(c.replace(times=tuple(float(t) for t in grid)))
        point = _fd.evaluate(private, dense, spec.fd_parameters(self.x),
                             _no_values, event_counters=spec.counters)
        if not point.ok:
            warnings.warn(f'The model curves could not be computed '
                          f'({point.error}); plotting the model at the data '
                          'times.', UserWarning, stacklevel=3)
            return None
        out = {}
        for e, result, blocks in zip(problem.experiments, point.results,
                                     spec.blocks):
            curves = {}
            for b in blocks:
                scale, offset = (point.values[v] if isinstance(v, str) else v
                                 for v in (b.scale, b.offset))
                curves[b.observable.name] = b.observable._apply(
                    result[b.observable.selection], scale, offset)
            t = result.time
            events = set()
            for counter in spec.counters:
                with np.errstate(all='ignore'):
                    jump = np.diff(result[counter]) != 0
                events.update(t[1:][jump].tolist())
            out[e.name] = (t, curves, sorted(events))
        return out

    def plot(self, *, n_points=PLOT_POINTS):
        """Plot the data against the fitted model, per experiment and
        observable.

        Each panel shows the data with ``sigma`` error bars and the model at
        the estimates on a dense time grid (``n_points`` times from ``t0``
        to the experiment condition's ``t_end``, plus the data times), with
        dotted vertical lines where an event counter (e.g. a FeedSpike's
        ``count_var``) changes: at the first output time after the event,
        so within ``(t_end - t0) / n_points`` of it. Without the problem
        (after :meth:`from_csv`) the model is drawn at the data times only,
        from :attr:`residuals`.

        Parameters
        ----------
        n_points : int
            Output times of the dense model curve (>= 2; default 200).

        Returns
        -------
        (matplotlib.figure.Figure, numpy.ndarray)
            The figure and its axes, shape ``(n_experiments,
            n_observables)`` (a panel without data is switched off).

        Raises
        ------
        TypeError, ValueError
            For an invalid ``n_points``.

        Warns
        -----
        UserWarning
            If the dense model runs fail (the model is then drawn at the
            data times).
        """
        n_points = _count('n_points', n_points, 2)
        res = self.residuals
        exps = list(dict.fromkeys(res['experiment']))
        obs = list(dict.fromkeys(res['observable']))
        curves = self._curves(n_points)
        _set_font()
        fig, axes = plt.subplots(len(exps), len(obs), squeeze=False,
                                 figsize=(4.2 * len(obs), 3.2 * len(exps)))
        time_label = _bf('Time') + f' [{self.units.get("time", "")}]'
        legend = True
        for i, e in enumerate(exps):
            for k, o in enumerate(obs):
                ax = axes[i, k]
                rows = res[(res['experiment'] == e) & (res['observable'] == o)]
                if rows.empty:
                    ax.axis('off')
                    continue
                t = rows['time'].to_numpy(dtype=float)
                ax.errorbar(t, rows['data'].to_numpy(dtype=float),
                            yerr=rows['sigma'].to_numpy(dtype=float),
                            fmt='o', ms=4, color='black', ecolor='gray',
                            elinewidth=0.8, capsize=2, label='data')
                if curves is not None:
                    tt, ys, events = curves[e]
                    ax.plot(tt, ys[o], color='tab:blue', label='model')
                else:
                    order = np.argsort(t, kind='stable')
                    ax.plot(t[order], rows['model'].to_numpy()[order],
                            color='tab:blue', marker='x', ms=4,
                            label='model')
                    events = ()
                for m, te_ in enumerate(events):
                    ax.axvline(te_, color='dimgray', ls=':', lw=1,
                               label='event' if m == 0 else None)
                ax.set_title(f'{e}: {o}', fontsize=11)
                ax.set_xlabel(time_label)
                ax.set_ylabel(o)
                _style_axes(ax)
                if legend:
                    ax.legend(fontsize=8, frameon=False)
                    legend = False
        fig.tight_layout()
        return fig, axes

    # --- CSV -----------------------------------------------------------------
    def to_csv(self, path):
        """Write the result to a long-format CSV.

        Columns ``quantity, key, parameter, index, value``: ``meta`` rows
        (``key``: field), per-parameter rows ``x``, ``lower``, ``upper``,
        ``active_mask``, ``jacobian`` entries (``index``: data row),
        ``fd_flags`` (``key``: column), and the ``residuals`` and ``starts``
        tables (``key``: column, ``index``: row). Values are JSON, so
        :meth:`from_csv` restores every stored field exactly and recomputes
        the derived ones.

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
            If the units are not JSON-safe.
        """
        dumps = json.dumps
        meta = {'nskinetics': _CSV_FORMAT, 'names': self.names,
                'targets': self.targets, 'experiments': self.experiments,
                'kinds': self.kinds, 'log': [bool(v) for v in self.log],
                'cost': self.cost, 'chi2': self.chi2, 'n_data': self.n_data,
                'n_failed_evals': self.n_failed_evals,
                'success': self.success, 'message': self.message,
                'method': self.method, 'loss': self.loss,
                'f_scale': self.f_scale,
                'absolute_sigma': self.absolute_sigma, 'seed': self.seed,
                'units': self.units, 'rank_rtol': self.rank_rtol,
                'starts_columns': list(self.starts.columns),
                'n_starts': len(self.starts)}
        rows = [('meta', key, '', '', dumps(value))
                for key, value in meta.items()]
        for q in ('x', 'lower', 'upper', 'active_mask'):
            arr = getattr(self, q)
            rows += [(q, '', name, '', dumps(_py(arr[j])))
                     for j, name in enumerate(self.names)]
        for j, name in enumerate(self.names):
            rows += [('jacobian', '', name, i, dumps(float(v)))
                     for i, v in enumerate(self.jacobian[:, j])]
        for name, row in self.fd_flags.iterrows():
            rows += [('fd_flags', c, name, '', dumps(_py(row[c])))
                     for c in self.fd_flags.columns]
        for q in ('residuals', 'starts'):
            frame = getattr(self, q)
            for c in frame.columns:
                rows += [(q, c, '', i, dumps(_py(v)))
                         for i, v in enumerate(frame[c].tolist())]
        pd.DataFrame(rows, columns=list(_CSV_COLUMNS)).to_csv(
            path, index=False, lineterminator='\n', encoding='utf-8')

    @classmethod
    def from_csv(cls, path):
        """Read a result written by :meth:`to_csv`.

        Parameters
        ----------
        path : str or path-like

        Returns
        -------
        FitResult
            Equal in every stored field, with the derived attributes
            recomputed (so they are equal too); ``problem`` is ``None``.

        Raises
        ------
        ValueError
            If the file was not written by :meth:`to_csv` or is incomplete
            (a missing meta key, or a missing or duplicated entry).
        """
        meta, groups = _read_long_csv(path, _CSV_FORMAT, required=(
            'names', 'targets', 'experiments', 'kinds', 'log', 'cost', 'chi2',
            'n_data', 'n_failed_evals', 'success', 'message', 'method',
            'loss', 'f_scale', 'absolute_sigma', 'seed', 'units', 'rank_rtol',
            'starts_columns', 'n_starts'))
        names, n_data = meta['names'], meta['n_data']
        k = len(names)
        par = {p: j for j, p in enumerate(names)}

        def scatter(quantity, shape, *axes, **kwargs):
            return _scatter(path, _CSV_FORMAT, groups, quantity, shape, axes,
                            **kwargs)

        x = scatter('x', (k,), ('parameter', par))
        lower = scatter('lower', (k,), ('parameter', par))
        upper = scatter('upper', (k,), ('parameter', par))
        active = scatter('active_mask', (k,), ('parameter', par), parse=int,
                         dtype=int)
        jac = scatter('jacobian', (n_data, k), ('index', None),
                      ('parameter', par))
        flag_cols = list(_fd.FLAG_COLUMNS)
        flag_rows = _read_table(path, groups, 'fd_flags', flag_cols,
                                names, 'parameter')
        fd_flags = _fd.flags_frame(
            [{c: flag_rows[c][j] for c in flag_cols} for j in range(k)], names)
        residuals = _typed_frame(
            _read_table(path, groups, 'residuals', list(_RESIDUAL_SCHEMA),
                        range(n_data), 'index'), _RESIDUAL_SCHEMA)
        start_cols = meta['starts_columns']
        schema = {c: _START_SCHEMA.get(c, float) for c in start_cols}
        starts = _typed_frame(
            _read_table(path, groups, 'starts', start_cols,
                        range(meta['n_starts']), 'index'), schema)
        return cls(names=names, targets=meta['targets'],
                   experiments=meta['experiments'], kinds=meta['kinds'],
                   log=np.array(meta['log'], dtype=bool), lower=lower,
                   upper=upper, x=x, cost=meta['cost'], chi2=meta['chi2'],
                   n_data=n_data, jacobian=jac, active_mask=active,
                   fd_flags=fd_flags, residuals=residuals, starts=starts,
                   n_failed_evals=meta['n_failed_evals'],
                   success=meta['success'], message=meta['message'],
                   method=meta['method'], loss=meta['loss'],
                   f_scale=meta['f_scale'],
                   absolute_sigma=meta['absolute_sigma'], seed=meta['seed'],
                   units=meta['units'], rank_rtol=meta['rank_rtol'])


def _read_table(path, groups, quantity, columns, keys, axis,
                kind=_CSV_FORMAT):
    """``{column: [values]}`` of a table stored as ``(key=column, axis=row
    key)`` rows of a long CSV written by ``kind``; ``ValueError`` unless
    every cell appears exactly once."""
    keys = [str(k) for k in keys]
    pos = {k: i for i, k in enumerate(keys)}
    g = groups.get(quantity, _EMPTY)
    expected = len(columns) * len(keys)
    if len(g) != expected:
        raise _incomplete(path, kind,
                          f'{quantity}: expected {expected} rows, found '
                          f'{len(g)}')
    out = {c: [None] * len(keys) for c in columns}
    seen = set()
    for c, key, value in zip(g['key'], g[axis], g['value']):
        if c not in out or key not in pos or (c, key) in seen:
            raise _incomplete(path, kind,
                              f'{quantity}: unexpected or duplicate cell '
                              f'({c}, {key})')
        seen.add((c, key))
        out[c][pos[key]] = json.loads(value)
    return out


def _no_values(results, values):
    """An ``fn`` for :func:`_fd.evaluate` when only the runs are wanted."""
    return np.zeros(0)


def _quantiles(q):
    """Validated percentiles (tuple of float in ``[0, 100]``, sorted,
    unique)."""
    if isinstance(q, (str, bytes)) or not hasattr(q, '__iter__'):
        q = [q]
    values = [_finite('q', v) for v in q]
    if not values or any(not 0.0 <= v <= 100.0 for v in values):
        raise ValueError(f'q must hold percentiles in [0, 100]; got {q!r}.')
    if len(set(values)) != len(values):
        raise ValueError(f'Duplicate percentiles in {q!r}.')
    return tuple(sorted(values))


def _quantile_name(value):
    """Column name of a percentile: ``5 -> 'q5'``, ``2.5 -> 'q2.5'``."""
    return f'q{value:g}'


def _warn_observable_outputs(problem, outputs):
    """Warn (for :meth:`FitResult.predict`) about outputs named like one of
    the problem's observables (shared or an experiment's own) that records
    another selection: outputs are model selections."""
    recorded = {}
    for obs in _problem_observables(problem):
        recorded.setdefault(obs.name, set()).add(obs.selection)
    for out in outputs:
        other = sorted(s for s in recorded.get(out, ()) if s != out)
        if other:
            sel = ' or '.join(map(repr, other))
            # user -> predict -> here
            warnings.warn(f'output {out!r} is read as the model selection '
                          f'{out!r}, not as the observable {out!r}, which '
                          f'records {sel}; pass {sel} for the observable\'s '
                          'model quantity.', UserWarning, stacklevel=3)


def _warn_lack_of_fit(fit, what, stacklevel):
    """Warn once that ``fit`` is flagged ``lack_of_fit``: ``what`` (profile
    intervals or prediction bands) assume that the model matches the data.
    Only this flag is propagated: ``no_covariance`` already raises in
    :meth:`FitResult.sample`, and ``at_bound`` is handled where it matters
    (the profile grid stops at the bound)."""
    flags = fit.flags
    hit = flags[flags['flag'] == 'lack_of_fit']
    if hit.empty:
        return
    value = float(hit['value'].iloc[0])
    stat = 'reduced chi2' if fit.loss == 'linear' else '2*cost/dof'
    warnings.warn(
        f'This fit is flagged lack_of_fit ({stat} = {value:.3g}): the {what} '
        'assume that the model matches the data within the given standard '
        'deviations. They are local to this optimum and may be far too '
        'narrow. See fit.summary(); if the optimum is a wrong local minimum, '
        "refit with more starts or method='hybrid'.",
        UserWarning, stacklevel=stacklevel)


# --- profile likelihood ------------------------------------------------------
def _fixed_spec(model, spec, j, value):
    """``spec`` with free parameter ``j`` fixed at the natural ``value``:
    a model target is pinned into the conditions it applies to, an
    observable scale/offset replaces its name in the blocks."""
    keep = [k for k in range(len(spec.names)) if k != j]
    name, target, applies = spec.names[j], spec.targets[j], spec.applies[j]
    conditions, blocks = spec.conditions, spec.blocks
    if target is not None:
        conditions = tuple(
            _fd.pin(model, c, {target: value})
            if applies is None or i in applies else c
            for i, c in enumerate(conditions))
    else:
        blocks = tuple(tuple(dataclasses.replace(
            b, scale=value if b.scale == name else b.scale,
            offset=value if b.offset == name else b.offset) for b in bs)
            for bs in blocks)
    return dataclasses.replace(
        spec, names=tuple(spec.names[k] for k in keep),
        targets=tuple(spec.targets[k] for k in keep), log=spec.log[keep],
        applies=tuple(spec.applies[k] for k in keep), lo=spec.lo[keep],
        hi=spec.hi[keep], conditions=conditions, blocks=blocks)


def _profile_point(model, spec, j, xj, start):
    """Re-optimize every free parameter but ``j`` with ``j`` fixed at
    ``xj`` (fitting space), from ``start`` (the others, fitting space)."""
    value = math.exp(xj) if spec.log[j] else xj
    if len(spec.names) == 1:
        evaluator = _Evaluator(model, spec)
        r, ok = evaluator.evaluate(np.array([xj]))
        return {'x': np.zeros(0), 'ok': ok, 'success': ok,
                'cost': spec.cost(r) if ok else math.nan,
                'chi2': float(np.sum(r * r)) if ok else math.nan,
                'nfev': 1, 'n_failed': evaluator.n_failed,
                'message': '' if ok else f'The point failed: '
                                         f'{evaluator.error}'}
    sub = _fixed_spec(model, spec, j, value)
    return _run_start(_Evaluator(model, sub), _clip(start, sub.lo, sub.hi),
                      {'max_nfev': None})


def _profile_slopes(fit):
    """Per free parameter ``j``, the tangent ``d(others)/dx_j`` of the
    profile path at the optimum (fitting space): the linearized
    re-optimization of the others, ``-lstsq(J_others, J_j)``, from the
    fit's weighted-residual Jacobian. At full rank it is ``Cov[others, j] /
    Cov[j, j]``; on a rank-deficient direction it is still the
    compensating move (``-1`` for ``k2`` of a product ``k1 * k2`` in log
    space), where the pseudo-inverse covariance is not. Zeros (no
    prediction) if the Jacobian is not finite."""
    k = len(fit.names)
    slopes = [np.zeros(k - 1) for _ in range(k)]
    if k < 2 or not np.all(np.isfinite(fit.jacobian)):
        return slopes
    jac = np.where(fit._unseen, 0.0, fit.jacobian)
    for j in range(k):
        try:
            with np.errstate(all='ignore'):
                s = -np.linalg.lstsq(np.delete(jac, j, axis=1), jac[:, j],
                                     rcond=None)[0]
        except np.linalg.LinAlgError:
            continue
        if np.all(np.isfinite(s)):
            slopes[j] = s
    return slopes


def _predicted_start(path, xj, slope):
    """Start for the others at ``xj``: extrapolated linearly along the
    ``path`` of re-optimized ``(x_j, others)`` pairs (in the order run),
    through its last two points (a secant), or along ``slope`` from its
    only point (the optimum)."""
    x1, s1 = path[-1]
    if len(path) == 1:
        return s1 + slope * (xj - x1)
    x0, s0 = path[-2]
    return s1 + (s1 - s0) * ((xj - x1) / (x1 - x0))


def _profile_side(model, item, *, spec, threshold):
    """One side of one profile: re-optimize at each grid point outward
    from the optimum, each started from a prediction of the others
    (along the fit's tangent ``slope`` for the first point, then along
    the secant through the last two re-optimized points); past the grid,
    continue with doubling steps until the threshold is crossed or the
    bound is reached, then bisect the last bracket (below, above the
    threshold) down to the grid step, so that an interval end is never
    interpolated over more than one grid step. Returns ``[(xj, record)]``
    in the order run."""
    j, xs, step, bound, extend, start, cost0, scale, x0, slope = item
    base = step
    out = []
    # the re-optimized points of this side so far, from the optimum out
    path = [(float(x0), np.asarray(start, dtype=float))]
    last_ok = None      # the last point that re-optimized (its x)

    def run(xj, warm):
        rec = _profile_point(model, spec, j, xj, warm)
        out.append((float(xj), rec))
        return rec

    def above(rec):
        return rec['ok'] and 2.0 * (rec['cost'] - cost0) / scale >= threshold

    crossed = False
    for xj in xs:
        rec = run(xj, _predicted_start(path, xj, slope))
        if rec['ok']:
            last_ok = xj
            path.append((float(xj), rec['x']))
        crossed = crossed or above(rec)
    xj = xs[-1] if xs else None
    while extend and not crossed and xj != bound:
        # no point so far crossed the threshold: the last ok one is the
        # inner end of the next bracket
        inner, inner_path = last_ok, path[-2:]
        step *= 2.0
        xj = (min(xj + step, bound) if bound > xj else max(xj - step, bound))
        rec = run(xj, _predicted_start(path, xj, slope))
        if rec['ok']:
            last_ok = xj
            path.append((float(xj), rec['x']))
        crossed = above(rec)
        if not crossed or inner is None:
            continue
        # bisect, each point predicted from the inner end's path
        outer = xj
        while abs(outer - inner) > base * (1.0 + 1e-9):
            mid = 0.5 * (inner + outer)
            rec = run(mid, _predicted_start(inner_path, mid, slope))
            if not rec['ok']:
                break
            if above(rec):
                outer = mid
            else:
                inner = mid
                inner_path = [inner_path[-1], (float(mid), rec['x'])]
    return out


def _profile_grid(x, se, lo, hi, n_points, span, active):
    """Per side (lower, upper): ``(points, step, bound, extend)`` in the
    fitting space. ``span * se`` (clipped to the bound) is split into
    equal steps; with a non-finite or zero ``se`` the grid runs to the
    bound. A side whose bound the estimate is at (``active``: the
    optimizer's active mask, -1 or 1) is empty."""
    n_lower = (n_points - 1) // 2
    sides = []
    for sign, n_side, bound in ((-1.0, n_lower, lo),
                                (1.0, n_points - 1 - n_lower, hi)):
        room = 0.0 if active == sign else abs(bound - x)
        if np.isfinite(se) and se > 0:
            extent = min(span * se, room)
        else:
            extent = room
        if extent <= 0 or n_side == 0:
            sides.append(([], 0.0, bound, False))
            continue
        step = extent / n_side
        points = [x + sign * step * k for k in range(1, n_side + 1)]
        if extent == room:
            points[-1] = bound
        sides.append((points, step, bound, extent < room))
    return sides


def profile_likelihood(fit, names, *, n_points=21, span=3.0, workers=1):
    """Profile likelihoods of fitted parameters.

    For each parameter, a grid of fixed values runs outward from the
    estimate on both sides; at each value every other free parameter is
    re-optimized (``least_squares``), started from a prediction: for the
    first value, the optimum moved along the fit's tangent (the
    linearized re-optimization from its Jacobian, ``Cov[others, j] /
    Cov[j, j]`` at full rank), then the line through the solutions at
    the two previous values. The profile is
    ``delta = 2 * (cost - cost_min) / s2`` with ``cost`` the problem's
    objective ``0.5 * sum(rho(r**2))`` -- for the linear loss the increase
    of the weighted sum of squared residuals (chi-square) -- and ``s2 = 1``
    for absolute sigmas or the fit's reduced chi-square otherwise.

    Parameters
    ----------
    fit : FitResult
        A result of :func:`fit_parameters` (with its problem; not one read
        by :meth:`FitResult.from_csv`).
    names : str or sequence of str
        Free parameters to profile.
    n_points : int
        Grid values per parameter, the estimate included (>= 3; default
        21, i.e. 10 per side).
    span : float
        Half-width of the grid in standard errors (fitting space; default
        3.0), clipped to the bounds. With an infinite or undefined SE the
        grid spans the whole range to each bound, so its step is the range
        to the bound divided by the points per side: an end is then
        interpolated over that coarse step and is not exact at a kink
        (such as a re-optimized parameter reaching its bound); raise
        ``n_points`` for a sharper end. A side that has not
        crossed the threshold at ``span`` SE continues outward with
        doubling steps until it crosses or reaches the bound; the last
        step's bracket is then bisected down to the grid step, so every
        interval end is interpolated over at most one grid step (the
        profile may then have more rows than ``n_points``). An estimate
        at a bound (``FitResult.active_mask``) has no grid on that side.
    workers : int
        Processes over which the sides of the profiles run (each side is
        one sequential chain); ``1`` (default) runs serially. Under the
        ``spawn`` start method (macOS, Windows) every observable
        ``transform`` must be picklable (a module-level function, not a
        lambda), and a script must call this under
        ``if __name__ == '__main__':``.

    Returns
    -------
    dict
        ``{name: ProfileResult}``, in the order of ``names``.

    Raises
    ------
    TypeError, ValueError
        For an invalid argument, an unknown or repeated name, a ``fit``
        without its problem, ``absolute_sigma=False`` without degrees of
        freedom, or a robust loss with ``absolute_sigma=False``; before any
        run.

    Warns
    -----
    UserWarning
        If a profile goes below the fit's optimum by more than
        ``PROFILE_DROP_WARN`` (0.01, in units of ``delta``): the fit had
        not reached the minimum, and the lowest cost found is used as
        ``cost_min``. Also once per interval end interpolated across a
        step in which ``delta`` rose by more than the threshold beyond the
        Wald parabola (:attr:`ProfileResult.abrupt_ends`): the profile may
        jump there and that end may be too narrow (see Notes), or the grid
        is too coarse for the profile's shape. Also, before any run, if the
        fit is flagged ``lack_of_fit`` (:attr:`FitResult.flags`): the
        intervals assume that the model matches the data, are local to
        this optimum and may be far too narrow.

    Notes
    -----
    The 95% interval is where ``delta <= chi2.ppf(0.95, 1) = 3.84``
    (``PROFILE_THRESHOLD``) for absolute sigmas. With
    ``absolute_sigma=False`` the noise level is estimated, so the threshold
    is ``F.ppf(0.95, 1, dof) = t.ppf(0.975, dof)**2`` (the square of the
    Wald interval's quantile; e.g. 6.61 at ``dof = 5``). Its ends are
    interpolated between the last grid value below and the first at or
    above the threshold, moving outward from the profile's lowest point
    (the estimate, when the fit is at the minimum), linearly in
    ``sqrt(delta)`` against the fitting-space value (exact for a
    quadratic profile). A side that does not cross within the bounds
    leaves the interval open there and the parameter is not identifiable.
    The converse does not hold: crossings are relative to every
    parameter's bounds, and a side can cross only because a re-optimized
    parameter reaches its bound (its ``'estimate:<other>'`` column equals
    that bound, e.g. ``k1`` of a product ``k1 * k2`` limited by the box of
    ``k2``). The interval then reflects the box constraints, not the
    data; read it together with the fit's flags (``rank_deficient``,
    ``high_correlation``, an infinite ``se``). Unlike the Wald interval
    this follows the actual objective, so it is reliable for skewed or
    flat profiles.

    The re-optimization is local. With fed-batch data the objective jumps
    where an event crosses a sample time, and a re-optimization can stay
    on the far side of such a jump although other values of the free
    parameters would avoid it; the predicted start follows the valley of
    the objective and avoids most of these, but not all. The profile then
    rises by the size of the jump within one step and the interval ends
    there, too narrow: such an end is flagged in
    :attr:`ProfileResult.abrupt_ends` with a warning. Check it by fixing
    the parameter just beyond that end in the conditions and refitting
    the others with ``method='hybrid'``, or raise ``n_points`` if the
    grid was merely too coarse.

    With a robust loss the threshold is not calibrated: the robust cost
    rises more slowly than chi-square, even for residuals within
    ``f_scale``, so the intervals come out wider than the linear-loss
    ones. For example, on 12 outlier-free Gaussian points with
    ``f_scale=1``, the profile intervals of a decay rate and its initial
    value were about 1.1 times as wide for ``'huber'`` (quadratic within
    ``f_scale``), 1.3 times for ``'soft_l1'`` and 1.6-1.8 times for
    ``'cauchy'`` as the linear-loss ones; the factor shrinks as
    ``f_scale`` grows. Use a robust-loss profile to see the shape of the
    objective and locate the optimum, and for intervals refit with the
    linear loss from the robust estimates (``fit_parameters(linear_problem,
    initial=robust_fit)``, with ``linear_problem`` the same
    :class:`CalibrationProblem` built with ``loss='linear'``, the
    default) and profile that fit.

    Each grid value costs a full re-optimization (about
    ``1 + 2 (n_free - 1)`` runs per experiment and iteration).
    """
    if not isinstance(fit, FitResult):
        raise TypeError(f'fit must be a FitResult; got {type(fit).__name__}.')
    if fit.problem is None:
        raise ValueError('This FitResult has no problem (it was read from a '
                         'file); refit to profile it.')
    if isinstance(names, str):
        names = [names]
    if not hasattr(names, '__iter__'):
        raise TypeError('names must be a free-parameter name or a sequence '
                        'of them.')
    names = list(names)
    unknown = [n for n in names if n not in fit.names]
    if not names or unknown:
        raise ValueError(f'names must be free parameters of the fit '
                         f'{fit.names}; got {names}.')
    CalibrationProblem._check_names('profile', names)
    n_points = _count('n_points', n_points, 3)
    span = _positive('span', span)
    workers = _ex._check_workers(workers)
    if fit.loss != 'linear' and not fit.absolute_sigma:
        raise ValueError(f'A fit with loss={fit.loss!r} needs absolute sigmas '
                         'to be profiled (see CalibrationProblem).')
    scale = 1.0 if fit.absolute_sigma else fit.reduced_chi2
    threshold = (PROFILE_THRESHOLD if fit.absolute_sigma
                 else _wald_quantile(False, fit.dof) ** 2)
    if not (np.isfinite(scale) and scale > 0):
        raise ValueError('With absolute_sigma=False the profile is scaled by '
                         'the reduced chi-square, which needs n_data > n_free '
                         'and a nonzero chi2.')
    # user -> profile_likelihood -> the helper
    _warn_lack_of_fit(fit, 'profile intervals', stacklevel=3)
    problem = fit.problem
    spec, private = problem._spec, problem._ready()
    x, se = fit.x, fit.se.to_numpy()
    slopes = _profile_slopes(fit)
    items = []
    for name in names:
        j = fit.names.index(name)
        others = np.delete(x, j)
        for points, step, bound, extend in _profile_grid(
                x[j], se[j], spec.lo[j], spec.hi[j], n_points, span,
                fit.active_mask[j]):
            items.append((j, points, step, bound, extend, others, fit.cost,
                          scale, x[j], slopes[j]))
    task = functools.partial(_profile_side, spec=spec, threshold=threshold)
    if workers == 1:
        sides = [task(private, item) for item in items]
    else:
        sides = _ex._map_with_model(
            private, task, items, workers=workers, chunksize=1,
            spawn_check=functools.partial(_check_transforms_picklable,
                                          problem))
    out = {}
    for i, name in enumerate(names):
        j = fit.names.index(name)
        others = [n for n in fit.names if n != name]
        center = {'side': 'center', 'x': float(x[j]), 'cost': fit.cost,
                  'chi2': fit.chi2, 'status': 'ok', 'success': fit.success,
                  'nfev': 0, 'n_failed_evals': 0, 'message': fit.message,
                  'others': np.delete(x, j)}
        rows = [center]
        for side, chain in zip(('lower', 'upper'), sides[2 * i:2 * i + 2]):
            for xj, rec in chain:
                rows.append({'side': side, 'x': xj, 'cost': rec['cost'],
                             'chi2': rec['chi2'],
                             'status': 'ok' if rec['ok'] else 'failed',
                             'success': bool(rec['success']),
                             'nfev': int(rec['nfev']),
                             'n_failed_evals': int(rec.get('n_failed', 0)),
                             'message': rec['message'],
                             'others': rec['x'] if rec['ok']
                             else np.full(len(others), math.nan)})
        rows.sort(key=lambda row: row['x'])
        log = np.delete(spec.log, j)
        columns = {c: [row[c] for row in rows] for c in _PROFILE_SCHEMA
                   if c != 'value'}
        columns['value'] = [math.exp(row['x']) if spec.log[j] else row['x']
                            for row in rows]
        for k, other in enumerate(others):
            columns[f'estimate:{other}'] = [
                float(np.exp(row['others'][k]) if log[k]
                      else row['others'][k]) for row in rows]
        schema = {**_PROFILE_SCHEMA,
                  **{f'estimate:{o}': float for o in others}}
        profile = _typed_frame({c: columns[c] for c in schema}, schema)
        result = ProfileResult(
            name=name, log=bool(spec.log[j]), x=float(x[j]),
            se=float(se[j]), lower=float(fit.lower[j]),
            upper=float(fit.upper[j]), cost=fit.cost, scale=scale,
            profile=profile, threshold=threshold, loss=fit.loss,
            units=dict(fit.units))
        drop = 2.0 * (fit.cost - result.reference_cost) / scale
        if drop > PROFILE_DROP_WARN:
            warnings.warn(
                f'The profile of {name!r} reached a cost below the fit\'s '
                f'optimum (delta = {-drop:.3g}): the fit had not found the '
                'minimum. Its lowest cost is used as the reference; refit '
                '(e.g. more starts) before trusting the intervals.',
                UserWarning, stacklevel=2)
        for side, step in zip(('lower', 'upper'), result._abrupt_steps):
            if step is None:
                continue
            a, da, b, db = step
            warnings.warn(
                f'The {side} end of the profile interval of {name!r} is not '
                f'resolved: delta rose from {da:.3g} to {db:.3g} within one '
                f'grid step ({a:.6g} to {b:.6g}), by more than the threshold '
                f'({threshold:.3g}) beyond the Wald parabola. The profile may '
                'jump there: the re-optimization of the other parameters '
                'may have stayed on the far side of a jump in the objective '
                '(e.g. an event crossing a sample time), and the end may be '
                'too narrow. Check it by fixing the parameter beyond that '
                'end in the conditions and refitting the others with '
                "method='hybrid'; if the profile is merely steep there "
                '(e.g. next to a bound), raise n_points.',
                UserWarning, stacklevel=2)
        out[name] = result
    return out


@dataclass(eq=False, repr=False)
class ProfileResult:
    """Profile likelihood of one parameter (from :func:`profile_likelihood`).

    Built from the fields below (what :meth:`to_csv` stores); the
    attributes are derived from them.

    Parameters
    ----------
    name : str
        The profiled free parameter.
    log : bool
        Whether it is fitted in log space (the grid is uniform there).
    x : float
        The fit's estimate, in the fitting space.
    se : float
        Its Wald standard error in the fitting space (``inf`` if not
        identifiable at the optimum).
    lower, upper : float
        Its bounds (natural units).
    cost : float
        The fit's objective ``0.5 * sum(rho(r**2))`` at the optimum.
    scale : float
        ``s2``: 1 for absolute sigmas, else the fit's reduced chi-square.
    profile : pandas.DataFrame
        One row per grid value, sorted by value: ``side`` (``'lower'``,
        ``'center'``, ``'upper'``), ``x`` (fitting space), ``value``
        (natural units), ``cost`` and ``chi2`` after re-optimizing the
        others (``NaN`` if it failed), ``status`` (``'ok'|'failed'``),
        ``success`` (the optimizer converged), ``nfev``, ``n_failed_evals``,
        ``message`` and ``'estimate:<name>'`` of each other parameter
        (natural units). A ``delta`` column is added (see Attributes).
    threshold : float
        Threshold on ``delta``: by default ``PROFILE_THRESHOLD``
        (``chi2.ppf(0.95, 1) = 3.84``); :func:`profile_likelihood` gives
        ``t.ppf(0.975, dof)**2`` with ``absolute_sigma=False``.
    loss : str
        The fit's loss.
    units : dict
        The model's units.

    Attributes
    ----------
    reference_cost : float
        ``cost_min``: the lowest of ``cost`` and the profile's costs.
    confidence_interval : (float, float)
        Interval in natural units where ``delta <= threshold``, searched
        outward on each side from the profile's lowest point (``delta =
        0``: the estimate when the fit is at the minimum), its ends
        interpolated between grid values (linearly in ``sqrt(delta)``,
        exact for a quadratic profile); an end whose side does not cross
        within the bounds is open: ``0`` (log) or ``-inf`` below, ``inf``
        above (as :attr:`FitResult.ci` reports an infinite SE). When the
        profile finds a lower cost than the fit's,
        :func:`profile_likelihood` warns, and the estimate may lie outside
        the interval.
    is_identifiable : bool
        Whether both sides cross the threshold within the bounds.
        Crossings are relative to every parameter's bounds: a side can
        cross only because a re-optimized parameter reaches its bound (its
        ``'estimate:<other>'`` column in ``profile`` equals that bound,
        e.g. ``k1`` of a product ``k1 * k2`` limited by the box of
        ``k2``). The interval then reflects the box constraints, not the
        data; read it together with the fit's flags (``rank_deficient``,
        ``high_correlation``, an infinite ``se``).
    abrupt_ends : (bool, bool)
        Per side (lower, upper), whether its interval end is interpolated
        across a step in which ``delta`` rose by more than ``threshold``
        beyond the rise of the Wald parabola ``((x - x_hat) / se)**2``
        over the same step (by more than ``threshold`` with a non-finite
        ``se``). Such an end is not resolved by the grid
        (:func:`profile_likelihood` warns): either the profile jumps
        there, e.g. because the re-optimization stayed on the far side of
        an event crossing a sample time, and the end may be too narrow, or
        the grid is too coarse for the profile's shape (e.g. a steep rise
        next to a bound). ``False`` for an open side.
    n_failed : int
        Grid values whose re-optimization failed.

    Raises
    ------
    ValueError
        If the profile lacks a column.
    """
    name: str
    log: bool
    x: float
    se: float
    lower: float
    upper: float
    cost: float
    scale: float
    profile: pd.DataFrame
    threshold: float = PROFILE_THRESHOLD
    loss: str = 'linear'
    units: dict = field(default_factory=dict)

    def __post_init__(self):
        self.name, self.log, self.loss = str(self.name), bool(self.log), \
            str(self.loss)
        for arg in ('x', 'se', 'lower', 'upper', 'cost', 'scale',
                    'threshold'):
            setattr(self, arg, float(getattr(self, arg)))
        missing = [c for c in _PROFILE_SCHEMA if c not in self.profile]
        if missing:
            raise ValueError(f'profile lacks the column(s) {missing}.')
        frame = self.profile.drop(columns='delta', errors='ignore')
        costs = frame['cost'].to_numpy(dtype=float)
        with np.errstate(all='ignore'):
            self.reference_cost = float(np.nanmin(np.append(costs,
                                                            self.cost)))
            delta = 2.0 * (costs - self.reference_cost) / self.scale
        frame.insert(frame.columns.get_loc('chi2') + 1, 'delta', delta)
        self.profile = frame.reset_index(drop=True)
        self.n_failed = int((frame['status'] != 'ok').sum())
        ends, closed = [], []
        x = frame['x'].to_numpy(dtype=float)
        ok = (frame['status'] == 'ok').to_numpy()
        # The search starts where delta = 0: the estimate when the fit is at
        # the minimum, else the profile's lowest point (a fit that had not
        # reached the minimum must not give an empty interval around it).
        x0 = self.x
        d0 = 2.0 * (self.cost - self.reference_cost) / self.scale
        lowest = np.flatnonzero(ok & np.isfinite(delta))
        if lowest.size and np.min(delta[lowest]) < d0:
            i0 = lowest[np.argmin(delta[lowest])]
            x0, d0 = x[i0], delta[i0]
        wald = np.isfinite(self.se) and self.se > 0
        steps = []
        for sign in (-1.0, 1.0):
            # outward from x0: x0 itself, then this side's values
            side = np.flatnonzero(ok & (sign * (x - x0) > 0))
            side = side[np.argsort(sign * (x[side] - x0), kind='stable')]
            prev_x, prev_d = x0, d0
            end = step = None
            for i in side:
                if delta[i] >= self.threshold:
                    # linear in sqrt(delta), exact for a quadratic profile
                    # (prev_d < threshold <= delta[i]: a positive divisor)
                    root = math.sqrt(max(prev_d, 0.0))
                    w = ((math.sqrt(self.threshold) - root)
                         / (math.sqrt(delta[i]) - root))
                    end = prev_x + w * (x[i] - prev_x)
                    # the crossing step: did delta rise by more than the
                    # threshold beyond the Wald parabola's rise over it?
                    rise = 0.0
                    if wald:
                        with np.errstate(all='ignore'):
                            rise = float(((x[i] - self.x) ** 2
                                          - (prev_x - self.x) ** 2)
                                         / self.se ** 2)
                        rise = rise if math.isfinite(rise) else 0.0
                    if delta[i] - prev_d > self.threshold + max(rise, 0.0):
                        a, b = (math.exp(v) if self.log else float(v)
                                for v in (prev_x, x[i]))
                        step = (a, float(prev_d), b, float(delta[i]))
                    break
                prev_x, prev_d = x[i], delta[i]
            closed.append(end is not None)
            steps.append(step)
            if end is None:
                end = sign * math.inf
            ends.append(math.exp(end) if self.log else end)
        self.confidence_interval = (ends[0], ends[1])
        self.is_identifiable = all(closed)
        # (value, delta) before and after each abrupt crossing step
        self._abrupt_steps = tuple(steps)
        self.abrupt_ends = tuple(step is not None for step in steps)

    def __repr__(self):
        lo, hi = self.confidence_interval
        return (f'ProfileResult({self.name!r}, estimate='
                f'{self.estimate:.6g}, confidence_interval=({lo:.6g}, '
                f'{hi:.6g}), is_identifiable={self.is_identifiable}, '
                f'{len(self.profile)} points)')

    @property
    def estimate(self):
        """The fit's estimate in natural units (float)."""
        return math.exp(self.x) if self.log else self.x

    def plot(self, ax=None):
        """Plot ``delta`` against the parameter value.

        The profile (failed grid values left out) with the threshold and
        the interval ends dashed, the estimate marked at its own ``delta``
        (0 when the fit is at the minimum), and the Wald approximation
        ``((x - x_hat) / se)**2`` dotted for comparison (if
        the SE is finite). Log parameters get a log value axis. The
        vertical axis is labelled ``Delta chi2`` for the linear loss with
        absolute sigmas, else ``2 Delta cost`` (divided by ``s^2``, the
        reduced chi-square, with ``absolute_sigma=False``).

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure otherwise.

        Returns
        -------
        (matplotlib.figure.Figure, matplotlib.axes.Axes)
        """
        _set_font()
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        frame = self.profile[self.profile['status'] == 'ok']
        ax.plot(frame['value'], frame['delta'], marker='o', ms=3,
                color='tab:blue', label='profile')
        if np.isfinite(self.se) and self.se > 0:
            x = frame['x'].to_numpy(dtype=float)
            ax.plot(frame['value'], ((x - self.x) / self.se) ** 2,
                    color='gray', ls=':', label='Wald')
        ax.axhline(self.threshold, color='black', ls='--', lw=1,
                   label=f'{PROFILE_LEVEL:.0%} threshold')
        for end in self.confidence_interval:
            if np.isfinite(end) and not (self.log and end <= 0):
                ax.axvline(end, color='tab:red', ls='--', lw=1)
        # at its own delta: above 0 if the profile found a lower cost
        ax.plot([self.estimate],
                [2.0 * (self.cost - self.reference_cost) / self.scale],
                marker='*', ms=9, color='tab:red', ls='none',
                label='estimate')
        ax.set_xlabel(self.name)
        if self.loss == 'linear' and self.scale == 1.0:
            ylabel = r'$\Delta\chi^2$'
        elif self.scale == 1.0:
            ylabel = r'$2\,\Delta\mathrm{cost}$'
        else:
            ylabel = r'$2\,\Delta\mathrm{cost}/s^2$'
        ax.set_ylabel(ylabel)
        _style_axes(ax)
        if self.log:
            ax.set_xscale('log')     # replaces the linear minor locator
            lo, hi = ax.get_xlim()
            if hi < 10.0 * lo:
                # Within a decade the log locator labels crowded minor ticks
                # in scientific notation; label a few plain values instead.
                ax.xaxis.set_major_locator(MaxNLocator(4))
                ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
                ax.xaxis.set_minor_locator(NullLocator())
        ax.legend(fontsize=8, frameon=False)
        return fig, ax

    def to_csv(self, path):
        """Write the result to a long-format CSV.

        Columns ``quantity, key, parameter, index, value``: ``meta`` rows
        (``key``: field) and the profile table (``key``: column, ``index``:
        row; the derived ``delta`` is not stored). Values are JSON, so
        :meth:`from_csv` restores the result exactly.

        Parameters
        ----------
        path : str or path-like
            Destination; an existing file is overwritten.

        Returns
        -------
        None
        """
        frame = self.profile.drop(columns='delta')
        meta = {'nskinetics': _PROFILE_FORMAT, 'name': self.name,
                'log': self.log, 'x': self.x, 'se': self.se,
                'lower': self.lower, 'upper': self.upper, 'cost': self.cost,
                'scale': self.scale, 'threshold': self.threshold,
                'loss': self.loss, 'units': self.units,
                'columns': list(frame.columns), 'n_rows': len(frame)}
        rows = [('meta', key, '', '', json.dumps(_py(value)))
                for key, value in meta.items()]
        for c in frame.columns:
            rows += [('profile', c, '', i, json.dumps(_py(v)))
                     for i, v in enumerate(frame[c].tolist())]
        pd.DataFrame(rows, columns=list(_CSV_COLUMNS)).to_csv(
            path, index=False, lineterminator='\n', encoding='utf-8')

    @classmethod
    def from_csv(cls, path):
        """Read a result written by :meth:`to_csv`.

        Parameters
        ----------
        path : str or path-like

        Returns
        -------
        ProfileResult
            Equal in every field and attribute.

        Raises
        ------
        ValueError
            If the file was not written by :meth:`to_csv` or is incomplete.
        """
        meta, groups = _read_long_csv(path, _PROFILE_FORMAT, required=(
            'name', 'log', 'x', 'se', 'lower', 'upper', 'cost', 'scale',
            'threshold', 'loss', 'units', 'columns', 'n_rows'))
        columns = meta['columns']
        schema = {c: _PROFILE_SCHEMA.get(c, float) for c in columns}
        profile = _typed_frame(
            _read_table(path, groups, 'profile', columns,
                        range(meta['n_rows']), 'index', kind=_PROFILE_FORMAT),
            schema)
        return cls(name=meta['name'], log=meta['log'], x=meta['x'],
                   se=meta['se'], lower=meta['lower'], upper=meta['upper'],
                   cost=meta['cost'], scale=meta['scale'], profile=profile,
                   threshold=meta['threshold'], loss=meta['loss'],
                   units=meta['units'])
