# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Event-aware local sensitivity analysis.

:func:`local_sensitivity` differentiates output time courses and scalar
metrics with respect to parameters and initial values by central (or
forward) finite differences in log-parameter space, each run restarting from
the restored baseline on a private copy of the model, and reports per
parameter whether a perturbation changed how often an event fired (see
:mod:`nskinetics.analysis._fd`). RoadRunner's own forward sensitivities are
not used: they drop the jumps events make in the state.
"""
import json
import math
import numbers
import warnings
from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from ..engine.conditions import Condition
from ..engine.kinetic_model import KineticModel
from ..engine import execution as _ex
from ..engine._plotting import _bf, _set_font, _style_axes
from . import _fd
from ._common import (_CSV_COLUMNS, _check_ids, _count, _id_list,
                      _incomplete, _read_long_csv, _scatter)
from .metrics import _check_metrics_condition

__all__ = ('LocalSensitivityResult', 'local_sensitivity')

_SCALINGS = ('log', 'semi', 'none')
#: The column of :meth:`LocalSensitivityResult.rank` each scaling sorts by.
_RANK_KEY = {'log': 'elasticity', 'semi': 'semi_elasticity',
             'none': 'derivative'}
_CSV_FORMAT = 'LocalSensitivityResult'
#: A log-scaled sensitivity of an output (``d ln y / d ln p``) is ``NaN``
#: where its estimated noise from the integrator's absolute tolerance,
#: ``FD_NOISE * atol_y / (|y| h)`` (``h``: the parameter's relative step;
#: ``atol_y``: the tolerance in the output's units, ``atol * |y0|`` for an
#: output starting at ``y0 != 0``, else ``atol``, since RoadRunner scales a
#: scalar tolerance by each state's start value; see
#: :data:`nskinetics.analysis._fd.FD_NOISE` and
#: :func:`~nskinetics.analysis._fd.effective_atol`), exceeds this: below
#: ``|y| = 30 atol_y / h`` (``|y| / |y0|`` below 3e-8 at the default
#: ``atol = 1e-12`` and ``h = 1e-3``; ``|y|`` below ``3e4 atol`` for an
#: output starting at 0) the value is noise, not a derivative. ``atol`` is
#: the recorded one: by default :func:`local_sensitivity` scales it to the
#: smallest peak of the outputs (and metric selections) that start at 0
#: (:func:`~nskinetics.analysis._fd.scale_atol_to_outputs`), so such an
#: output is masked only where it is small relative to the smallest such
#: peak.
LOG_NOISE_MAX = 0.1


def _scale(raw, nominal, values, scaling, floor=None):
    """Scaled sensitivities from ``raw = dy/dp``.

    ``raw`` has the parameter on its last axis, ``nominal`` the shape of
    ``raw`` without it and ``values`` one entry per parameter. ``'none'``:
    ``dy/dp``; ``'semi'``: ``p dy/dp``; ``'log'``: ``(p / y) dy/dp``, ``NaN``
    where ``y`` is 0 and, given ``floor`` (``|y|`` cut-offs that broadcast
    against ``raw``, e.g. one per output and parameter for an array
    ``(n_t, n_out, n_p)``; ``NaN`` for none), where ``|y|`` is below it.
    One function for the constructor and :meth:`from_csv`, so a round trip
    is exact.

    A log-scaled value is only as good as ``y``: where ``|y|`` sits at the
    integrator's noise level (a depleted substrate, ``Final('[S]')`` of a
    residual substrate) the difference of ``y`` is integration noise of
    size about ``FD_NOISE * atol_y / h`` and ``(p / y) dy/dp`` is noise
    divided by a tiny number (see :data:`LOG_NOISE_MAX`).
    """
    with np.errstate(all='ignore'):
        if scaling == 'none':
            return np.array(raw, dtype=float, copy=True)
        semi = raw * values
        if scaling == 'semi':
            return semi
        y = np.asarray(nominal, dtype=float)[..., None]
        out = np.where(y != 0, semi / y, np.nan)
        if floor is not None:
            out = np.where(np.abs(y) < np.asarray(floor, dtype=float),
                           np.nan, out)
        return out


@dataclass(eq=False, repr=False)
class LocalSensitivityResult:
    """Outcome of :func:`local_sensitivity`.

    Parameters
    ----------
    time : numpy.ndarray
        Output times, shape ``(n_t,)``, model time units.
    outputs : list of str
        Differentiated selections (e.g. ``'[P]'``).
    parameters : list of str
        Perturbed ids, in order.
    values : numpy.ndarray
        Nominal parameter values, shape ``(n_p,)``, model units.
    nominal : numpy.ndarray
        Nominal output trajectories ``y``, shape ``(n_t, n_out)``.
    raw : numpy.ndarray
        ``dy/dp``, shape ``(n_t, n_out, n_p)``: output units per parameter
        unit.
    metric_names : list of str
        Differentiated metrics.
    metric_values : numpy.ndarray
        Nominal metric values ``m``, shape ``(n_m,)``.
    metric_derivatives : numpy.ndarray
        ``dm/dp``, shape ``(n_m, n_p)``.
    flags : pandas.DataFrame
        Finite-difference diagnostics, one row per parameter (index
        ``parameter``): ``value``, ``log`` (differenced in ``ln p``),
        ``step`` (the step of the reported difference: relative in log
        space, else absolute), ``n_refinements``, ``event_count_changed``
        (an event counter still differed between the perturbed and nominal
        runs at some output time at the smallest step: a firing added,
        removed or moved across an output time), ``reliable`` (no such
        change and no failed run), ``scheme`` (the difference used:
        ``'central'``, ``'forward'`` or ``'backward'``) and ``error`` (why a
        column is ``NaN``, or ``None``).
    scaling : {'log', 'semi', 'none'}
        Scaling of :attr:`scaled`.
    method : {'central', 'forward'}
    rel_step : float
        Initial step in log space.
    counters : tuple of str
        Event counters that were watched.
    units : dict
        The model's units.
    condition : Condition or None
        The base condition.
    tolerances : (float, float) or None
        ``(atol, rtol)`` of the integrator the runs used (``None``: unknown,
        or an integrator without tolerances).
    max_refinements : int
        Maximum number of step refinements that was allowed.

    Attributes
    ----------
    scaled : numpy.ndarray
        Shape ``(n_t, n_out, n_p)``: ``(p / y) dy/dp = d ln y / d ln p``
        (``'log'``), ``p dy/dp = dy / d ln p`` (``'semi'``) or ``dy/dp``
        (``'none'``). Log-scaled entries are ``NaN`` where ``y`` is 0 and,
        when :attr:`tolerances` are known, where ``|y|`` is at the
        integrator's noise level: below ``30 atol_y / h`` (``h``: the
        parameter's relative step in :attr:`flags`; ``atol_y``: the
        absolute tolerance in the output's units, ``atol * |y0|`` for an
        output whose first row ``y0`` is nonzero, as RoadRunner scales a
        scalar tolerance by each state's start value, else ``atol``), where
        the estimated noise ``3 atol_y / (|y| h)`` of ``d ln y / d ln p``
        exceeds 0.1: ``|y| / |y0|`` below 3e-8 at ``atol = 1e-12`` and
        ``rel_step = 1e-3``, or, for an output starting at 0, ``|y|``
        below ``3e4`` times the ``atol`` in use (recorded in
        :attr:`tolerances`; 3e-8 when it is 1e-12). Above that cut-off the
        noise falls as ``1 / |y|`` (about 1% at ``300 atol_y / h``), so
        treat log sensitivities of outputs below about ``1000 atol_y / h``
        (1e-6 relative to ``|y0|`` at the defaults; more after step
        refinement) with care, and prefer ``'semi'`` or ``'none'`` there.
        ``flags['reliable']`` reflects event counts and failed runs only,
        not this per-entry accuracy.

    Raises
    ------
    ValueError
        If ``scaling`` is not one of ``'log'``, ``'semi'``, ``'none'``.
    """
    time: np.ndarray
    outputs: list
    parameters: list
    values: np.ndarray
    nominal: np.ndarray
    raw: np.ndarray
    metric_names: list
    metric_values: np.ndarray
    metric_derivatives: np.ndarray
    flags: pd.DataFrame
    scaling: str = 'log'
    method: str = 'central'
    rel_step: float = 1e-3
    counters: tuple = ()
    units: dict = field(default_factory=dict)
    condition: Condition | None = None
    tolerances: tuple | None = None
    max_refinements: int = 2

    def __post_init__(self):
        if self.scaling not in _SCALINGS:
            raise ValueError(f'scaling must be one of {_SCALINGS}; got '
                             f'{self.scaling!r}.')
        for name in ('time', 'values', 'nominal', 'raw', 'metric_values',
                     'metric_derivatives'):
            # A copy, so later changes to the caller's arrays cannot leave
            # the derived `scaled` stale.
            setattr(self, name, np.array(getattr(self, name), dtype=float))
        self.outputs = list(self.outputs)
        self.parameters = list(self.parameters)
        self.metric_names = list(self.metric_names)
        self.counters = tuple(self.counters)
        self.units = dict(self.units)
        if self.tolerances is not None:
            self.tolerances = tuple(float(v) for v in self.tolerances)
        self.max_refinements = int(self.max_refinements)
        self.scaled = _scale(self.raw, self.nominal, self.values, self.scaling,
                             self._log_floor())

    def __repr__(self):
        n_bad = int((~self.flags['reliable']).sum())
        return (f'LocalSensitivityResult(outputs={self.outputs}, metrics='
                f'{self.metric_names}, parameters={self.parameters}, '
                f'n_t={len(self.time)}, scaling={self.scaling!r}, '
                f'unreliable={n_bad})')

    def _log_floor(self):
        """Per output and parameter (shape ``(n_out, n_p)``), the ``|y|``
        below which a log-scaled output sensitivity is noise
        (``FD_NOISE * atol_y / (LOG_NOISE_MAX * h)``, with ``atol_y`` the
        absolute tolerance in the output's units from its first nominal
        row; ``NaN`` -- no cut-off -- without tolerances or a relative
        step). From the recorded tolerances, nominal values and steps only,
        so :meth:`from_csv` restores it exactly."""
        n_out = len(self.outputs)
        floor = np.full((n_out, len(self.parameters)), np.nan)
        if (self.tolerances is None or not len(self.flags) or not n_out
                or not len(self.nominal)):
            return floor
        atol = _fd.effective_atol(self.tolerances[0], self.nominal[0])
        flags = self.flags
        for j, (value, log, step) in enumerate(zip(
                flags['value'], flags['log'], flags['step'])):
            h = _fd.relative_step(value, log, step)
            if h > 0:                                   # False for NaN
                with np.errstate(over='ignore'):    # inf at a subnormal h
                    floor[:, j] = _fd.FD_NOISE * atol / (LOG_NOISE_MAX * h)
        return floor

    # --- views ----------------------------------------------------------------
    @property
    def metric_sensitivities(self):
        """Metric sensitivities as a long ``pandas.DataFrame``.

        One row per (metric, parameter) with the columns ``metric``,
        ``parameter``, ``metric_value`` (nominal metric ``m``),
        ``derivative`` (``dm/dp``), ``semi_elasticity`` (``p dm/dp``,
        defined where ``m`` is 0), ``elasticity``
        (``d ln m / d ln p = (p/m) dm/dp``; ``NaN`` where ``m`` is 0) and the
        parameter's finite-difference flags ``reliable``,
        ``event_count_changed`` and ``n_refinements`` (see :attr:`flags`):
        check them before using a row.

        The flags reflect event counts and failed runs only, not per-entry
        accuracy: the elasticity of a metric at the integrator's noise
        level -- e.g. ``Final('[S]')`` of a residual substrate below about
        ``1000 atol_y / rel_step`` (``atol_y``: the absolute tolerance in
        the metric's units, ``atol`` times the substrate's start value;
        1e-6 of the start value at the defaults ``atol = 1e-12``,
        ``rel_step = 1e-3``) -- is integration noise divided by a tiny
        ``m``, not a derivative, and is not masked here; use the
        ``derivative`` or ``semi_elasticity`` column there.
        """
        elasticity = _scale(self.metric_derivatives, self.metric_values,
                            self.values, 'log')
        semi = _scale(self.metric_derivatives, self.metric_values,
                      self.values, 'semi')
        n_m, n_p = len(self.metric_names), len(self.parameters)
        flags = self.flags
        return pd.DataFrame({
            'metric': pd.Series(np.repeat(self.metric_names, n_p),
                                dtype=object),
            'parameter': pd.Series(self.parameters * n_m, dtype=object),
            'metric_value': np.repeat(self.metric_values, n_p).astype(float),
            'derivative': self.metric_derivatives.reshape(-1).astype(float),
            'semi_elasticity': semi.reshape(-1).astype(float),
            'elasticity': elasticity.reshape(-1).astype(float),
            'reliable': np.tile(flags['reliable'].to_numpy(dtype=bool), n_m),
            'event_count_changed': np.tile(
                flags['event_count_changed'].to_numpy(dtype=bool), n_m),
            'n_refinements': np.tile(
                flags['n_refinements'].to_numpy(dtype=np.int64), n_m)})

    def _output_index(self, output):
        if output not in self.outputs:
            raise KeyError(f'{output!r} is not an output; outputs: '
                           f'{self.outputs}.')
        return self.outputs.index(output)

    def time_course(self, output):
        """Scaled sensitivities of one output over time.

        Parameters
        ----------
        output : str
            One of :attr:`outputs`.

        Returns
        -------
        pandas.DataFrame
            Index ``time`` (model time units), one column per parameter,
            values from :attr:`scaled` (log-scaled entries of an output at
            the integrator's noise level are ``NaN``; see :attr:`scaled`);
            ``.plot()`` draws them. The columns carry no reliability
            information: check :attr:`flags` (``reliable``,
            ``event_count_changed``) before trusting a column, or use
            :meth:`plot_time_course`, which draws unreliable parameters
            dashed.

        Raises
        ------
        KeyError
            If ``output`` is not an output.
        """
        k = self._output_index(output)
        return pd.DataFrame(self.scaled[:, k, :].copy(),
                            index=pd.Index(self.time.copy(), name='time'),
                            columns=pd.Index(self.parameters,
                                             name='parameter'))

    def _measure(self, name, at=None, arg='at'):
        """``(dy/dp per parameter, y)`` for an output (at ``at``, default
        the last time; linear interpolation in between) or a metric.
        ``arg`` names the time argument in error messages."""
        if name in self.metric_names and name not in self.outputs:
            if at is not None:
                raise ValueError(f'{arg}= applies to outputs, not to '
                                 'metrics.')
            k = self.metric_names.index(name)
            return self.metric_derivatives[k].copy(), self.metric_values[k]
        if name not in self.outputs:
            raise KeyError(f'{name!r} is neither an output nor a metric; '
                           f'outputs: {self.outputs}, metrics: '
                           f'{self.metric_names}.')
        k = self.outputs.index(name)
        if at is None:
            return self.raw[-1, k].copy(), self.nominal[-1, k]
        if isinstance(at, bool) or not isinstance(at, numbers.Real):
            raise TypeError(f'{arg} must be a time (number); got {at!r}.')
        t = float(at)
        if not math.isfinite(t):
            raise ValueError(f'{arg} must be a finite time; got {at!r}.')
        if not self.time[0] <= t <= self.time[-1]:
            raise ValueError(f'{arg}={at!r} is outside the output range '
                             f'[{self.time[0]!r}, {self.time[-1]!r}].')
        deriv = np.array([np.interp(t, self.time, self.raw[:, k, j])
                          for j in range(len(self.parameters))])
        return deriv, float(np.interp(t, self.time, self.nominal[:, k]))

    def rank(self, output_or_metric, at=None):
        """Parameters ranked by their effect on an output or a metric.

        Parameters
        ----------
        output_or_metric : str
            An output (its value at ``at``) or a metric name.
        at : float, optional
            Time (model time units) for an output; defaults to the last
            output time. Values between output times are interpolated
            linearly. Not allowed for a metric.

        Returns
        -------
        pandas.DataFrame
            Tornado-ready, one row per parameter sorted by decreasing
            magnitude of the scaled measure (``elasticity`` for
            ``scaling='log'``, ``semi_elasticity`` for ``'semi'``,
            ``derivative`` for ``'none'``; ``NaN`` last), with the columns
            ``parameter``, ``value`` (nominal ``p``), ``derivative``
            (``dy/dp``), ``semi_elasticity`` (``p dy/dp``), ``elasticity``
            (``(p/y) dy/dp``; ``NaN`` where ``y`` is 0, and for an output
            where ``|y|`` is at the integrator's noise level as in
            :attr:`scaled`) and the flags ``reliable``,
            ``event_count_changed``, ``n_refinements`` inline: read them
            before trusting a row. They reflect event counts and failed
            runs only: the elasticity of a metric at noise level (e.g. a
            residual substrate) is noise, not a derivative, and is not
            masked.

        Raises
        ------
        KeyError
            If the name is neither an output nor a metric.
        ValueError
            If ``at`` is not finite, outside the output range, or given for
            a metric.
        """
        deriv, y = self._measure(output_or_metric, at)
        with np.errstate(all='ignore'):
            semi = deriv * self.values
            elasticity = semi / y if y != 0 else np.full_like(semi, np.nan)
        if output_or_metric in self.outputs:
            floor = self._log_floor()[self.outputs.index(output_or_metric)]
            elasticity = np.where(abs(y) < floor, np.nan, elasticity)
        flags = self.flags
        table = pd.DataFrame({
            'parameter': pd.Series(self.parameters, dtype=object),
            'value': self.values.astype(float),
            'derivative': deriv.astype(float),
            'semi_elasticity': semi.astype(float),
            'elasticity': elasticity.astype(float),
            'reliable': flags['reliable'].to_numpy(dtype=bool),
            'event_count_changed':
                flags['event_count_changed'].to_numpy(dtype=bool),
            'n_refinements': flags['n_refinements'].to_numpy(dtype=np.int64)})
        key = np.abs(table[_RANK_KEY[self.scaling]].to_numpy())
        # Decreasing magnitude, ties in parameter order, NaN last.
        order = np.argsort(np.where(np.isnan(key), np.inf, -key),
                           kind='stable')
        return table.iloc[order].reset_index(drop=True)

    # --- plots ----------------------------------------------------------------
    def _ylabel(self, output):
        return {'log': f'd ln {output} / d ln p',
                'semi': f'd {output} / d ln p',
                'none': f'd {output} / dp'}[self.scaling]

    def plot_time_course(self, output, top=8, ax=None):
        """Plot the scaled sensitivities of one output over time.

        Parameters
        ----------
        output : str
            One of :attr:`outputs`.
        top : int
            Plot the ``top`` parameters with the largest peak magnitude
            (default 8).
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure otherwise.

        Returns
        -------
        (matplotlib.figure.Figure, matplotlib.axes.Axes)

        Raises
        ------
        KeyError
            If ``output`` is not an output.
        TypeError
            If ``top`` is not an integer.
        ValueError
            If ``top < 1``.

        Notes
        -----
        Parameters flagged unreliable are drawn dashed and labelled
        ``'(unreliable)'``.
        """
        top = _count('top', top, 1)
        course = self.time_course(output)
        values = np.abs(course.to_numpy())
        values[np.isnan(values)] = -np.inf
        peak = values.max(axis=0) if len(values) else np.zeros(
            len(self.parameters))
        order = np.argsort(-peak, kind='stable')[:top]
        _set_font()
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        reliable = self.flags['reliable'].to_numpy(dtype=bool)
        for j in order:
            label = self.parameters[j]
            style = '-'
            if not reliable[j]:
                label += ' (unreliable)'
                style = '--'
            ax.plot(self.time, course.iloc[:, j].to_numpy(), style,
                    label=label)
        ax.axhline(0.0, color='gray', linewidth=0.6)
        ax.set_xlabel(_bf('Time') + f" [{self.units.get('time', '')}]")
        ax.set_ylabel(self._ylabel(output))
        _style_axes(ax)
        ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left',
                  edgecolor='white')
        return fig, ax

    def plot_heatmap(self, t=None, ax=None):
        """Heat map of scaled sensitivities: outputs and metrics x parameters.

        Parameters
        ----------
        t : float, optional
            Time (model time units) of the output rows; defaults to the last
            output time. Metric rows are time-independent.
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure otherwise.

        Returns
        -------
        (matplotlib.figure.Figure, matplotlib.axes.Axes)

        Raises
        ------
        TypeError
            If ``t`` is not a number.
        ValueError
            If ``t`` is not finite or is outside the output range.

        Notes
        -----
        Cells use the result's scaling (``semi_elasticity`` / ``elasticity``
        / ``derivative`` as in :meth:`rank`), on a symmetric diverging colour
        scale; ``NaN`` cells are blank. Unreliable parameters are marked
        with ``*`` in the column labels.
        """
        rows, labels = [], []
        for name in self.outputs:
            rows.append(self._scaled_measure(name, t, arg='t'))
            labels.append(name)
        for name in self.metric_names:
            rows.append(self._scaled_measure(name, None))
            labels.append(name)
        matrix = np.array(rows, dtype=float).reshape(len(rows),
                                                     len(self.parameters))
        finite = np.abs(matrix[np.isfinite(matrix)])
        vmax = float(finite.max()) if finite.size and finite.max() > 0 else 1.0
        _set_font()
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        image = ax.imshow(np.ma.masked_invalid(matrix), cmap='RdBu_r',
                          vmin=-vmax, vmax=vmax, aspect='auto')
        reliable = self.flags['reliable'].to_numpy(dtype=bool)
        ax.set_xticks(range(len(self.parameters)))
        ax.set_xticklabels([p if ok else f'{p}*' for p, ok
                            in zip(self.parameters, reliable)], rotation=90)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels)
        when = self.time[-1] if t is None else float(t)
        ax.set_title(f"{_RANK_KEY[self.scaling].replace('_', ' ')} at t = "
                     f"{when:g} {self.units.get('time', '')}".rstrip())
        fig.colorbar(image, ax=ax)
        return fig, ax

    def _scaled_measure(self, name, at, arg='at'):
        deriv, y = self._measure(name, at, arg)
        floor = (self._log_floor()[self.outputs.index(name)]
                 if name in self.outputs else None)
        return _scale(deriv, y, self.values, self.scaling, floor)

    # --- IO -------------------------------------------------------------------
    def to_csv(self, path):
        """Write the result to a long-format CSV.

        The file has the columns ``quantity, key, parameter, index, value``.
        ``quantity`` is ``meta`` (``key``: the field name), ``time``
        (``index``: row), ``nominal`` (``key``: output, ``index``: row),
        ``raw`` (``key``: output, ``parameter``, ``index``: row),
        ``parameter_value``, ``flag`` (``key``: flag column),
        ``metric_value`` (``key``: metric) or ``metric_derivative``
        (``key``: metric, ``parameter``). Every ``value`` is JSON (floats in
        their shortest round-trip form, ``NaN``/``Infinity`` for non-finite
        ones), so :meth:`from_csv` restores the result exactly; the scaled
        arrays are recomputed from ``raw`` by the same function.

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
        meta = {'nskinetics': _CSV_FORMAT, 'outputs': self.outputs,
                'parameters': self.parameters,
                'metric_names': self.metric_names, 'scaling': self.scaling,
                'method': self.method, 'rel_step': self.rel_step,
                'counters': list(self.counters), 'units': self.units,
                'condition': (None if self.condition is None
                              else self.condition.to_dict()),
                'n_time': len(self.time),
                'flag_columns': list(self.flags.columns),
                'tolerances': (None if self.tolerances is None
                               else list(self.tolerances)),
                'max_refinements': self.max_refinements}
        rows = [('meta', key, '', '', json.dumps(value, allow_nan=False))
                for key, value in meta.items()]
        dumps = json.dumps
        rows += [('time', '', '', i, dumps(float(v)))
                 for i, v in enumerate(self.time)]
        for k, out in enumerate(self.outputs):
            rows += [('nominal', out, '', i, dumps(float(v)))
                     for i, v in enumerate(self.nominal[:, k])]
        for k, out in enumerate(self.outputs):
            for j, par in enumerate(self.parameters):
                rows += [('raw', out, par, i, dumps(float(v)))
                         for i, v in enumerate(self.raw[:, k, j])]
        for j, par in enumerate(self.parameters):
            rows.append(('parameter_value', '', par, '',
                         dumps(float(self.values[j]))))
            for col in self.flags.columns:
                value = self.flags[col].iloc[j]
                if isinstance(value, (np.bool_, bool)):
                    value = bool(value)
                elif isinstance(value, np.integer):
                    value = int(value)
                elif isinstance(value, np.floating):
                    value = float(value)
                rows.append(('flag', col, par, '', dumps(value)))
        for m, name in enumerate(self.metric_names):
            rows.append(('metric_value', name, '', '',
                         dumps(float(self.metric_values[m]))))
            for j, par in enumerate(self.parameters):
                rows.append(('metric_derivative', name, par, '',
                             dumps(float(self.metric_derivatives[m, j]))))
        pd.DataFrame(rows, columns=list(_CSV_COLUMNS)).to_csv(
            path, index=False, lineterminator='\n', encoding='utf-8')

    @classmethod
    def from_csv(cls, path):
        """Read a result written by :meth:`to_csv`.

        Parameters
        ----------
        path : str or path-like
            CSV written by :meth:`to_csv`.

        Returns
        -------
        LocalSensitivityResult
            Equal to the one written: same arrays (exact values, ``NaN`` in
            place), flags table (columns, dtypes, values), names, scaling,
            settings (method, steps, counters, tolerances,
            ``max_refinements``), units and condition.

        Raises
        ------
        ValueError
            If the file was not written by :meth:`to_csv`, or is incomplete
            (a missing ``meta`` key, or an entry of an array or of the flags
            missing or duplicated, e.g. a truncated or edited file).
        """
        meta, groups = _read_long_csv(path, _CSV_FORMAT, required=(
            'outputs', 'parameters', 'metric_names', 'scaling', 'method',
            'rel_step', 'counters', 'units', 'condition', 'n_time',
            'flag_columns', 'tolerances', 'max_refinements'))
        outputs, params = meta['outputs'], meta['parameters']
        metrics, n_t = meta['metric_names'], meta['n_time']
        out_idx = {o: k for k, o in enumerate(outputs)}
        par_idx = {p: j for j, p in enumerate(params)}
        met_idx = {m: k for k, m in enumerate(metrics)}
        n_o, n_p, n_m = len(outputs), len(params), len(metrics)

        def scatter(quantity, shape, *axes, **kwargs):
            return _scatter(path, _CSV_FORMAT, groups, quantity, shape, axes,
                            **kwargs)

        time = scatter('time', (n_t,), ('index', None))
        nominal = scatter('nominal', (n_t, n_o), ('index', None),
                          ('key', out_idx))
        raw = scatter('raw', (n_t, n_o, n_p), ('index', None),
                      ('key', out_idx), ('parameter', par_idx))
        values = scatter('parameter_value', (n_p,), ('parameter', par_idx))
        metric_values = scatter('metric_value', (n_m,), ('key', met_idx))
        metric_derivatives = scatter('metric_derivative', (n_m, n_p),
                                     ('key', met_idx),
                                     ('parameter', par_idx))
        columns = meta['flag_columns']
        missing = [c for c in _fd.FLAG_COLUMNS if c not in columns]
        if missing:
            raise _incomplete(path, _CSV_FORMAT,
                              f'flag_columns lacks {missing}')
        cells = scatter('flag', (len(columns), n_p),
                        ('key', {c: i for i, c in enumerate(columns)}),
                        ('parameter', par_idx), parse=lambda v: v,
                        dtype=object)
        flag_rows = [{c: cells[i, j] for i, c in enumerate(columns)}
                     for j in range(n_p)]
        flags = _fd.flags_frame(flag_rows, params)[columns]
        condition = (None if meta['condition'] is None
                     else Condition.from_dict(meta['condition']))
        return cls(time=time, outputs=outputs, parameters=params,
                   values=values, nominal=nominal, raw=raw,
                   metric_names=metrics, metric_values=metric_values,
                   metric_derivatives=metric_derivatives, flags=flags,
                   scaling=meta['scaling'], method=meta['method'],
                   rel_step=meta['rel_step'],
                   counters=tuple(meta['counters']), units=meta['units'],
                   condition=condition, tolerances=meta['tolerances'],
                   max_refinements=meta['max_refinements'])


# --- public API -----------------------------------------------------------------
def _sensitivity_vector(results, values, outputs, metrics):
    """Stack the output trajectories (row-major, ``[t, output]``) and the
    metric values of one run into one vector (the function differenced)."""
    result = results[0]
    trajectories = (result[outputs] if outputs
                    else np.empty((len(result), 0)))
    values = []
    for name, metric in metrics.items():
        with np.errstate(all='ignore'):
            values.append(_ex._scalar(name, metric(result)))
    return np.concatenate([trajectories.reshape(-1),
                           np.asarray(values, dtype=float)])


def local_sensitivity(model, condition, params, *, outputs=None, metrics=None,
                      method='central', rel_step=1e-3, scaled='log',
                      event_counters=None, workers=1, abs_step=None,
                      tolerances=None, max_refinements=2):
    """Event-aware local sensitivities of trajectories and metrics.

    Each parameter is perturbed by ``p * exp(+-rel_step)`` (central
    differences in ``ln p``; ``forward`` uses the ``+`` side only) and the
    condition is re-run from the restored baseline on a private copy of the
    model. Perturbations are merged into ``condition`` with
    :meth:`~nskinetics.Condition.with_parameters` /
    :meth:`~nskinetics.Condition.with_initial`, so
    its own values stay in force. Event counters are watched: if a
    perturbation changes a counter at any output time (a firing added,
    removed, or moved across an output time), the step is refined (divided
    by 10, up to ``max_refinements`` times) and a parameter whose counters
    still differ is flagged (see :attr:`LocalSensitivityResult.flags`).

    Parameters
    ----------
    model : KineticModel
        The model; never modified (all runs use a private copy).
    condition : Condition
        The nominal run (``reset=True``). Times are model time units.
    params : sequence of str
        Ids to perturb: global parameters (e.g. ``'k'``), initial values
        (``'[S]'`` concentration, ``'S'`` amount, ``'init(x)'``),
        compartments (``'env'``)... A value the condition sets is the
        nominal value; any other is the value the condition's run starts
        from. Every run, the nominal one included, sets each id explicitly
        (see Notes).
    outputs : sequence of str, optional
        Selections whose time courses are differentiated. Defaults to the
        condition's selections (or the model's floating species
        concentrations) without ``'time'``. May be empty if ``metrics`` are
        given.
    metrics : mapping or sequence, optional
        ``{name: metric}`` or metric objects with a ``name`` (see
        :class:`Metric`); their sensitivities are in
        :attr:`LocalSensitivityResult.metric_sensitivities`.
    method : {'central', 'forward'}
        Finite-difference scheme (default ``'central'``): second-order
        accurate (truncation error ~``h^2``, relative), or first-order
        (~``h / 2``) at half the runs. Use ``'central'`` when the
        sensitivities feed :func:`identifiability`, whose rank test needs
        derivatives accurate to well below ``1e-6`` (see its Notes).
    rel_step : float
        Step in ``ln p`` (in ``(0, 1]``); default ``1e-3``. A larger step
        raises the truncation error, which matters for
        :func:`identifiability`; keep it at most ``1e-3`` there.
    scaled : {'log', 'semi', 'none'}
        Scaling of :attr:`LocalSensitivityResult.scaled` and of
        :meth:`~LocalSensitivityResult.rank` (default ``'log'``):
        ``d ln y / d ln p``, ``dy / d ln p`` (use where ``y`` passes
        through 0) or ``dy/dp``.
    event_counters : sequence of str, optional
        Counters to watch. ``None`` (default) detects every variable that
        only events change and that an event sets to ``v + 1`` (FeedSpike
        ``count_var`` and user counters; a reacting species pulsed by
        ``S = S + 1`` is not one, and a concentration species in a
        compartment whose size changes is left out with a warning, since
        its amount is not a count -- count with a parameter, or list it
        here; one in a constant compartment is watched through its
        concentration ``'[n]'``, the form to pass here too); ``()``
        disables the check.
    workers : int
        Processes for the runs (default ``1``: serial; see
        :func:`nskinetics.run_many`).
    abs_step : float or mapping, optional
        Absolute steps (model units) for parameters differenced in linear
        space: those whose nominal value is ``<= 0`` (scalar or mapping
        entry), and any parameter listed in a mapping. Default for a linear
        parameter: ``rel_step * |p|``, or ``rel_step`` for ``p = 0``. That
        last default is in the parameter's own units and unrelated to the
        model's scale (``1e-3`` mol/L is far too large a step for a
        micromolar product), so a warning is issued for each parameter at
        0 without an ``abs_step``.
    tolerances : (float, float), optional
        ``(atol, rtol)`` for the private copy's integrator, used as given.
        ``None`` (default) tightens the copy to at most ``(1e-12, 1e-10)``
        (tighter model values are kept), and then ``atol`` further, by the
        smallest peak magnitude below 1 of the ``outputs``, and of the
        columns the ``metrics`` read, that start at 0 (one extra run):
        RoadRunner gives such a state a fixed absolute
        tolerance ``atol`` for the whole run, too coarse for a product at
        the micro- or nanomolar scale. The tolerances used are recorded in
        :attr:`LocalSensitivityResult.tolerances`. Tighter tolerances need
        more integrator steps, so the copy's ``maximum_num_steps`` is then
        raised 10-fold; if a run still fails where ``model.run(condition)``
        succeeds, pass looser ``tolerances=``. A ``UserWarning`` is issued
        if the relative tolerance in use exceeds ``rel_step / 100``. That
        check catches gross mismatches only: below it, errors of a few
        percent remain possible (``rtol = 1e-6`` with ``rel_step = 1e-3``
        gives about 1% on a plain decay, an ``rtol`` just under
        ``rel_step / 100`` about 3%), so keep ``rtol`` near the default
        ``1e-10`` when accuracy matters.
    max_refinements : int
        Maximum number of step refinements on a count change (default 2).

    Returns
    -------
    LocalSensitivityResult

    Raises
    ------
    TypeError, ValueError
        For invalid arguments, ids that do not fit the model, or a metric
        that cannot be evaluated on the condition (e.g. :class:`At` with
        ``t`` outside its time span); the message names the argument.
        Raised before any run (the caller's model is never touched).
    KineticSimulationError
        If the nominal run (or the run that reads start values) fails; the
        message names the integrator tolerances in use and the remedies. A
        failed perturbed run makes only that parameter's column ``NaN`` and
        unreliable.

    Warns
    -----
    UserWarning
        If the integrator's relative tolerance exceeds ``rel_step / 100``,
        or ``1/100`` of the relative size of a smaller step a reported
        difference used (a refined step, or one set by ``abs_step``; gross
        mismatches only, see ``tolerances``); if a perturbed id is read by
        an SBML initial assignment (see Notes); if a counter candidate is
        left out (see ``event_counters``); if a parameter at 0 has no
        ``abs_step`` (its default step, ``rel_step`` in its own units, is
        unrelated to the model's scale); if one side of a central
        difference failed to run and a one-sided difference was used.

    Notes
    -----
    Sensitivities depend on the condition (start values, feeding
    parameters, output grid): they are local. Every listed id is held at
    its nominal value in every run, so each column is a partial derivative
    with the other listed ids fixed, and a column can change when another
    id is added. Changing a compartment keeps the species amounts unless a
    concentration ``'[S]'`` is also listed, in which case the concentration
    is held: on a toy with rate ``k*S*env``, ``d[P]/d env`` is -9.50 with
    ``params=['env']`` and 0 with ``['env', '[S]']``.

    A listed id is set after the reset has evaluated the model's initial
    assignments, so they are not re-evaluated: for ``S = S_init`` (an SBML
    initial assignment), perturbing ``'S_init'`` leaves ``S`` at its
    nominal start and gives ``dP/dS_init = 0`` (a warning says so).
    Perturb ``'init(S_init)'`` instead, which re-evaluates the assignment.

    Read the flags before using a value: ``event_count_changed=True`` means
    that even at the smallest step the perturbation moved a counted event
    across an output time or the end of the run (or added or removed one),
    where the derivative does not exist (on the affected rows the
    difference would measure the jump). The firings each counter counts by
    every output time (counter minus its start value, so perturbing a
    counter's own start value is no change) are compared with the nominal
    run's; a crossing is refined first, and
    flagged only if it persists after ``max_refinements``, i.e. if the
    output time lies closer to the firing time than the shift that the
    smallest step (``rel_step / 10**max_refinements``) causes -- about
    1e-5 of the firing time when it scales with ``1/p``, as for a rate
    constant. The reported column is then the one-sided
    difference on the side whose counters agree with the nominal run, if
    exactly one does, at the largest step at which it agreed
    (``flags['step']``). Stage switches that do not increment a counter are
    not detected automatically; list such variables in
    ``event_counters``. A parameter at 0 is differenced in linear space
    (central differences then also run it at ``-step``); its log
    elasticity is 0. Choose its ``abs_step`` on the quantity's scale,
    e.g. ``1e-3`` of the concentrations it competes with: the default,
    ``rel_step`` in its own units, can be orders of magnitude too large
    (wrong derivatives, even of the wrong sign, still flagged reliable).
    If one side of a central difference fails to run
    (e.g. a negative value that breaks the integration) while the other
    runs with the nominal firings, that side's one-sided difference is
    used (``flags['scheme']`` is ``'forward'`` or ``'backward'``, the
    column stays reliable, and a warning names the parameter); if the
    running side counted other firings, the step is refined first, as for
    any count change. The column is ``NaN`` only if both sides fail, or if
    the running side never agrees with the nominal firings (then
    ``event_count_changed`` is set).

    ``reliable`` reflects event counts and failed runs, not the accuracy of
    each entry. A log-scaled sensitivity (or elasticity) of a quantity at
    the integrator's noise level is noise, not a derivative: its
    difference is integration error of about ``3 atol_y / h`` divided by
    a tiny ``y``, where ``atol_y`` is the absolute tolerance in the
    output's units. RoadRunner multiplies a scalar ``absolute_tolerance``
    by each state's start value, so ``atol_y = atol |y0|`` for an output
    starting at ``y0 != 0`` (and ``atol`` for one starting at 0, which
    the default ``tolerances`` scale to the smallest peak of the outputs
    and metric selections that start at 0): the noise level is relative
    to the start value, or to the smallest such peak. Log-scaled output
    entries below ``|y| = 30 atol_y / h`` (``|y| / |y0|`` below 3e-8 at
    the defaults ``atol = 1e-12``, ``rel_step = 1e-3``; ``|y|`` below
    ``3e4`` times the ``atol`` in use, recorded in
    :attr:`LocalSensitivityResult.tolerances`, for an output starting at
    0) are ``NaN`` (see :attr:`LocalSensitivityResult.scaled`); entries
    up to about ``1000 atol_y / h`` (1e-6 relative), and metric
    elasticities of a near-zero metric such as ``Final('[S]')`` of a
    residual substrate, are kept but noisy: use ``scaled='semi'`` or
    ``'none'``, or leave those times out.

    Metrics are differenced as computed on the output grid, so a metric
    that is piecewise constant in the parameters -- :class:`EventCount`, or
    a :class:`CallableMetric` that returns a grid time -- gives 0 or a jump
    rather than a derivative; use :func:`morris` / :func:`sobol` for such
    metrics, or a smooth basis (:class:`TimeAtMax` refines the time of the
    maximum between grid points for this reason).

    Integration error enters the differences divided by the step. At
    CVODE's default ``rtol = 1e-6`` and ``rel_step = 1e-3`` it reaches ~1%
    of the derivative even on a plain decay (``S => P; (k1 + k2)*S``),
    hence the tightened default, which keeps it below 1e-6 at a modest
    cost; refined steps (``1e-5`` after two refinements) amplify it
    100-fold.

    Cost: two runs per parameter (one for ``method='forward'``), the
    nominal run, one more run to read start values the condition does not
    set, one to scale ``atol`` (default ``tolerances`` only), and two more
    per step refinement. Each call also first builds a private copy of the
    model (:meth:`~nskinetics.KineticModel.copy`), which takes about a
    second for the shipped yeast model, far longer than one of its runs:
    for a few parameters that rebuild dominates, so a loop of calls (e.g.
    sensitivities at many operating points) pays mostly for rebuilding.

    Examples
    --------
    >>> ls = local_sensitivity(km, cond, ['k', 'threshold'],     # doctest: +SKIP
    ...                        outputs=['[P]'], metrics={'P': Final('[P]')})
    >>> ls.rank('P')                                              # doctest: +SKIP
    """
    if not isinstance(model, KineticModel):
        raise TypeError(f'model must be a KineticModel; got '
                        f'{type(model).__name__}.')
    if not isinstance(condition, Condition):
        raise TypeError(f'condition must be a Condition; got '
                        f'{type(condition).__name__}.')
    params = _id_list('params', params)
    if method not in _fd.METHODS:
        raise ValueError(f'method must be one of {_fd.METHODS}; got '
                         f'{method!r}.')
    if scaled not in _SCALINGS:
        raise ValueError(f'scaled must be one of {_SCALINGS}; got '
                         f'{scaled!r}.')
    rel_step = _fd.check_step(rel_step)
    abs_step = _fd.check_abs_step(abs_step, params)
    _fd.check_tolerances(tolerances)
    max_refinements = _count('max_refinements', max_refinements, 0)
    counters = (None if event_counters is None
                else _fd.check_counters(event_counters))
    # Ids, labelled by the argument they came from; then batch validation on
    # the caller's model (touches no state): reset semantics, selections.
    _check_ids(model, 'params', params)
    _check_ids(model, 'metrics',
               _ex._metric_selections(_ex._normalize_metrics(metrics)),
               settable=False)
    _, metrics, metric_sels, workers = _ex._validate_batch(
        model, [condition], metrics, workers, True, 'record')
    # e.g. At(t) outside the condition's time span: before any run
    _check_metrics_condition(metrics, condition)
    if outputs is None:
        outputs = [s for s in model._run_selections(condition) if s != 'time']
    outputs = _id_list('outputs', outputs, allow_empty=True)
    if 'time' in outputs:
        raise ValueError("outputs cannot include 'time'.")
    if not outputs and not metrics:
        raise ValueError('Nothing to differentiate: give outputs or metrics.')
    clash = [m for m in metrics if m in outputs]
    if clash:
        raise ValueError(f'Metric name(s) {clash} clash with output names.')
    _check_ids(model, 'outputs', outputs, settable=False)
    if counters:
        _check_ids(model, 'event_counters', counters, settable=False)
    # Each perturbed condition must be constructible (e.g. '[S]' with an 'S'
    # amount already in the condition is not); the values are irrelevant.
    try:
        _fd.pin(model, condition, dict.fromkeys(params, 1.0))
    except ValueError as e:
        raise ValueError(
            f'params: {e} List each species once, in the same form (amount '
            "id or '[concentration]') as the condition's own initial "
            'values.') from e
    _fd.warn_initial_assignment_inputs(model, params, stacklevel=3)

    private = _fd.private_model(model, tolerances=tolerances,
                                rel_step=rel_step)
    # Every column recorded: the outputs and the selections the metrics read.
    extra = tuple(dict.fromkeys([*outputs, *metric_sels]))
    if tolerances is None:
        # RoadRunner integrates a state that starts at 0 to a fixed absolute
        # tolerance: scale atol to the peak of such outputs and metric
        # selections (one run), so their differences -- and the recorded
        # noise model -- hold at any concentration scale, and a metric's
        # derivative does not depend on whether its column is also an
        # output. Explicit tolerances are kept as given.
        _fd.scale_atol_to_outputs(private, [condition], extra)
    if counters is None:
        counters = private._fd_counters         # detected once, on the copy
        _fd.warn_unwatched_counters(private._fd_unwatched, stacklevel=3)
    values = _fd.nominal_values(private, condition, params)
    # A parameter at 0 without an abs_step gets rel_step in its own units,
    # whatever the model's scale (_fd stays silent: it runs in optimizers).
    defaulted = [p for p, v in values.items() if v == 0
                 and not (abs_step is not None
                          and (not isinstance(abs_step, Mapping)
                               or p in abs_step))]
    if defaulted:
        warnings.warn(
            f'{defaulted} are at 0, so they are differenced in linear space '
            f'with the default absolute step rel_step={rel_step:g} in their '
            "own units, which is unrelated to the model's scale: a step far "
            "above the quantity's typical size (e.g. 1e-3 mol/L against "
            'micromolar concentrations) gives wrong derivatives, and '
            'reliable=True does not check this. Pass abs_step={id: step} on '
            "the quantity's scale.", UserWarning, stacklevel=2)
    fd_params = _fd.make_parameters(values, rel_step=rel_step,
                                    abs_step=abs_step)
    fn = lambda results, values: _sensitivity_vector(  # noqa: E731
        results, values, outputs, metrics)
    jac = _fd.fd_jacobian(private, [condition], fd_params, fn, method=method,
                          event_counters=counters,
                          max_refinements=max_refinements, workers=workers,
                          extra_selections=extra)
    # A smaller relative step than rel_step (which private_model checked) --
    # a refined step, or a small abs_step -- amplifies integration error:
    # warn for the differences whose relative step fell below 100 rtol.
    noisy = _fd.noisy_parameters(jac.flags, jac.tolerances, checked=rel_step)
    if noisy:
        warnings.warn(
            f'The finite differences of {noisy} use a step (refined, or set '
            'by abs_step) whose relative size is below 100 x the integrator '
            f'relative tolerance ({jac.tolerances[1]:g}), so integration '
            'error may dominate them; pass tighter tolerances=(atol, rtol)'
            ' or a larger abs_step.', UserWarning, stacklevel=2)
    # Central differences that fell back to one side because the other
    # side's run failed (the engine does not warn; it runs in optimizers).
    one_sided = [name for name, row in jac.flags.iterrows()
                 if method == 'central' and row['scheme'] != 'central'
                 and not row['event_count_changed']]
    if one_sided:
        warnings.warn(
            f'One side of the central difference of {one_sided} failed to '
            "run; a one-sided difference was used (flags['scheme']).",
            UserWarning, stacklevel=2)
    result = jac.results[0]
    n_t, n_out = len(result), len(outputs)
    natural = jac.natural
    n_traj = n_t * n_out
    return LocalSensitivityResult(
        time=result.time, outputs=outputs, parameters=params,
        values=np.array([p.value for p in fd_params], dtype=float),
        nominal=jac.f0[:n_traj].reshape(n_t, n_out),
        raw=natural[:n_traj].reshape(n_t, n_out, len(params)),
        metric_names=list(metrics), metric_values=jac.f0[n_traj:].copy(),
        metric_derivatives=natural[n_traj:].copy(), flags=jac.flags,
        scaling=scaled, method=method, rel_step=rel_step,
        counters=jac.counters, units=dict(model._units),
        condition=condition, tolerances=jac.tolerances,
        max_refinements=max_refinements)
