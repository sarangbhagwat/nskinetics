# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Event-aware finite differences on a private model (private engine).

The one finite-difference (FD) engine behind
:func:`~nskinetics.analysis.local_sensitivity`, identifiability (the Fisher
information matrix) and calibration (the weighted-residual Jacobian, the
covariance of a fit and profile likelihoods).

Why finite differences: RoadRunner's forward sensitivities drop the jumps
that events make in the state, so derivatives across a feed spike come out
wrong (and all zero on the shipped model). Here every evaluation is a full
:meth:`KineticModel.run` from the restored baseline, so event jumps are
part of the trajectory being differenced. What FD cannot handle is a
perturbation that changes *how many times* an event fires -- the output then
jumps and ``(f+ - f-) / 2h`` is meaningless. The same holds at any output
time that a perturbation moves a firing across: on that row the difference
measures the jump. Event counters (variables only events change,
incremented by exactly one, ``v = v + 1``: a FeedSpike's ``count_var`` or a
user counter) are therefore recorded on every run, and the firings they
count by each output time (counter minus its start value) are compared with
the nominal run *at every output time*; when a perturbed run counted
different firings by some output time, the step is divided by 10 (up to
``max_refinements`` times) and, if they still differ, the parameter is
flagged ``event_count_changed=True, reliable=False`` and the one-sided
difference on the side whose firings agree with the nominal run is kept
(when exactly one side does), taken at the largest step at which that side
agreed (integration error grows as the step shrinks). A central difference
one of whose sides failed to run also falls back to the other side's
one-sided difference.

Integration error enters every difference divided by the step, so the
private copy's integrator is tightened by default (:data:`FD_TOLERANCES`),
and its step budget raised (:data:`STEP_BUDGET_FACTOR`); its absolute
tolerance can further be scaled to the outputs that start at 0
(:func:`scale_atol_to_outputs`).

Layout:

- :class:`FDParameter`: one perturbed quantity -- a global parameter id or an
  initial-value selection such as ``'[S]'`` -- with its nominal value, the
  space it is differenced in (``ln p`` or ``p``), the step in that space and
  the conditions it applies to (all, or a subset, e.g. one experiment's local
  parameter in calibration).
- :func:`fd_jacobian`: the Jacobian of an arbitrary vector function of the
  results of several :class:`~nskinetics.Condition` objects. Perturbations
  are merged into each condition with :meth:`Condition.with_parameters` /
  :meth:`Condition.with_initial` (placement as in :func:`nskinetics.scan`),
  so the condition's own values survive. Runs go through F1's executor
  (serially on the given model, or one model per worker process).
- :func:`evaluate`: the same runs and ``fn`` at the nominal point only (an
  :class:`FDPoint`), recording failures instead of raising -- calibration's
  residual function -- which :func:`fd_jacobian` accepts as ``nominal=`` so
  the nominal runs are not repeated.
- :func:`private_model`, :func:`scale_atol_to_outputs`,
  :func:`nominal_values`, :func:`make_parameters`,
  :func:`detect_event_counters`: the set-up helpers.
- :func:`noisy_parameters`, :func:`effective_atol`: the parameters whose
  step is too small for the tolerances, and the absolute tolerance in each
  output's units.
- :func:`error_level`, :func:`truncation_error`: the estimated error of the
  differences relative to the derivatives (integration noise and
  truncation), the accuracy of a rank test on a Jacobian.
"""
import functools
import math
import numbers
import warnings
from collections.abc import Mapping
from dataclasses import dataclass
from typing import NamedTuple

import libsbml
import numpy as np
import pandas as pd

from ..engine.conditions import Condition
from ..engine import execution as _ex
from ..exceptions import KineticSimulationError
from ._common import _count, _finite

__all__ = ()

METHODS = ('central', 'forward')
#: ``(atol, rtol)`` the private copy's integrator is tightened to (each at
#: most) when no ``tolerances`` are given. With CVODE's defaults (rtol
#: 1e-6) the global error of a run is ~1e-5 relative, which a step of 1e-3
#: turns into ~1% errors in the derivatives; at rtol 1e-10, below 1e-6.
FD_TOLERANCES = (1e-12, 1e-10)
#: Smallest factor by which :func:`scale_atol_to_outputs` tightens the
#: absolute tolerance to the peak of an output that starts at 0.
ATOL_MIN_FACTOR = 1e-15
#: Largest ``atol`` the finite-difference error level (:func:`error_level`)
#: is valid for: it models ``rtol`` only, assuming RoadRunner's
#: start-value-scaled ``atol`` is at most this (the default
#: :data:`FD_TOLERANCES` atol is 1e-12, scaled to the outputs). On
#: ``S => P; (k1 + k2)*S`` with ``S0 = 10``, ``rtol = 1e-10``,
#: ``sigma=('relative', 0.1)`` and the run ended at ``t_end = 10`` (61
#: output times), an ``atol`` up to 3e-7 keeps rank 1, while 1e-6 gives
#: rank 2 with se ~2e3 (~20 at ``rel = 1e-3``: se scales with ``rel``).
#: The effect is erratic in ``atol`` and run length (1e-5 keeps rank 1 at
#: ``t_end = 10`` but not at 15, and 3e-7 gives rank 2 at 30); 1e-7 kept
#: rank 1 at every run length tried (5, 10, 15 and 20). identifiability
#: warns above it, and records no error level.
ATOL_MAX = 1e-8
#: Factor by which the private copy's ``maximum_num_steps`` (CVODE's step
#: budget per output interval) is raised when its relative tolerance is
#: tightened below the model's own: tighter tolerances need more internal
#: steps, and a larger budget never changes a run that already succeeds.
STEP_BUDGET_FACTOR = 10
#: A difference at relative step ``h`` is flagged as noise-dominated when the
#: integrator's ``rtol`` exceeds ``h / TOLERANCE_RATIO``. This catches gross
#: mismatches only: just below the threshold, errors of a few percent remain
#: possible (e.g. ~1% at ``rtol = 1e-6`` with ``h = 1e-3`` on a plain decay,
#: ~3% at ``rtol`` just under ``h / 100``); hence the tighter default.
TOLERANCE_RATIO = 100
#: Finite-difference noise per unit of integration error: a difference at
#: step ``h`` (in the differencing space) of an output ``y`` is taken to
#: carry noise up to about ``FD_NOISE * (atol_y + rtol * |y|) / h``, with
#: ``atol_y`` the absolute tolerance in the output's units
#: (:func:`effective_atol`: RoadRunner scales a scalar ``atol`` by each
#: state's start value). Measured 0.1-3 x ``rtol / h`` on toys at rtol
#: 1e-10; near a depleted species the ``atol`` term dominates, and the
#: error of ``d ln S / d ln k`` on ``S => P; k*S`` depends on ``S / S0``,
#: not on ``S`` (the same at ``S0 = 10`` and ``1e-5``: about 0.02-0.1 at
#: ``S / S0 = 8e-10`` and 0.6-5 at ``4e-11``, with ``atol = 1e-12`` and
#: ``h = 1e-3``; 4 and 80 estimated, so the estimate is conservative).
#: Used for the log-scaling mask of local sensitivities and the
#: finite-difference error level (:func:`error_level`) behind the rank
#: warning of identifiability.
FD_NOISE = 3.0
#: Truncation error of a difference at relative step ``h``, relative to the
#: derivative: about ``TRUNCATION_ONE_SIDED * h`` for a one-sided difference
#: (first order; the leading Taylor coefficient) and ``TRUNCATION_CENTRAL *
#: h**2`` for a central one. Unlike integration noise it does not cancel in
#: an exact non-identifiability that is not linear in ``ln p``: on ``S => P;
#: (k1 + k2)*S`` it lifts the null singular value to 0.03-0.1 ``h`` of the
#: largest (forward) and 0.01-0.26 ``h**2`` (central, the largest value
#: measured, under a relative noise model). A product such as ``k1*k2``
#: cancels at any step. Calibrated on smooth toys, where the higher
#: derivatives in ``ln p`` are of the order of the first; at a depletion
#: front or a switch a step in ``ln p`` shifts the front in time, and the
#: truncation error is orders of magnitude larger (630 x the estimate on
#: the ``m1``/``m2`` null of Monod ``(m1 + m2) * S / (Ks + S) * X`` at
#: ``h = 1e-3``). Used by :func:`error_level`.
TRUNCATION_ONE_SIDED = 0.5
TRUNCATION_CENTRAL = 0.3
#: Columns of the per-parameter flags table, in order, with their dtypes.
FLAG_COLUMNS = {'value': float, 'log': bool, 'step': float,
                'n_refinements': np.int64, 'event_count_changed': bool,
                'reliable': bool, 'scheme': object, 'error': object}


# --- parameters -----------------------------------------------------------------
@dataclass(frozen=True)
class FDParameter:
    """One perturbed quantity of a finite-difference Jacobian.

    Attributes
    ----------
    name : str
        Column name of the Jacobian (e.g. ``'k'`` or ``'[S]@FB1'``).
    target : str or None
        What is set on the model: a global parameter id, a compartment or an
        initial-value selection such as ``'[S]'`` (placed into a condition's
        ``parameters`` or ``initial`` like a :func:`nskinetics.scan` key).
        ``None`` for a quantity that only ``fn`` uses (e.g. an observable
        scale in calibration): perturbing it re-evaluates ``fn`` on the
        nominal results with the new value in ``values``, without any run.
    value : float
        Nominal value, in the model's own units. Every run of an applicable
        condition -- the nominal one included -- sets ``target`` explicitly.
    log : bool
        Difference in ``x = ln(value)`` (requires ``value > 0``) instead of
        ``x = value``. Defaults to ``True``.
    step : float
        Step ``h`` in that space: a relative step for ``log=True``
        (``value * exp(+-h)``), an absolute one (model units) otherwise.
    conditions : tuple of int or None
        Indices of the conditions the quantity applies to; ``None`` (default)
        means all. Ignored when ``target`` is ``None``.
    """
    name: str
    target: str | None
    value: float
    log: bool = True
    step: float = 1e-3
    conditions: tuple | None = None

    def __post_init__(self):
        for arg in ('name', 'target'):
            v = getattr(self, arg)
            if (arg == 'target' and v is None) or (isinstance(v, str) and v):
                continue
            raise TypeError(f'FDParameter.{arg} must be a non-empty string'
                            f"{' or None' if arg == 'target' else ''}; got "
                            f'{v!r}.')
        put = functools.partial(object.__setattr__, self)
        put('value', _finite(f'value of {self.name!r}', self.value))
        put('step', _finite(f'step of {self.name!r}', self.step))
        if not isinstance(self.log, (bool, np.bool_)):
            raise TypeError(f'log of {self.name!r} must be a bool.')
        put('log', bool(self.log))
        if self.step <= 0:
            raise ValueError(f'step of {self.name!r} must be > 0; got '
                             f'{self.step!r}.')
        if self.log and self.value <= 0:
            raise ValueError(f'{self.name!r} = {self.value!r} cannot be '
                             'differenced in log space; use log=False with an '
                             'absolute step.')
        if self.conditions is not None:
            idx = tuple(sorted(set(self.conditions)))
            if not idx or any(isinstance(i, bool)
                              or not isinstance(i, numbers.Integral) or i < 0
                              for i in idx):
                raise ValueError(f'conditions of {self.name!r} must be a '
                                 'non-empty sequence of condition indices.')
            put('conditions', tuple(int(i) for i in idx))

    def applies_to(self, i):
        """Whether the quantity is set in condition ``i`` (never for a
        ``target`` of ``None``)."""
        if self.target is None:
            return False
        return self.conditions is None or i in self.conditions

    def at(self, offset):
        """Natural value at ``offset`` in the differencing space."""
        if self.log:
            return self.value * math.exp(offset)
        return self.value + offset


def relative_step(value, log, step):
    """The step of a difference relative to the parameter value: ``step``
    itself in log space, ``step / |value|`` in linear space (``NaN`` at
    ``value = 0``, where no relative step exists)."""
    if log:
        return float(step)
    return float(step) / abs(value) if value != 0 else math.nan


def truncation_error(flags):
    """Estimated truncation error of the differences in ``flags`` (an FD
    flags table), relative to the derivatives.

    Returns ``(one_sided, central, names)``: the largest
    ``TRUNCATION_ONE_SIDED * h`` over the rows whose ``scheme`` is
    ``'forward'`` or ``'backward'`` (first-order accurate), the largest
    ``TRUNCATION_CENTRAL * h**2`` over the central rows, each 0.0 without
    such a row, and the names of the one-sided rows; ``h`` is each row's
    :func:`relative_step` (rows without one, at ``value = 0``, are
    skipped). Needs no tolerances.
    """
    one_sided, central, names = 0.0, 0.0, []
    for name, row in flags.iterrows():
        h = relative_step(row['value'], row['log'], row['step'])
        first_order = row['scheme'] in ('forward', 'backward')
        if first_order:
            names.append(name)
        if not h > 0:                                   # NaN at value = 0
            continue
        if first_order:
            one_sided = max(one_sided, TRUNCATION_ONE_SIDED * h)
        else:
            central = max(central, TRUNCATION_CENTRAL * h ** 2)
    return one_sided, central, names


class ErrorLevel(NamedTuple):
    """Estimated error of finite differences relative to the derivatives
    (see :func:`error_level`).

    Attributes
    ----------
    noise : float
        Integration noise ``FD_NOISE * rtol / h_min`` (0.0 with unknown
        tolerances or no relative step).
    one_sided : float
        Truncation error of the reliable one-sided differences
        (``TRUNCATION_ONE_SIDED * h``, largest; 0.0 without one).
    central : float
        Truncation error of the reliable central differences
        (``TRUNCATION_CENTRAL * h**2``, largest; 0.0 without one).
    one_sided_names : list of str
        The reliable rows differenced on one side.
    h_min : float
        The smallest relative step of all rows (``NaN`` without one).
    """
    noise: float
    one_sided: float
    central: float
    one_sided_names: list
    h_min: float

    @property
    def total(self):
        """``noise + max(one_sided, central)``."""
        return self.noise + max(self.one_sided, self.central)


def error_level(flags, tolerances):
    """Estimated error of the differences in ``flags`` (an FD flags table,
    e.g. :attr:`FDJacobian.flags` or a selection of its rows) relative to
    the derivatives, given the integrator ``tolerances`` (``(atol, rtol)``
    or ``None`` if unknown).

    It is the level at which an exact non-identifiability shows up in the
    singular values of a Jacobian built from these differences, relative
    to the largest: integration noise ``FD_NOISE * rtol / h_min``, with
    ``h_min`` the smallest relative step (:func:`relative_step`) of all
    rows, plus the truncation error of the reliable rows
    (:func:`truncation_error`; an unreliable row's error is dominated by an
    event jump instead). The estimates are conservative on the sum toys
    measured (see :data:`FD_NOISE`, :data:`TRUNCATION_ONE_SIDED`); the
    total level was up to 24 x low on toys with large higher ``ln p``
    derivatives (see ``identifiability.FD_ERROR_MARGIN``), and far too low
    at depletion fronts or switches. The level counts ``rtol`` only and
    assumes an ``atol`` of at most :data:`ATOL_MAX` (``1e-8``; the default
    is 1e-12, scaled to the outputs); identifiability warns about a larger
    one.

    Returns
    -------
    ErrorLevel
    """
    steps = [relative_step(v, log, h) for v, log, h in
             zip(flags['value'], flags['log'], flags['step'])]
    steps = [h for h in steps if h > 0]          # drops NaN (value = 0)
    h_min = min(steps) if steps else math.nan
    noise = 0.0
    if tolerances is not None and h_min > 0:
        with np.errstate(over='ignore'):         # inf at a subnormal step
            noise = FD_NOISE * float(tolerances[1]) / h_min
    reliable = flags['reliable'].to_numpy(dtype=bool)
    one_sided, central, names = truncation_error(flags[reliable])
    return ErrorLevel(noise, one_sided, central, names, h_min)


def check_step(rel_step):
    """``rel_step`` as a float in ``(0, 1]``, or ``TypeError``/``ValueError``."""
    rel_step = _finite('rel_step', rel_step)
    if not 0 < rel_step <= 1:
        raise ValueError(f'rel_step must be in (0, 1]; got {rel_step!r}.')
    return rel_step


def check_abs_step(abs_step, names):
    """Validate ``abs_step`` (None, a positive float, or ``{name: float}``
    for names among ``names``); return it with floats."""
    if abs_step is None:
        return None
    if isinstance(abs_step, Mapping):
        unknown = [k for k in abs_step if k not in names]
        if unknown:
            raise ValueError(f'abs_step names {unknown} that are not among '
                             f'the parameters {list(names)}.')
        out = {k: _finite(f'abs_step[{k!r}]', v) for k, v in abs_step.items()}
        bad = [k for k, v in out.items() if v <= 0]
        if bad:
            raise ValueError(f'abs_step must be > 0; got {bad}.')
        return out
    out = _finite('abs_step', abs_step)
    if out <= 0:
        raise ValueError(f'abs_step must be > 0; got {abs_step!r}.')
    return out


def make_parameters(values, *, rel_step, abs_step=None):
    """:class:`FDParameter` objects for ``{id: nominal value}``.

    A positive value is differenced in log space with step ``rel_step``
    (``value * exp(+-rel_step)``). A value ``<= 0`` (log space impossible),
    or any id listed in a mapping ``abs_step``, is differenced in linear
    space with an absolute step: ``abs_step`` (its entry, or the scalar),
    else ``rel_step * |value|``, else -- for a value of exactly 0 --
    ``rel_step`` itself in the id's own units.

    Parameters
    ----------
    values : mapping of str to float
    rel_step : float
    abs_step : None, float or mapping, optional

    Returns
    -------
    list of FDParameter
    """
    out = []
    for name, value in values.items():
        listed = isinstance(abs_step, Mapping) and name in abs_step
        if value > 0 and not listed:
            out.append(FDParameter(name, name, value, True, rel_step))
            continue
        if listed:
            step = abs_step[name]
        elif abs_step is not None and not isinstance(abs_step, Mapping):
            step = abs_step
        else:
            step = rel_step * abs(value) if value != 0 else rel_step
        out.append(FDParameter(name, name, value, False, step))
    return out


# --- model set-up ---------------------------------------------------------------
def check_tolerances(tolerances):
    """``None`` or a validated ``(atol, rtol)`` pair of positive floats."""
    if tolerances is None:
        return None
    if (isinstance(tolerances, (str, bytes, Mapping))
            or not hasattr(tolerances, '__len__') or len(tolerances) != 2):
        raise TypeError(f'tolerances must be an (atol, rtol) pair; got '
                        f'{tolerances!r}.')
    atol = _finite('tolerances[0] (atol)', tolerances[0])
    rtol = _finite('tolerances[1] (rtol)', tolerances[1])
    if atol <= 0 or rtol <= 0:
        raise ValueError(f'tolerances must be positive; got {tolerances!r}.')
    return (atol, rtol)


def _has_tolerances(model):
    settings = model._te.getIntegrator().getSettings()
    return 'relative_tolerance' in settings and 'absolute_tolerance' in settings


def private_model(model, *, tolerances=None, rel_step=None, stacklevel=3):
    """A private copy of ``model`` for finite-difference runs.

    Parameters
    ----------
    model : KineticModel
        The caller's model; never modified.
    tolerances : (float, float), optional
        ``(atol, rtol)`` written to the copy's integrator. ``None``
        (default) tightens each to at most :data:`FD_TOLERANCES`
        (``(1e-12, 1e-10)``), keeping any tighter value of the model; an
        integrator without tolerances (fixed-step) is left as it is. On
        that default path a caller that knows its outputs should then
        tighten ``atol`` further with :func:`scale_atol_to_outputs` (by the
        smallest peak magnitude of the outputs that start at 0, since
        RoadRunner gives such a state a fixed absolute tolerance ``atol``
        for the whole run), as
        :func:`~nskinetics.analysis.local_sensitivity` does.
    rel_step : float, optional
        The FD step in log space; when given, a ``UserWarning`` is issued if
        the copy's relative tolerance exceeds ``rel_step / 100`` (integration
        error would then dominate the differences).
    stacklevel : int
        For the warning.

    Returns
    -------
    KineticModel
        ``model.copy()`` with the tolerances applied. When the copy's
        relative tolerance ends up below the model's own, its integrator's
        ``maximum_num_steps`` (if it has a positive one) is also raised by
        :data:`STEP_BUDGET_FACTOR` (10), since tighter tolerances need more
        internal steps; this changes no run that succeeds anyway. Its
        settings (global parameter values) as built are kept on it for
        :func:`fd_jacobian`, which writes them back before its runs, so a
        value an earlier run's custom ``reset_func`` wrote cannot leak into
        later runs. Its event counters are detected once and kept on it as
        ``_fd_counters`` (the default of :func:`evaluate` and
        :func:`fd_jacobian`), with the left-out candidates in
        ``_fd_unwatched`` (see :func:`detect_event_counters`; no warning is
        issued here). These records are private attributes, which
        ``copy()`` and a blueprint build do not carry: a model built that
        way (e.g. a process-pool worker of ``_map_with_model``) gets them
        at its first :func:`evaluate` or :func:`fd_jacobian` call, from its
        state at that moment.

    Raises
    ------
    ValueError
        If ``tolerances`` are given but the model's integrator has none
        (e.g. a fixed-step integrator); raised before copying.
    """
    tol = check_tolerances(tolerances)
    if tol is not None and not _has_tolerances(model):
        raise ValueError(
            f"tolerances= needs an integrator with tolerances; "
            f"{model._te.getIntegrator().getName()!r} has none.")
    private = model.copy()
    if _has_tolerances(private):
        integrator = private._te.getIntegrator()
        current = tolerances_of(private)
        if tol is None:
            tol = (min(current[0], FD_TOLERANCES[0]),
                   min(current[1], FD_TOLERANCES[1]))
        integrator.setValue('absolute_tolerance', tol[0])
        integrator.setValue('relative_tolerance', tol[1])
        if (tol[1] < current[1]
                and 'maximum_num_steps' in integrator.getSettings()):
            budget = int(integrator.getValue('maximum_num_steps'))
            if budget > 0:      # <= 0 has a special meaning in SUNDIALS
                integrator.setValue('maximum_num_steps',
                                    STEP_BUDGET_FACTOR * budget)
    # The settings as built (a private attribute: copy() does not carry it).
    private._fd_settings = private._user_parameter_values()
    # The event counters, detected once (private too, like _fd_settings):
    # evaluate() and fd_jacobian() default to them instead of re-reading the
    # SBML on every call. Left-out candidates are warned about by the caller.
    private._fd_counters, private._fd_unwatched = _scan_counters(private)
    if rel_step is not None:
        check_tolerance(private, rel_step, stacklevel=stacklevel + 1)
    return private


def tolerances_of(model):
    """``(atol, rtol)`` of ``model``'s integrator as floats, or ``None`` for
    an integrator without tolerances."""
    if not _has_tolerances(model):
        return None
    integrator = model._te.getIntegrator()
    return (float(integrator.getValue('absolute_tolerance')),
            float(integrator.getValue('relative_tolerance')))


def scale_atol_to_outputs(model, conditions, selections, *,
                          min_factor=ATOL_MIN_FACTOR):
    """Tighten the absolute tolerance of a private copy to the scale of the
    outputs that start at 0.

    RoadRunner multiplies a scalar ``absolute_tolerance`` by each state's
    start value, but a state that starts at 0 gets a fixed absolute
    tolerance ``atol`` (in its own units) for the whole run. An output that
    starts at 0 and stays small (a product at the micromolar scale, a trace
    by-product) is then integrated only to ``atol`` absolute, unless other
    states hold the step size down -- up to 10% error in its log
    sensitivities on a decoupled nanomolar species. This runs each
    condition once and, if a selected output starts at exactly 0 and peaks
    at ``m_k = max_t |y_k| < 1`` (model units), multiplies the scalar
    ``atol`` by ``m = min_k m_k`` (over outputs and conditions; at least
    ``min_factor``). Such outputs are then integrated to ``atol`` relative
    to their peak, as :func:`effective_atol` assumes.

    Parameters
    ----------
    model : KineticModel
        A private copy (e.g. from :func:`private_model`); its integrator's
        ``absolute_tolerance`` is changed in place. A scalar setting, so a
        blueprint build (a process worker) carries it.
    conditions : sequence of Condition
        The runs whose outputs are measured (e.g. several experiments).
    selections : sequence of str
        The output selections to consider (every recorded column, e.g. the
        outputs and the selections that metrics read).
    min_factor : float
        Lower clip of the factor ``m`` (default
        :data:`ATOL_MIN_FACTOR`, ``1e-15``).

    Returns
    -------
    (float, float) or None
        The ``(atol, rtol)`` in force afterwards; ``None`` for an integrator
        without tolerances, which is left as it is.

    Notes
    -----
    Costs one run per condition (none without selections or tolerances).
    If a run fails, the tolerance is left unchanged: the analysis's own
    nominal run then reports the failure. An output that starts nonzero is
    left to RoadRunner's start-value scaling; one that stays 0, or is not
    finite, is ignored.
    """
    tol = tolerances_of(model)
    selections = list(dict.fromkeys(selections))
    if tol is None or not selections:
        return tol
    factor = 1.0
    for cond in conditions:
        try:
            result = model.run(cond, extra_selections=tuple(selections))
        except KineticSimulationError:
            return tol
        y = np.asarray(result[selections], dtype=float).reshape(
            len(result), len(selections))
        if not len(y):
            continue
        with np.errstate(all='ignore'):
            peak = np.abs(y).max(axis=0)
            zero = (y[0] == 0) & np.isfinite(peak) & (peak > 0)
        if zero.any():
            factor = min(factor, float(peak[zero].min()))
    factor = max(factor, min_factor)
    if factor < 1:
        model._te.getIntegrator().setValue('absolute_tolerance',
                                           tol[0] * factor)
    return tolerances_of(model)


def failure_hint(model):
    """The remedies appended to a failed nominal or start-value run of an
    analysis: the integrator tolerances in use and how to change them (``''``
    for an integrator without tolerances)."""
    tol = tolerances_of(model)
    if tol is None:
        return ''
    return (f' The analysis runs used integrator tolerances (atol, rtol) = '
            f'{tol}; finite-difference analyses tighten them to at most '
            f'{FD_TOLERANCES} by default, which needs more integrator steps. '
            'If the model runs at its own tolerances, pass tolerances=(atol, '
            "rtol) to loosen them, or raise the integrator's "
            'maximum_num_steps (or use a denser output grid).')


def check_tolerance(model, rel_step, *, stacklevel=2):
    """Warn (``UserWarning``) if ``model``'s integrator relative tolerance
    exceeds ``rel_step / 100``; integrators without tolerances are skipped.
    Only gross mismatches are caught (see :data:`TOLERANCE_RATIO`)."""
    tol = tolerances_of(model)
    if tol is None:
        return
    rtol = tol[1]
    limit = rel_step / TOLERANCE_RATIO
    if rtol > limit:
        warnings.warn(
            f'The integrator relative tolerance ({rtol:g}) exceeds '
            f'rel_step/100 ({limit:g}), so integration error may dominate '
            f'the finite differences; pass tolerances=(atol, rtol) with '
            f'rtol <= {limit:g}, or a larger rel_step (a larger step '
            'raises the truncation error, which matters for '
            'identifiability).',
            UserWarning, stacklevel=stacklevel)


def noisy_parameters(flags, tolerances, *, checked=None):
    """Names in ``flags`` (an FD flags table) whose reported difference has
    a relative step below ``TOLERANCE_RATIO * rtol``: integration error,
    amplified by ``rtol / step``, may dominate them (e.g. after step
    refinement, or with a small ``abs_step``). ``[]`` without tolerances.

    ``checked`` (optional) is a relative step already checked against the
    tolerance (e.g. ``rel_step`` by :func:`private_model`): rows whose
    relative step is at least that are skipped, so a mismatch is reported
    once.
    """
    if tolerances is None:
        return []
    rtol = tolerances[1]
    out = []
    for name, row in flags.iterrows():
        h = relative_step(row['value'], row['log'], row['step'])
        if checked is not None and not h < checked * (1 - 1e-9):
            continue                            # checked already, or NaN
        if h < TOLERANCE_RATIO * rtol:          # False for NaN
            out.append(name)
    return out


def effective_atol(atol, y0):
    """The integrator's absolute tolerance in the units of each output.

    RoadRunner multiplies a scalar ``absolute_tolerance`` by the value of
    each state element at the start of a run (a species amount, or a
    rate-rule variable); a species starting at 0 gets ``atol`` times its
    compartment volume (``atol`` as a concentration), a non-species state
    starting at 0 gets ``atol`` itself. So the integration error of an
    output that starts at ``y0 != 0`` scales with ``|y0|``, not with
    ``|y|`` alone. For an output that starts at 0, ``atol`` is an upper
    bound that is only tight when nothing else holds the integrator's
    step size down; by default
    :func:`~nskinetics.analysis.local_sensitivity` scales ``atol``
    to the peak of such outputs (:func:`scale_atol_to_outputs`), so the
    bound is then relative to their own size.

    Parameters
    ----------
    atol : float
        The scalar absolute tolerance of the integrator.
    y0 : array_like
        Start values of the outputs (model units), e.g. the first row of
        the nominal run.

    Returns
    -------
    numpy.ndarray
        ``atol * |y0|`` where ``y0`` is finite and nonzero, else ``atol``,
        elementwise (the shape of ``y0``), in the outputs' units.
    """
    y0 = np.abs(np.asarray(y0, dtype=float))
    with np.errstate(all='ignore'):
        use = np.isfinite(y0) & (y0 > 0)
        return np.where(use, atol * np.where(use, y0, 1.0), atol)


def _is_increment(ast, var):
    """Whether the AST is exactly ``var + 1`` (or ``1 + var``).

    Uses the predicates (``getCharacter``, ``isName``, ``isNumber``) rather
    than ``getType()``: once RoadRunner's SWIG module is loaded,
    ``ASTNode.getType()`` returns an opaque SWIG pointer, not the enum int,
    so comparing it with ``libsbml.AST_PLUS`` is always false.
    """
    if (ast is None or not ast.isOperator() or ast.getCharacter() != '+'
            or ast.getNumChildren() != 2):
        return False
    a, b = ast.getChild(0), ast.getChild(1)
    for name, number in ((a, b), (b, a)):
        if (name.isName() and name.getName() == var
                and number.isNumber() and number.getValue() == 1):
            return True
    return False


def _changed_outside_events(sbml_model):
    """Ids that change other than by events: non-boundary species that are
    a reactant or product of a reaction, and rate- or assignment-rule
    variables."""
    ids = set()
    for reaction in sbml_model.getListOfReactions():
        for ref in (*reaction.getListOfReactants(),
                    *reaction.getListOfProducts()):
            species = sbml_model.getSpecies(ref.getSpecies())
            if species is not None and not species.getBoundaryCondition():
                ids.add(ref.getSpecies())
    for rule in sbml_model.getListOfRules():
        if rule.isRate() or rule.isAssignment():
            ids.add(rule.getVariable())
    return ids


def _variable_compartments(sbml_model):
    """Ids of compartments whose size can change: declared non-constant, or
    the variable of a rule or an event assignment."""
    ids = {c.getId() for c in sbml_model.getListOfCompartments()}
    out = {c.getId() for c in sbml_model.getListOfCompartments()
           if not c.getConstant()}
    out |= {r.getVariable() for r in sbml_model.getListOfRules()} & ids
    for event in sbml_model.getListOfEvents():
        out |= {ea.getVariable()
                for ea in event.getListOfEventAssignments()} & ids
    return out


def _scan_counters(model):
    """``(counters, unwatched)``: the event counters of ``model`` (see
    :func:`detect_event_counters`), and the candidates left out because
    they are concentration species in a compartment whose size changes."""
    if not model._sbml_event_targets():
        return (), ()
    # Keep the document referenced while its model is used: libsbml frees a
    # temporary document's children with it (dangling ASTs otherwise).
    document = libsbml.readSBMLFromString(model._get_origin_sbml())
    sbml_model = document.getModel()
    if sbml_model is None:
        return (), ()
    continuous = _changed_outside_events(sbml_model)
    variable = _variable_compartments(sbml_model)
    counters, unwatched, seen = [], [], set()
    for event in sbml_model.getListOfEvents():
        for ea in event.getListOfEventAssignments():
            var = ea.getVariable()
            if (not var or var in seen or var in continuous
                    or not _is_increment(ea.getMath(), var)):
                continue
            seen.add(var)
            species = sbml_model.getSpecies(var)
            if species is None or species.getHasOnlySubstanceUnits():
                counters.append(var)        # a parameter, or an amount
            elif species.getCompartment() in variable:
                # A concentration species in a growing compartment: v = v + 1
                # adds 1 to the concentration, so neither its amount (count x
                # the volume at each firing) nor its concentration (diluted
                # as the volume changes) is a count.
                unwatched.append(var)
            else:
                # In a constant compartment its concentration rises by
                # exactly 1 per firing, whatever the (possibly perturbed)
                # size; its amount, count x size, would not.
                counters.append(f'[{var}]')
    return tuple(counters), tuple(unwatched)


def warn_unwatched_counters(ids, *, stacklevel=3):
    """Warn (``UserWarning``) that the counter candidates ``ids`` (see
    :func:`detect_event_counters`) are not watched; nothing for none."""
    if not ids:
        return
    warnings.warn(
        f'{list(ids)} are incremented by events (v = v + 1) but are '
        'concentration species in a compartment whose size changes, so their '
        'recorded amount is the count times the volume at each firing, not a '
        'count; they are not watched as event counters. Count with a '
        "parameter (as a FeedSpike's count_var does), or pass "
        'event_counters= explicitly.', UserWarning, stacklevel=stacklevel)


def detect_event_counters(model, *, warn=True, stacklevel=2):
    """Event counters of ``model``: variables that only events change, and
    that an event sets to exactly ``v + 1``.

    Read from the model's origin SBML, so it covers compiled
    :class:`~nskinetics.Event` / :class:`~nskinetics.FeedSpike` objects
    (e.g. a FeedSpike's ``count_var``) and events written in the
    Antimony/SBML source alike; pending (not yet compiled) events never fire
    and are ignored. Touches no model state.

    Parameters
    ----------
    model : KineticModel
    warn : bool
        Warn about candidates that are left out (see Notes).
    stacklevel : int
        For the warning.

    Returns
    -------
    tuple of str
        Counter selections in document order: plain ids, except a
        concentration species counter, which is read as its concentration
        ``'[n]'`` (see Notes).

    Warns
    -----
    UserWarning
        If ``warn`` and a candidate is left out because it is a species in
        a compartment whose size changes (see Notes).

    Notes
    -----
    A variable that also changes outside events -- a (non-boundary) species
    that takes part in a reaction, or a rate- or assignment-rule variable --
    is not a counter even if an event adds 1 to it (e.g. a pulse feed
    ``S = S + 1`` of a substrate): its final value differs between runs
    whether or not an event count changed. Events that change state
    without counting (e.g. a stage switch ``is_aerobic = 0``) are not
    detected either; pass ``event_counters=`` explicitly to watch such
    variables.

    A species with concentration semantics (``hasOnlySubstanceUnits``
    false, the Antimony default) in a compartment whose size changes (a
    non-constant compartment, or one set by a rule or an event) is left out
    with a warning: ``n = n + 1`` adds 1 to its concentration, so its
    amount is the count times the volume at each firing time, which moves
    with any parameter that moves the volume or the firing time. Count
    with a parameter instead, or list it in ``event_counters=``. In a
    constant compartment the same species is watched through its
    concentration ``'[n]'``, which rises by exactly 1 per firing whatever
    the compartment's size (its amount would change with a perturbed
    size); a species with ``hasOnlySubstanceUnits`` true is watched
    through its plain id, its amount.
    """
    counters, unwatched = _scan_counters(model)
    if warn:
        warn_unwatched_counters(unwatched, stacklevel=stacklevel + 1)
    return counters


def _names_in(ast, out):
    """Collect the names (``isName()`` nodes) of an AST into ``out``."""
    if ast is None:
        return
    if ast.isName():
        out.add(ast.getName())
    for i in range(ast.getNumChildren()):
        _names_in(ast.getChild(i), out)


def initial_assignment_inputs(model):
    """Which ids the model's SBML initial assignments read.

    Parameters
    ----------
    model : KineticModel

    Returns
    -------
    dict
        ``{id: (symbol, ...)}``: each id read by the math of an initial
        assignment, with the symbols (species, compartments, parameters)
        whose initial values those assignments compute, in document order.

    Notes
    -----
    Reads the origin SBML (see :func:`detect_event_counters`); touches no
    model state. An analysis that sets such an id after the reset (through
    :attr:`Condition.parameters` or :attr:`Condition.initial`) does not
    re-evaluate the assignments; setting ``'init(<id>)'`` does.
    """
    # Keep the document referenced while its model is used (see above).
    document = libsbml.readSBMLFromString(model._get_origin_sbml())
    sbml_model = document.getModel()
    if sbml_model is None:
        return {}
    readers = {}
    for ia in sbml_model.getListOfInitialAssignments():
        names = set()
        _names_in(ia.getMath(), names)
        symbol = ia.getSymbol()
        for name in sorted(names):
            mine = readers.setdefault(name, [])
            if symbol not in mine:
                mine.append(symbol)
    return {name: tuple(symbols) for name, symbols in readers.items()}


def warn_initial_assignment_inputs(model, targets, *, stacklevel=3):
    """Warn (``UserWarning``) about ``targets`` that SBML initial
    assignments read.

    Such a target is set after the reset has evaluated the assignments, so
    an analysis that varies it holds the assigned initial values fixed: its
    effect through them is missing (a sensitivity of 0, an inert parameter
    in screening). ``'init(<id>)'`` targets are re-evaluated and skipped.

    Parameters
    ----------
    model : KineticModel
    targets : iterable of str
        Ids the analysis varies (``'[S]'`` is matched as ``S``).
    stacklevel : int
        For the warning.

    Returns
    -------
    list of str
        The targets warned about.
    """
    plain = {t: t[1:-1] if t.startswith('[') and t.endswith(']') else t
             for t in targets if not t.startswith('init(')}
    if not plain:
        return []
    readers = initial_assignment_inputs(model)
    hit = {t: readers[i] for t, i in plain.items() if i in readers}
    if hit:
        detail = '; '.join(f'{t!r} is read by the initial assignment(s) of '
                           f'{list(s)}' for t, s in hit.items())
        suggest = ', '.join(repr(f'init({plain[t]})') for t in hit)
        warnings.warn(
            f'{detail}. It is set after the reset has computed those initial '
            f'values, so the result holds them fixed and misses its effect '
            f'through them; vary {suggest} instead to include it.',
            UserWarning, stacklevel=stacklevel)
    return list(hit)


def check_counters(event_counters):
    """``event_counters`` as a tuple of unique non-empty strings."""
    if isinstance(event_counters, str):
        event_counters = (event_counters,)
    if isinstance(event_counters, (bytes, Mapping)) or not hasattr(
            event_counters, '__iter__'):
        raise TypeError('event_counters must be a sequence of ids; got '
                        f'{event_counters!r}.')
    out = []
    for c in event_counters:
        if not isinstance(c, str) or not c:
            raise TypeError(f'event_counters must hold non-empty strings; got '
                            f'{c!r}.')
        if c not in out:
            out.append(c)
    return tuple(out)


def nominal_values(model, condition, targets):
    """The value of each target in a run of ``condition``.

    A target the condition sets (in ``initial`` or ``parameters``) has that
    value. Any other is read from the first row (``t0``, after the reset and
    the condition's values are applied) of one run of ``condition`` on
    ``model`` that records it, so it is exactly what that run starts from --
    whatever the reset, a custom ``reset_func`` or the condition's own
    entries (a volume, an ``init(x)``) made of it.

    Parameters
    ----------
    model : KineticModel
        A private model (it is run).
    condition : Condition
    targets : sequence of str

    Returns
    -------
    dict
        ``{target: float}`` in the order of ``targets``.

    Raises
    ------
    KineticSimulationError
        If that run fails; the message names the integrator tolerances in
        use and the remedies (see :func:`failure_hint`).
    """
    values, missing = {}, []
    for t in targets:
        if t in condition.initial:
            values[t] = float(condition.initial[t])
        elif t in condition.parameters:
            values[t] = float(condition.parameters[t])
        else:
            missing.append(t)
    if missing:
        try:
            first = model.run(condition, extra_selections=missing).df.iloc[0]
        except KineticSimulationError as e:
            label = condition.label
            raise KineticSimulationError(
                f'The start-value run of the condition'
                f'{f" ({label!r})" if label else ""}, which reads the '
                f'nominal values of {missing}, failed: '
                f'{str(e).rstrip(". ")}.' + failure_hint(model)) from e
        for t in missing:
            values[t] = float(first[t])
    return {t: values[t] for t in targets}


def pin(model, condition, values, placement=None):
    """``condition`` with ``values`` merged in (``with_parameters`` /
    ``with_initial``), each id placed like a :func:`nskinetics.scan` key
    (``placement``: precomputed ``{id: 'initial'|'parameters'}``)."""
    initial, parameters = {}, {}
    for target, value in values.items():
        where = (placement[target] if placement is not None
                 else _ex._placement_for(model, condition, target))
        (initial if where == 'initial' else parameters)[target] = value
    if parameters:
        condition = condition.with_parameters(parameters)
    if initial:
        condition = condition.with_initial(initial)
    return condition


def flags_frame(rows, names):
    """The per-parameter flags table: one row per name, the columns and
    dtypes of :data:`FLAG_COLUMNS`, index named ``'parameter'``."""
    data = {c: pd.Series([row[c] for row in rows], dtype=dtype)
            for c, dtype in FLAG_COLUMNS.items()}
    frame = pd.DataFrame(data)
    frame.index = pd.Index(list(names), name='parameter', dtype=object)
    return frame


# --- the engine -----------------------------------------------------------------
@dataclass(eq=False)
class FDJacobian:
    """Outcome of :func:`fd_jacobian`.

    Attributes
    ----------
    jacobian : numpy.ndarray
        ``(n_out, n_params)``: ``df / dx_j`` with ``x_j = ln p_j`` for a
        log-space parameter, ``p_j`` otherwise (the calibration Jacobian in
        transformed space).
    f0 : numpy.ndarray
        ``(n_out,)``: ``fn`` at the nominal point.
    parameters : tuple of FDParameter
        As given (``step`` is the initial step; the one used is in
        ``flags``).
    flags : pandas.DataFrame
        One row per parameter (index ``parameter``): ``value``, ``log``,
        ``step`` (the step of the reported difference, in the differencing
        space), ``n_refinements``, ``event_count_changed``, ``reliable``,
        ``scheme`` (``'central'``, ``'forward'`` or ``'backward'``: the
        difference actually used; one-sided for a central difference after
        a count change or a failed run on the other side) and ``error``
        (``None``, or why the column is ``NaN``: every run it needed
        failed, or ``fn`` failed on them).
    counters : tuple of str
        Event counters watched.
    conditions : list of Condition
        The nominal conditions (every applicable target pinned).
    results : list of SimulationResult
        The nominal runs.
    tolerances : (float, float) or None
        ``(atol, rtol)`` of the model's integrator during the runs
        (``None`` for an integrator without tolerances).
    """
    jacobian: np.ndarray
    f0: np.ndarray
    parameters: tuple
    flags: pd.DataFrame
    counters: tuple
    conditions: list
    results: list
    tolerances: tuple | None = None

    @property
    def natural(self):
        """``df / dp`` (``numpy.ndarray``, ``(n_out, n_params)``): the
        Jacobian with log-space columns divided by the parameter value."""
        scale = np.array([p.value if p.log else 1.0 for p in self.parameters])
        with np.errstate(all='ignore'):
            return self.jacobian / scale


@dataclass(eq=False)
class FDPoint:
    """Outcome of :func:`evaluate`: the runs and ``fn`` at the nominal point.

    Pass it to :func:`fd_jacobian` as ``nominal=`` (same model, conditions,
    parameters, ``fn`` and options) to skip the nominal runs there.

    Attributes
    ----------
    f : numpy.ndarray or None
        ``fn(results, values)`` (1-D), or ``None`` if a run failed.
    values : dict
        ``{name: natural value}`` of every parameter, as ``fn`` got it.
    parameters : tuple of FDParameter
    conditions : list of Condition
        The conditions run (every applicable target pinned).
    results : list
        One :class:`~nskinetics.SimulationResult` per condition (``None``
        where the run failed).
    status : list of str
        Per condition, ``'ok'`` or ``'failed'``.
    errors : list
        Per condition, ``None`` or the error message.
    counters : tuple of str
        Event counters recorded.
    counts : list
        Per condition, the recorded counter columns on the output grid
        (``(n_t, n_counters)`` array; ``None`` for a failed run or without
        counters).
    extra_selections : tuple of str
        Columns recorded beyond the conditions' selections (counters
        included).
    settings : dict
        Global parameter values every run started from.
    tolerances : (float, float) or None
        ``(atol, rtol)`` of the integrator during the runs.
    """
    f: np.ndarray | None
    values: dict
    parameters: tuple
    conditions: list
    results: list
    status: list
    errors: list
    counters: tuple
    counts: list
    extra_selections: tuple
    settings: dict
    tolerances: tuple | None = None

    @property
    def ok(self):
        """Whether every run succeeded (``f`` is then set)."""
        return self.f is not None

    @property
    def error(self):
        """``None``, or the errors of the failed runs in one message."""
        failed = [f'condition {i}: {err}' for i, (status, err)
                  in enumerate(zip(self.status, self.errors))
                  if status != 'ok']
        return '; '.join(failed) or None


@dataclass
class _Side:
    """One perturbed evaluation: ``f``, or the error, and whether every
    watched counter counted the same firings as the nominal run by every
    output time (see :func:`_same_trace`)."""
    f: np.ndarray | None
    error: str | None
    same: bool


@dataclass(eq=False)
class _Setup:
    """The validated arguments shared by :func:`evaluate` and
    :func:`fd_jacobian`, with the nominal assignment and the placement of
    every target per condition (computed once)."""
    conditions: list
    parameters: tuple
    counters: tuple
    extra: tuple
    settings: dict
    workers: int
    assign: list
    place: list
    nominal: list

    def runner(self, model):
        """``run(conditions) -> [(status, error, result)]`` on ``model``."""
        return functools.partial(_run_batch, model, extra=self.extra,
                                 workers=self.workers, settings=self.settings)

    def perturbed(self, model, j, value):
        """``[(i, condition)]``: the conditions parameter ``j`` applies to,
        with it at ``value`` (the others pinned at their nominal values)."""
        p = self.parameters[j]
        return [(i, pin(model, c, {**self.assign[i], p.target: value},
                        self.place[i]))
                for i, c in enumerate(self.conditions) if p.applies_to(i)]


def _setup(model, conditions, parameters, *, event_counters, workers,
           extra_selections, settings):
    """Validate the shared arguments before any run (``TypeError`` /
    ``ValueError``), pin the nominal conditions and write the settings back
    on ``model``."""
    if isinstance(conditions, Condition):
        conditions = [conditions]
    conditions = list(conditions)
    if not conditions:
        raise ValueError('At least one condition is needed.')
    parameters = tuple(parameters)
    for p in parameters:
        if not isinstance(p, FDParameter):
            raise TypeError(f'parameters must be FDParameter objects; got '
                            f'{type(p).__name__}.')
    names = [p.name for p in parameters]
    if len(set(names)) != len(names):
        raise ValueError(f'Duplicate parameter names in {names}.')
    n_c = len(conditions)
    for p in parameters:
        if p.conditions is not None and p.conditions[-1] >= n_c:
            raise ValueError(f'{p.name!r} applies to condition '
                             f'{p.conditions[-1]}, but there are {n_c}.')
    workers = _ex._check_workers(workers)
    # A model that did not come from private_model (a copy() or a blueprint
    # build, e.g. a process-pool worker) has no record yet: the first call's
    # snapshot and detection are kept on it below, once the arguments are
    # valid, so later calls start from the same settings (a value an earlier
    # run's reset_func wrote cannot leak) and do not re-read the SBML.
    cache_settings = cache_counters = False
    if settings is None:
        settings = getattr(model, '_fd_settings', None)
        if settings is None:
            settings = model._user_parameter_values()
            cache_settings = True
    elif not isinstance(settings, Mapping):
        raise TypeError(f'settings must be a mapping; got '
                        f'{type(settings).__name__}.')
    settings = dict(settings)
    if event_counters is not None:
        counters = check_counters(event_counters)
    else:
        # Detected once by private_model (or a first call); else now.
        counters = getattr(model, '_fd_counters', None)
        if counters is None:
            counters, unwatched = _scan_counters(model)
            warn_unwatched_counters(unwatched, stacklevel=4)
            cache_counters = True
    extra = tuple(dict.fromkeys([*extra_selections, *counters]))

    # Nominal assignments and placements per condition (computed once).
    assign, place = [], []
    for i, c in enumerate(conditions):
        mine = [p for p in parameters if p.applies_to(i)]
        targets = [p.target for p in mine]
        dup = sorted({t for t in targets if targets.count(t) > 1})
        if dup:
            raise ValueError(f'Condition {i} has several parameters setting '
                             f'{dup}.')
        assign.append({p.target: p.value for p in mine})
        place.append({t: _ex._placement_for(model, c, t) for t in targets})
    nominal = [pin(model, c, assign[i], place[i])
               for i, c in enumerate(conditions)]
    # Order independence and ids, as for a batch; then the extra columns.
    _ex._validate_batch(model, nominal, None, workers, True, 'record')
    if extra:
        model._validate_ids(extra, settable=False)
    # Every argument is valid: record what was derived here (never an
    # explicit settings= or event_counters=; see above).
    if cache_settings:
        model._fd_settings = dict(settings)
    if cache_counters:
        model._fd_counters, model._fd_unwatched = counters, unwatched
    # Start from the settings as built: serial runs restore them before each
    # run, and a process pool builds its models from this state.
    _ex._restore_settings(model, settings)
    return _Setup(conditions=conditions, parameters=parameters,
                  counters=counters, extra=extra, settings=settings,
                  workers=workers, assign=assign, place=place,
                  nominal=nominal)


def _run_batch(model, conditions, *, extra, workers, settings):
    """``[(status, error, result)]`` for ``conditions``: serially on
    ``model`` (``workers=1``), else on F1's process pool."""
    fn = functools.partial(_ex._execute, metrics={}, extra_selections=extra,
                           keep_results=True, on_error='record')
    if workers == 1:
        out = []
        for condition in conditions:
            _ex._restore_settings(model, settings)
            out.append(fn(model, condition))
    else:
        out = _ex._map_with_model(model, fn, conditions, workers=workers)
    return [(status, error, result) for status, error, _, result in out]


def _vector(fn, results, values, size=None):
    """``fn(results, values)`` as a 1-D float array (of length ``size``)."""
    with np.errstate(all='ignore'):
        f = np.asarray(fn(results, values), dtype=float)
    if f.ndim != 1:
        raise ValueError(f'fn must return a 1-D array; got shape {f.shape}.')
    if size is not None and f.size != size:
        raise ValueError(f'fn returned {f.size} values; the nominal run gave '
                         f'{size}.')
    return f


def _counter_trace(result, counters):
    """The watched counters' recorded columns on the output grid
    (``(n_t, n_counters)`` float array), or ``None`` without counters."""
    if not counters:
        return None
    return result.df[list(counters)].to_numpy(dtype=float)


def _same_trace(trace, nominal):
    """Whether a perturbed run counted the same firings as the nominal run
    by every output time (``True`` without counters).

    Each run's counters are taken relative to its own first row (counter
    minus its start value), so perturbing a counter's start value is no
    change, while a firing added, removed or moved across an output time
    is. They are compared up to round-off (``1e-9`` of their magnitude).
    """
    if nominal is None:
        return True
    if trace.shape != nominal.shape:
        return False
    with np.errstate(all='ignore'):
        a, b = trace - trace[:1], nominal - nominal[:1]
        scale = 1.0 + max(np.abs(trace).max(initial=0.0),
                          np.abs(nominal).max(initial=0.0))
        return bool(np.allclose(a, b, rtol=0.0, atol=1e-9 * scale))


def _point(model, setup, fn, outputs):
    """The :class:`FDPoint` of the nominal ``outputs`` of ``setup``."""
    parameters = setup.parameters
    values = {p.name: p.value for p in parameters}
    status = [status for status, _, _ in outputs]
    errors = [err for _, err, _ in outputs]
    results = [result for _, _, result in outputs]
    ok = all(s == 'ok' for s in status)
    f = _vector(fn, results, values) if ok else None
    counts = [None if r is None else _counter_trace(r, setup.counters)
              for r in results]
    return FDPoint(f=f, values=values, parameters=parameters,
                   conditions=list(setup.nominal), results=results,
                   status=status, errors=errors, counters=setup.counters,
                   counts=counts, extra_selections=setup.extra,
                   settings=dict(setup.settings),
                   tolerances=tolerances_of(model))


def evaluate(model, conditions, parameters, fn, *, settings=None,
             event_counters=None, extra_selections=(), workers=1):
    """Run the nominal point of a finite-difference problem and evaluate
    ``fn`` there, without differencing.

    The runs are exactly the nominal runs of :func:`fd_jacobian` (every
    applicable target pinned at its nominal value, settings written back
    before each serial run, counters and extra columns recorded), so a
    function value and a Jacobian computed from the same arguments agree --
    calibration's residual function and its Jacobian.

    Parameters
    ----------
    model : KineticModel
        A *private* model (see :func:`fd_jacobian`).
    conditions : sequence of Condition
        The base conditions.
    parameters : sequence of FDParameter
        The point: each parameter's ``value`` (``log`` and ``step`` are not
        used here).
    fn : callable
        ``fn(results, values) -> 1-D array``, as for :func:`fd_jacobian`;
        evaluated under ``numpy.errstate(all='ignore')`` only if every run
        succeeded.
    settings, event_counters, extra_selections, workers
        As for :func:`fd_jacobian`.

    Returns
    -------
    FDPoint
        A failed run is recorded (``status``, ``errors``; ``f`` is
        ``None``), never raised, so a caller can apply a penalty.

    Raises
    ------
    TypeError, ValueError
        For invalid arguments, before any run.
    Exception
        Whatever ``fn`` raises on the results.
    """
    if not callable(fn):
        raise TypeError('fn must be callable.')
    setup = _setup(model, conditions, parameters,
                   event_counters=event_counters, workers=workers,
                   extra_selections=extra_selections, settings=settings)
    return _point(model, setup, fn, setup.runner(model)(setup.nominal))


def _check_nominal(nominal, setup, model):
    """``ValueError`` unless ``nominal`` is a successful :class:`FDPoint` of
    the same problem as ``setup``."""
    if not isinstance(nominal, FDPoint):
        raise TypeError(f'nominal must be an FDPoint (from evaluate); got '
                        f'{type(nominal).__name__}.')
    if not nominal.ok:
        raise ValueError(f'nominal is a failed evaluation ({nominal.error}); '
                         'a Jacobian needs successful nominal runs.')
    key = lambda ps: [(p.name, p.target, p.value, p.conditions)  # noqa: E731
                      for p in ps]
    mismatch = [what for what, same in (
        ('parameters', key(nominal.parameters) == key(setup.parameters)),
        ('conditions', nominal.conditions == setup.nominal),
        ('event counters', nominal.counters == setup.counters),
        ('extra_selections', nominal.extra_selections == setup.extra),
        ('settings', nominal.settings == setup.settings),
        ('tolerances', nominal.tolerances == tolerances_of(model)))
        if not same]
    if mismatch:
        raise ValueError(f'nominal was evaluated for different {mismatch}; '
                         'pass the same arguments to evaluate and '
                         'fd_jacobian.')


def fd_jacobian(model, conditions, parameters, fn, *, method='central',
                event_counters=None, max_refinements=2, workers=1,
                extra_selections=(), settings=None, nominal=None):
    """Event-aware finite-difference Jacobian of ``fn`` over conditions.

    Parameters
    ----------
    model : KineticModel
        A *private* model (e.g. from :func:`private_model`): runs use it
        directly with ``workers=1`` (each from its baseline), or its
        blueprint for ``workers`` processes. Its settings are first written
        back to ``settings`` (and again before every serial run). Its
        definition is not changed.
    conditions : sequence of Condition
        The base conditions (``reset=True``; not ``mode='roadrunner'``).
    parameters : sequence of FDParameter
        What to perturb. Every applicable target is pinned at its nominal
        value in every run; parameter ``j`` is then moved alone.
    fn : callable
        ``fn(results, values) -> 1-D array``: ``results`` holds one
        :class:`~nskinetics.SimulationResult` per condition (conditions a
        perturbation does not apply to get their nominal result) and
        ``values`` is ``{name: natural value}`` of every parameter at that
        evaluation point (what ``fn`` needs for ``target=None``
        parameters). Evaluated in this process under
        ``numpy.errstate(all='ignore')``; it should depend on the results
        only through recorded columns.
    method : {'central', 'forward'}
    event_counters : sequence of str, optional
        Counters whose firings must agree with the nominal run at every
        output time; ``()`` disables the check. ``None`` (default): the
        counters detected when :func:`private_model` built ``model``
        (``model._fd_counters``); on a model that did not come from
        :func:`private_model` (a ``copy()`` or a blueprint build, e.g. a
        process-pool worker), :func:`detect_event_counters` at the first
        call of this function or :func:`evaluate`, recorded on the model
        and reused by later calls.
    max_refinements : int
        How many times the step may be divided by 10 when counters differ.
    workers : int
        Processes (``1``: serial, no pool).
    extra_selections : sequence of str
        Columns ``fn`` needs beyond the conditions' selections.
    settings : dict, optional
        ``{global parameter: value}`` every run starts from. Defaults to the
        settings of the private model as built by :func:`private_model`,
        so a value written by an earlier run's custom ``reset_func`` (e.g.
        for other ``reset_kwargs``) does not leak into these runs. On a
        model that did not come from :func:`private_model` (a ``copy()`` or
        a blueprint build, e.g. a process-pool worker), the settings are
        recorded on the model at the first call of this function or
        :func:`evaluate`, and later calls start from that record; take that
        first call before any other run on the model. A caller that
        changes global settings deliberately between calls must pass
        ``settings=``.
    nominal : FDPoint, optional
        The result of :func:`evaluate` with the same model, conditions,
        parameters, ``fn`` and options (e.g. the residual evaluation of an
        optimizer at this point): its runs and ``f`` are used and the
        nominal runs are not repeated. ``None`` (default) runs them, in one
        batch with the first perturbations.

    Returns
    -------
    FDJacobian

    Raises
    ------
    TypeError, ValueError
        For invalid arguments (a failed or mismatched ``nominal``
        included), before any run.
    KineticSimulationError
        If a nominal run fails; the message names the integrator tolerances
        in use and the remedies (see :func:`failure_hint`).
    Exception
        Whatever ``fn`` raises on the nominal results. A failed perturbed
        run (or ``fn`` failing on it) does not raise: with
        ``method='central'`` the column is the one-sided difference on the
        other side if that side ran and its firings agree with the nominal
        run (``scheme`` says which side, the column stays reliable); if
        they do not agree, the step is refined as below (both sides re-run)
        and, should the running side still disagree, the column is flagged
        ``event_count_changed``. Otherwise it is ``NaN`` with
        ``reliable=False`` and the error in ``flags``; with
        ``method='forward'`` it is ``NaN``.

    Notes
    -----
    The watched counters are compared as the firings they count by each
    output time (each run's counter minus its own start value, up to
    round-off), so perturbing a counter's own start value is no change. A
    perturbed run that counted different firings by any output time -- a
    firing added or removed, or moved across an output time, where the
    difference on that row measures the jump -- has the step of that
    parameter divided by 10 (both sides re-run) up to ``max_refinements``
    times. If they still differ and exactly one side agrees with the
    nominal run, the column is the one-sided difference on that side at
    the *largest* step at which it agreed (``flags['step']``; with
    ``method='forward'`` the backward runs are made for the whole ladder of
    steps in one batch), since integration error grows as the step shrinks.
    The engine never warns (it runs inside optimizers); callers read
    ``flags``.

    Steps are not limited to a parameter's search bounds: a run may be made
    up to ``step`` outside them (``value * exp(+-step)`` in log space,
    ``value +- step`` in linear space), so choose ``FDParameter.step``
    accordingly near a bound.
    """
    if method not in METHODS:
        raise ValueError(f'method must be one of {METHODS}; got {method!r}.')
    if not callable(fn):
        raise TypeError('fn must be callable.')
    parameters = tuple(parameters)
    if not parameters:
        raise ValueError('fd_jacobian needs at least one parameter.')
    max_refinements = _count('max_refinements', max_refinements, 0)
    setup = _setup(model, conditions, parameters,
                   event_counters=event_counters, workers=workers,
                   extra_selections=extra_selections, settings=settings)
    if nominal is not None:
        _check_nominal(nominal, setup, model)
    parameters, counters = setup.parameters, setup.counters
    names = [p.name for p in parameters]
    n_c = len(setup.conditions)
    run = setup.runner(model)
    signs = (1, -1) if method == 'central' else (1,)
    point = {p.name: p.value for p in parameters}

    def jobs(key):
        """``[(key, i, condition)]``: the runs of ``key = (j, sign, step)``
        (parameter ``j`` moved by ``sign * step``)."""
        j, sign, step = key
        return [(key, i, c) for i, c in setup.perturbed(
            model, j, parameters[j].at(sign * step))]

    def perturb(keys):
        """Run the perturbations ``keys``; ``{key: _Side}``."""
        batch = [job for key in keys for job in jobs(key)]
        return evaluate_sides(keys, batch,
                              run([job[2] for job in batch]) if batch else [])

    def evaluate_sides(keys, batch, outputs):
        """``{key: _Side}`` from the ``outputs`` of ``batch``. A key without
        runs (a ``target=None`` parameter) re-evaluates ``fn`` on the
        nominal results."""
        runs = {key: {} for key in keys}
        for (key, i, _), out in zip(batch, outputs):
            runs[key][i] = out
        new = {}
        for key, by_condition in runs.items():
            j, sign, step = key
            errors = [f'condition {i}: {err}'
                      for i, (status, err, _) in by_condition.items()
                      if status != 'ok']
            if errors:
                new[key] = _Side(None, '; '.join(errors), False)
                continue
            results = [by_condition[i][2] if i in by_condition
                       else nominal_r[i] for i in range(n_c)]
            values = {**point, parameters[j].name: parameters[j].at(
                sign * step)}
            try:
                f = _vector(fn, results, values, f0.size)
            except Exception as e:  # recorded: this column becomes NaN
                new[key] = _Side(None, f'fn: {_ex._describe(e)}', False)
                continue
            # The firings counted must match the nominal run on the whole
            # output grid: a firing moved across an output time is a jump
            # there.
            same = all(_same_trace(_counter_trace(by_condition[i][2],
                                                  counters), counts[i])
                       for i in by_condition)
            new[key] = _Side(f, None, same)
        return new

    # sides[j, sign]: the evaluation at the current step h[j];
    # best[j, sign]: (side, step) at the largest step whose counters agreed
    # with the nominal run (steps only shrink, so the first one recorded).
    h = [p.step for p in parameters]
    sides, best = {}, {}

    def record(new):
        for (j, sign, step), side in new.items():
            sides[j, sign] = side
            if side.error is None and side.same:
                best.setdefault((j, sign), (side, step))

    # The nominal runs (unless given) and the first perturbations, in one
    # batch (one process pool with workers > 1).
    keys = [(j, s, h[j]) for j in range(len(parameters)) for s in signs]
    batch = [job for key in keys for job in jobs(key)]
    if nominal is None:
        outputs = run([*setup.nominal, *(job[2] for job in batch)])
        for i, (status, err, _) in enumerate(outputs[:n_c]):
            if status != 'ok':
                label = setup.conditions[i].label
                raise KineticSimulationError(
                    f'The nominal run of condition {i}'
                    f'{f" ({label!r})" if label else ""} failed: '
                    f'{str(err).rstrip(". ")}.' + failure_hint(model))
        nominal = _point(model, setup, fn, outputs[:n_c])
        outputs = outputs[n_c:]
    else:
        outputs = run([job[2] for job in batch]) if batch else []
    nominal_r, f0, counts = nominal.results, nominal.f, nominal.counts
    record(evaluate_sides(keys, batch, outputs))
    n_ref = [0] * len(parameters)

    def changed(j):
        """Whether a side that ran counted other firings than the nominal
        run: both sides with ``method='forward'``'s one; with
        ``'central'``, also one side when the other failed (the failing
        side may run at the refined step, and the running side's one-sided
        fallback needs agreeing firings)."""
        mine = [sides[j, s] for s in signs]
        ran = [side for side in mine if side.error is None]
        return (bool(ran) and not all(side.same for side in ran)
                and (len(ran) == len(mine) or method == 'central'))

    # Refinement: step / 10 while a counter differs at some output time.
    for _ in range(max_refinements):
        pending = [j for j in range(len(parameters)) if changed(j)]
        if not pending:
            break
        for j in pending:
            h[j] /= 10
            n_ref[j] += 1
        record(perturb([(j, s, h[j]) for j in pending for s in signs]))

    # Still differing: keep a one-sided difference on the side that agrees
    # with the nominal counters. Forward differences need backward runs for
    # that: the whole ladder of steps tried, in one batch.
    flagged = [j for j in range(len(parameters)) if changed(j)]
    if method == 'forward' and flagged:
        ladder = []
        for j in flagged:
            step = parameters[j].step
            for _ in range(n_ref[j] + 1):
                ladder.append((j, -1, step))
                step /= 10
        record(perturb(ladder))         # largest step first: best is kept

    jac = np.full((f0.size, len(parameters)), np.nan)
    rows = []
    for j, p in enumerate(parameters):
        plus = sides[j, 1]
        minus = sides.get((j, -1))
        error = plus.error or (minus.error if method == 'central' else None)
        scheme, step, side = method, h[j], None
        if (method == 'central' and (plus.error is None)
                != (minus.error is None)):
            # One side failed (e.g. a parameter at 0 run below zero): the
            # other side's one-sided difference at the same step, if its
            # firings agree with the nominal run (if not, changed() refined
            # the step and the column is in ``flagged``).
            sign, ok = (1, plus) if plus.error is None else (-1, minus)
            if ok.same:
                side, error = ok, None
                scheme = 'forward' if sign == 1 else 'backward'
        if j in flagged:
            if method == 'central' and plus.same != minus.same:
                chosen = 1 if plus.same else -1
            else:
                agree = [s for s in (1, -1) if (j, s) in best]
                chosen = agree[0] if len(agree) == 1 else None
            if chosen is not None:
                # That side ran (with agreeing firings) at that step: a
                # failure of the other side, or at another step, does not
                # void the column.
                side, step = best[j, chosen]
                scheme = 'forward' if chosen == 1 else 'backward'
                error = None
        if error is None:
            with np.errstate(all='ignore'):
                if side is not None:
                    jac[:, j] = ((side.f - f0) if scheme == 'forward'
                                 else (f0 - side.f)) / step
                elif scheme == 'central':
                    jac[:, j] = (plus.f - minus.f) / (2 * step)
                else:
                    jac[:, j] = (plus.f - f0) / step
        rows.append({'value': p.value, 'log': p.log, 'step': step,
                     'n_refinements': n_ref[j],
                     'event_count_changed': j in flagged,
                     'reliable': error is None and j not in flagged,
                     'scheme': scheme, 'error': error})
    return FDJacobian(jacobian=jac, f0=f0.copy(), parameters=parameters,
                      flags=flags_frame(rows, names), counters=counters,
                      conditions=list(nominal.conditions),
                      results=list(nominal.results),
                      tolerances=nominal.tolerances)
