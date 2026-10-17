# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Practical identifiability from local sensitivities.

:func:`identifiability` turns a :class:`LocalSensitivityResult` and a noise
model (absolute, or relative with a floor at a small fraction of each
output's peak) into the Fisher information matrix (FIM) of the measured
outputs and the diagnostics built on it: singular values, condition number,
numerical rank and null-space directions, Cramer-Rao standard errors, and
the Brun et al. (2001) importance (``delta_msqr``) and collinearity
(``gamma_K``) measures that pick an identifiable parameter subset.

The SVD core (:func:`_svd_analysis`) and the collinearity helpers take any
weighted Jacobian, so calibration (the covariance of a fit) reuses them.

References
----------
Brun, R., Reichert, P., Kuensch, H. R. (2001). Practical identifiability
analysis of large environmental simulation models. Water Resources
Research 37(4), 1015-1030.
"""
import itertools
import json
import math
import warnings
from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from ..engine._plotting import _set_font, _style_axes
from . import _fd
from ._common import (_CSV_COLUMNS, _EMPTY, _count, _finite, _id_list,
                      _read_long_csv, _real, _scatter)
from .sensitivity import LocalSensitivityResult

__all__ = ('IdentifiabilityResult', 'identifiability')

#: Default relative threshold: singular values ``<= rank_rtol * s_max``
#: count as zero.
RANK_RTOL = 1e-6
#: A parameter whose loading on the null space exceeds this has an infinite
#: Cramer-Rao bound (it moves along an unidentifiable direction). Parameters
#: on a null direction load O(0.1-1); finite-difference noise spills about
#: noise / gap (~1e-8 to 1e-4) onto the others, which must stay estimable.
NULL_LOADING = 1e-3
#: A relative Cramer-Rao standard error (``rel_se``; for a ``'log'``
#: parameter, the SE of ``ln p``) above this counts as practically
#: non-identifiable: a 95% interval of ``ln p`` wider than about +-2, i.e. a
#: factor of about 7 either way. Calibration's "relative SE > 1" flag uses
#: the same constant.
REL_SE_MAX = 1.0
#: Default lower bound of a relative noise model, as a fraction of each
#: output's peak: ``sigma_i = rel * max(|y_i|, RELATIVE_FLOOR * max_t |y|)``
#: (a detection limit; see :func:`identifiability`).
RELATIVE_FLOOR = 1e-3
#: The rank warning of :func:`identifiability` fires when the estimated
#: finite-difference error level (:func:`_fd.error_level`, scaled by the
#: realized weight amplification of a lowered ``relative_floor``) exceeds
#: ``rank_rtol`` and the smallest singular value the rank test retained,
#: relative to ``s_max``, is at most this many times that level (see
#: :func:`_rank_test_doubtful`). The level was measured at 1e-8 to 5e-7
#: on smooth toys at the defaults of local_sensitivity (the upper end
#: under a relative sigma) and 1e-6 to 2e-3 at CVODE's own rtol 1e-6.
#: Exact null directions were measured at 0.2-0.85 x the level (sum toys:
#: forward 0.21, CVODE rtol 0.76, central rel_step 3e-3 0.85), 3 x on
#: Arrhenius ``A``/``Ea`` and 24 x on ``(sqrt(a - 1) + b) * S`` at
#: ``a = 1.01``, ``rel_step=2e-3``; well-posed decay toys sit at 245 x or
#: more (forward) and 4e3 x (``rel_step=1e-2``). On the depleted sum toy
#: of concepts.rst (``(k1 + k2) * S`` to ``t = 30``, ``('relative',
#: 0.1)``), floors of 1e-4 and 1e-5 lift the null to 1.9e-6 and 1.3e-5 of
#: ``s_max`` (full rank), and ``rel = 0.01`` at 1e-7 reads both
#: identifiable; the scaled level flags all three. The level assumes
#: smooth dependence on ``ln p``: at a depletion front the truncation
#: error is far larger (the ``m1``/``m2`` null of Monod
#: ``(m1 + m2) * S / (Ks + S) * X`` sits at 3.8e-4 of ``s_max``, 630 x the
#: level at the defaults, with ``gamma`` about 490), so a ratio above the
#: margin is no proof of identifiability; the step-convergence check of
#: :func:`identifiability`'s Notes is.
FD_ERROR_MARGIN = 100
_SPACES = ('log', 'linear')
_PLOTS = ('delta_msqr', 'singular_values')
_CSV_FORMAT = 'IdentifiabilityResult'
#: Subsets per vectorized eigenvalue batch in the collinearity search.
_CHUNK = 10000
_FLOAT_MAX = float(np.finfo(float).max)


# --- numerical cores (private; also used by calibration) ---------------------
@dataclass(frozen=True, eq=False)
class _SVD:
    """SVD summary of a weighted Jacobian ``J`` (see :func:`_svd_analysis`).

    Attributes
    ----------
    singular_values : numpy.ndarray
        ``(n_p,)``, descending, padded with zeros when ``J`` has fewer rows
        than columns. The FIM ``J^T J`` has eigenvalues
        ``singular_values ** 2``.
    vt : numpy.ndarray
        ``(n_p, n_p)``; row ``k`` is the right singular vector of
        ``singular_values[k]``, signed so its largest-magnitude entry is
        positive.
    numerical_rank : int
        Number of singular values above ``rank_rtol * s_max``.
    condition_number : float
        ``s_max / s_min`` of ``J`` (the FIM's is its square; ``inf`` if
        ``s_min`` is 0 or ``J`` is 0).
    covariance : numpy.ndarray
        ``(n_p, n_p)``: the pseudo-inverse of ``J^T J`` restricted to the
        numerical range (the Cramer-Rao covariance where it exists).
    se : numpy.ndarray
        ``(n_p,)``: ``sqrt(diag(covariance))``, ``inf`` for a parameter with
        a loading above ``NULL_LOADING`` (``1e-3``) on the null space.
    null_loading : numpy.ndarray
        ``(n_p,)``: norm of each parameter's projection onto the null space
        (0 for a full-rank ``J``).
    """
    singular_values: np.ndarray
    vt: np.ndarray
    numerical_rank: int
    condition_number: float
    covariance: np.ndarray
    se: np.ndarray
    null_loading: np.ndarray


def _svd_analysis(jac, *, rank_rtol=RANK_RTOL):
    """Rank, null space and Cramer-Rao covariance of ``J^T J``.

    Parameters
    ----------
    jac : array_like
        ``(n_rows, n_p)`` weighted Jacobian (residuals or outputs divided by
        their standard deviations, differentiated in each parameter's own
        space); finite.
    rank_rtol : float
        Relative threshold for the numerical rank. It assumes derivatives
        accurate to well below it, as central differences at a small step
        and tight tolerances usually are (see :func:`identifiability`'s
        Notes); a one-sided finite-difference Jacobian (e.g. near a bound)
        is not, so check ``_fd.error_level(flags, tolerances).total``
        (noise plus truncation, as identifiability does; scaled by the
        realized weight amplification of the rows, if any) against it
        and the smallest retained singular value, with
        :func:`_rank_test_doubtful`.

    Returns
    -------
    _SVD

    Raises
    ------
    ValueError
        If ``jac`` is not 2-D with at least one column, or not finite.
    """
    jac = np.asarray(jac, dtype=float)
    if jac.ndim != 2 or jac.shape[1] == 0:
        raise ValueError(f'The Jacobian must be 2-D with at least one column; '
                         f'got shape {jac.shape}.')
    if not np.all(np.isfinite(jac)):
        raise ValueError('The Jacobian must be finite.')
    n_rows, n_p = jac.shape
    s = np.zeros(n_p)
    if n_rows:
        # Full V only when it is needed (fewer rows than columns); U is never
        # used, so it is kept small.
        _, sv, vt = np.linalg.svd(jac, full_matrices=n_rows < n_p)
        s[:sv.size] = sv
    else:
        vt = np.eye(n_p)
    pivot = np.argmax(np.abs(vt), axis=1)
    signs = np.sign(vt[np.arange(n_p), pivot])
    signs[signs == 0] = 1.0
    vt = vt * signs[:, None]
    s_max = s[0]
    keep = (s > rank_rtol * s_max) & (s > 0)
    rank = int(keep.sum())
    with np.errstate(all='ignore'):
        condition = float(s_max / s[-1]) if s[-1] > 0 else math.inf
        v = vt.T
        kept = v[:, keep]
        covariance = (kept / s[keep] ** 2) @ kept.T
        null_loading = np.sqrt(np.sum(v[:, ~keep] ** 2, axis=1))
        se = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    se = np.where(null_loading > NULL_LOADING, math.inf, se)
    return _SVD(singular_values=s, vt=vt, numerical_rank=rank,
                condition_number=condition, covariance=covariance, se=se,
                null_loading=null_loading)


def _rel_se(values, space, se):
    """``(se_natural, rel_se)`` from standard errors in each parameter's
    space.

    ``se_natural`` is ``|p| * se`` for a ``'log'`` column (first order;
    ``p > 0`` there) and ``se`` for a ``'linear'`` one;
    ``rel_se = se_natural / |p|``, ``NaN`` at ``p = 0``.
    """
    values = np.asarray(values, dtype=float)
    se = np.asarray(se, dtype=float)
    log = np.array([s == 'log' for s in space], dtype=bool)
    with np.errstate(all='ignore'):
        natural = np.where(log, np.abs(values) * se, se)
        rel = np.where(values != 0, natural / np.abs(values), np.nan)
    return natural, rel


def _identifiable(values, space, se):
    """The practical verdict per parameter: a finite Cramer-Rao ``se`` (no
    loading on a null direction) and ``rel_se <= REL_SE_MAX`` (a
    ``'linear'`` parameter at ``p = 0``, whose ``rel_se`` is ``NaN``,
    passes that check). The one definition behind
    :attr:`IdentifiabilityResult.table` and
    :meth:`IdentifiabilityResult.best_subset`."""
    se = np.asarray(se, dtype=float)
    _, rel = _rel_se(values, space, se)
    with np.errstate(all='ignore'):
        return np.isfinite(se) & (np.isnan(rel) | (rel <= REL_SE_MAX))


def _relative_sigma(y, rel, peak, floor=RELATIVE_FLOOR):
    """Bounded relative noise ``rel * max(|y|, floor * peak)``: a detection
    limit at ``floor`` times the output's ``peak`` (``max_t |y|``).
    Overflow gives ``inf`` (no warning): check the result is finite and
    positive.

    For calibration (proportional Observable noise), the sigma inside
    least squares must not depend on the parameters being fitted: take
    ``y`` and ``peak`` from the data (which over-weights low observations,
    noticeably at a large ``rel``), or from the prediction at a previous
    iterate, held fixed for each solve (IRLS style). A sigma from the
    prediction inside ``residuals(theta)`` biases plain ``least_squares``
    toward larger predictions (it needs a likelihood with a
    ``sum(log(sigma))`` term), and its Jacobian is not the Gauss-Newton
    FIM that :func:`_svd_analysis` assumes. Expose a floor as
    ``relative_floor`` (default ``RELATIVE_FLOOR``) so both features share
    one name; a combined model's absolute term needs no floor."""
    with np.errstate(all='ignore'):
        return rel * np.maximum(np.abs(np.asarray(y, dtype=float)),
                                floor * peak)


def _rank_test_doubtful(level, rank_rtol, singular_values, numerical_rank):
    """Whether the estimate of the finite-difference error casts doubt on
    the rank test: the error ``level`` (relative to ``s_max``: e.g.
    ``_fd.error_level(flags, tolerances).total``, scaled by the realized
    weight amplification of a lowered ``relative_floor`` as
    :func:`identifiability` records it) exceeds ``rank_rtol``, and the
    smallest singular value the rank test retained,
    ``singular_values[numerical_rank - 1] / singular_values[0]``, is at
    most ``FD_ERROR_MARGIN`` times it, so it could be a null direction
    lifted by that error.

    Singular values already counted as zero (``se = inf``; exact zeros of
    a parameter without effect, or the zero padding of a Jacobian with
    fewer rows than columns, among them) are the conservative verdict and
    raise no doubt; nor do fewer than two parameters or an all-zero
    Jacobian. ``False`` is no proof: the estimate misses depletion fronts,
    switches and a loose ``atol``. The one gate behind the rank warning
    and :attr:`IdentifiabilityResult.rank_test_doubtful`."""
    s = np.asarray(singular_values, dtype=float)
    if len(s) < 2 or not numerical_rank > 0:
        return False
    if not level > rank_rtol:                    # also for NaN
        return False
    with np.errstate(all='ignore'):
        ratio = s[numerical_rank - 1] / s[0]     # s[0] > 0 at rank >= 1
    return bool(ratio <= FD_ERROR_MARGIN * level)


def _amplified(level, amp):
    """The error ``level`` under row weights amplified up to ``amp`` times
    (``amp >= 1``, see :func:`identifiability`): ``level * amp``, clipped
    to the largest float so that a vanishing ``relative_floor`` cannot
    overflow it; ``level`` itself for ``amp <= 1`` or ``level <= 0``."""
    if not (amp > 1 and level > 0):
        return level
    return min(level * amp, _FLOAT_MAX)


def _brun(jac):
    """``(delta_msqr, column norms, normalized Gram matrix)`` of ``jac``.

    ``delta_msqr_j = ||s_j|| / sqrt(n_rows)`` (Brun et al. 2001); the
    normalized Gram matrix is ``S~^T S~`` with unit-length columns ``S~``
    (rows/columns of zero-norm parameters are ``NaN``).
    """
    jac = np.asarray(jac, dtype=float)
    n_rows = jac.shape[0]
    norms = np.linalg.norm(jac, axis=0)
    with np.errstate(all='ignore'):
        delta = norms / math.sqrt(n_rows) if n_rows else np.zeros_like(norms)
        gram = jac.T @ jac
        scale = np.where(norms > 0, norms, np.nan)
        normalized = gram / np.outer(scale, scale)
    return delta, norms, normalized


def _gamma(normalized, zero, subsets):
    """Collinearity indices ``gamma_K = 1 / sqrt(min eig(S~_K^T S~_K))`` of
    the index ``subsets`` (``(n_sub, k)`` int array); ``inf`` for a subset
    holding a zero-norm parameter or a (numerically) singular one."""
    subsets = np.asarray(subsets, dtype=int)
    n_sub, k = subsets.shape
    out = np.empty(n_sub)
    for start in range(0, n_sub, _CHUNK):
        part = subsets[start:start + _CHUNK]
        bad = zero[part].any(axis=1)
        sub = normalized[part[:, :, None], part[:, None, :]]
        sub[bad] = np.eye(k)
        lam = np.linalg.eigvalsh(sub)[:, 0]
        with np.errstate(all='ignore'):
            g = 1.0 / np.sqrt(np.clip(lam, 0.0, None))
        g[bad] = math.inf
        out[start:start + len(part)] = g
    return out


def _rho(gram, n_rows, subsets):
    """Brun's determinant measure ``rho_K = det(S_K^T S_K / n)^(1/2k)``
    (so ``rho`` of one parameter is its ``delta_msqr``)."""
    subsets = np.asarray(subsets, dtype=int)
    n_sub, k = subsets.shape
    out = np.empty(n_sub)
    scale = 1.0 / n_rows if n_rows else 0.0
    for start in range(0, n_sub, _CHUNK):
        part = subsets[start:start + _CHUNK]
        sub = gram[part[:, :, None], part[:, None, :]] * scale
        lam = np.clip(np.linalg.eigvalsh(sub), 0.0, None)
        with np.errstate(all='ignore'):
            out[start:start + len(part)] = np.exp(
                np.sum(np.log(lam), axis=1) / (2 * k))
    return out


def _greedy_collinear(normalized, zero, k):
    """The most collinear size-``k`` subset grown from each parameter: start
    from it and repeatedly add the parameter that maximizes ``gamma`` (ties:
    lowest index). Parameters without any effect (``zero``), which make
    every subset infinitely collinear, are added only when no other is left.
    Distinct subsets, as sorted index tuples, in seed order."""
    n_p = normalized.shape[0]
    found = {}
    for seed in range(n_p):
        subset = [seed]
        while len(subset) < k:
            rest = [c for c in range(n_p) if c not in subset]
            candidates = [c for c in rest if not zero[c]] or rest
            subs = np.array([sorted([*subset, c]) for c in candidates])
            g = _gamma(normalized, zero, subs)
            subset = list(subs[int(np.argmax(g))])
        found.setdefault(tuple(int(i) for i in subset), None)
    return list(found)


# --- result ------------------------------------------------------------------
def _check_rank_rtol(rank_rtol):
    rank_rtol = _finite('rank_rtol', rank_rtol)
    if not 0 < rank_rtol < 1:
        raise ValueError(f'rank_rtol must be in (0, 1); got {rank_rtol!r}.')
    return rank_rtol


def _check_relative_floor(relative_floor):
    relative_floor = _finite('relative_floor', relative_floor)
    if not 0 < relative_floor <= 1:
        raise ValueError(f'relative_floor must be in (0, 1]; got '
                         f'{relative_floor!r}.')
    return relative_floor


def _check_fd_error_level(level):
    level = _finite('fd_error_level', level)
    if level < 0:
        raise ValueError(f'fd_error_level must be >= 0; got {level!r}.')
    return level


def _sigma_to_json(sigma):
    """JSON-safe form of a sigma record: ``{'relative': rel}`` for
    ``('relative', rel)``, ``{'mapping': {output: entry}}`` for a mapping,
    with a relative entry ``('relative', r)`` as ``{'relative': r}``."""
    if isinstance(sigma, tuple):
        return {'relative': sigma[1]}
    if isinstance(sigma, Mapping):
        return {'mapping': {out: ({'relative': v[1]} if isinstance(v, tuple)
                                  else v) for out, v in sigma.items()}}
    return {'absolute': sigma}


def _sigma_from_json(data):
    if data is None:
        return None
    if 'relative' in data:
        return ('relative', data['relative'])
    if 'mapping' in data:
        return {out: (('relative', v['relative']) if isinstance(v, dict)
                      else v) for out, v in data['mapping'].items()}
    return data['absolute']


@dataclass(eq=False, repr=False)
class IdentifiabilityResult:
    """Outcome of :func:`identifiability`.

    Constructed from a weighted, log-scaled sensitivity matrix; every other
    attribute is derived from it (so :meth:`from_csv` restores it exactly).

    Parameters
    ----------
    parameters : list of str
        Parameter ids, in column order.
    values : numpy.ndarray
        Nominal values ``p``, shape ``(n_p,)``, model units.
    space : list of {'log', 'linear'}
        Space of each column: ``'log'`` (``x = ln p``) for ``p > 0``,
        ``'linear'`` (``x = p``) otherwise.
    sensitivity : numpy.ndarray
        ``S``, shape ``(n_rows, n_p)``: ``(dy_i / dx_j) / sigma_i`` for each
        measured row ``i`` (an output at a time). Dimensionless for ``'log'``
        columns.
    row_outputs : list of str
        The output of each row.
    row_times : numpy.ndarray
        The time of each row (model time units).
    reliable : numpy.ndarray
        ``(n_p,)`` bool: the sensitivity column passed the finite-difference
        checks (see :attr:`LocalSensitivityResult.flags`).
    sigma : dict or tuple, optional
        The noise model used, for the record: ``{output: entry}`` with an
        absolute sigma (a number, or a list of one value per row time in
        increasing time) or ``('relative', rel)`` per entry; or
        ``('relative', rel)`` for every output.
    rank_rtol : float
        Singular values ``<= rank_rtol * s_max`` count as zero.
    relative_floor : float, optional
        Detection limit of the relative noise model, as a fraction of each
        output's peak: ``|y|`` counts as at least
        ``relative_floor * max_t |y|``, so
        ``sigma >= rel * relative_floor * max_t |y|`` (see
        :func:`identifiability`), for the record; in ``(0, 1]``. ``None``
        (default) when no measured output has relative-sigma rows (all are
        absolute, or every relative output was left out as always 0), or
        for a matrix built otherwise.
    fd_error_level : float, optional
        The estimated finite-difference error level of the sensitivity
        columns relative to the largest singular value (integration noise
        plus truncation error, estimated from ``sens.tolerances`` and
        ``sens.flags``; under a ``relative_floor`` below the default,
        scaled by the realized weight amplification of the rows, see
        :ref:`identifiability-noise`); ``>= 0``. It sets
        :attr:`rank_test_doubtful`. ``None`` (default) when unknown: for a
        matrix built otherwise, or from :func:`identifiability` when the
        sensitivities were computed with an ``atol`` above 1e-8, which the
        estimate does not model, or with a relative step so small that
        the estimate overflows.
    units : dict
        The model's units.

    Attributes
    ----------
    fim : pandas.DataFrame
        Fisher information ``S^T S`` (``n_p x n_p``) in each parameter's own
        space.
    singular_values : numpy.ndarray
        Singular values of ``S`` (descending, ``n_p`` of them, zero padded);
        the FIM's eigenvalues are their squares.
    condition_number : float
        ``s_max / s_min`` of ``S`` (the FIM's is its square); ``inf`` if
        ``s_min`` is 0.
    numerical_rank : int
        Numerical rank of ``S``: the number of singular values above
        ``rank_rtol * s_max``.
    directions : pandas.DataFrame
        Right singular vectors (index ``parameter``, columns ``v1 ... vn`` in
        the order of :attr:`singular_values`), each signed so its largest
        loading is positive.
    null_directions : pandas.DataFrame
        The columns of :attr:`directions` beyond :attr:`numerical_rank`:
        parameter combinations (in ``ln p`` for ``'log'`` columns) the data
        cannot see, e.g. ``(1, -1) / sqrt(2)`` when only ``k1 * k2``
        matters.
    covariance : pandas.DataFrame
        Cramer-Rao covariance: the pseudo-inverse of the FIM on its
        numerical range (the inverse at full rank).
    se : pandas.Series
        Cramer-Rao standard errors in each parameter's own space (for
        ``'log'``: of ``ln p``, i.e. relative); ``inf`` for a parameter on a
        null direction.
    delta_msqr : pandas.Series
        Brun et al. (2001) importance ``sqrt(mean_i S_ij^2)``.
    rank_test_doubtful : bool or None
        ``True`` when :attr:`fd_error_level` exceeds ``rank_rtol`` and the
        smallest retained singular value is within 100 times it relative
        to the largest, so it may be a null direction lifted by
        finite-difference error (:func:`identifiability` warns). ``False``
        when no doubt was raised (always for one parameter), which is no
        proof; ``None`` when :attr:`fd_error_level` is unknown, for two or
        more parameters. See :ref:`identifiability-noise`.

    Raises
    ------
    ValueError
        For duplicate parameters; a ``space`` that does not hold one of
        ``'log'``, ``'linear'`` per parameter; inconsistent shapes of
        ``sensitivity``, ``values``, ``reliable``, ``row_outputs`` and
        ``row_times``; a ``sensitivity`` matrix that is not finite or has no
        column; non-finite ``values``, or a ``'log'`` column at
        ``p <= 0``; ``rank_rtol`` outside ``(0, 1)``; ``relative_floor``
        outside ``(0, 1]``; or a non-finite or negative ``fd_error_level``.
    TypeError
        If ``rank_rtol``, ``relative_floor`` or ``fd_error_level`` is not a
        real number.
    """
    parameters: list
    values: np.ndarray
    space: list
    sensitivity: np.ndarray
    row_outputs: list
    row_times: np.ndarray
    reliable: np.ndarray
    sigma: object = None
    rank_rtol: float = RANK_RTOL
    relative_floor: float | None = None
    fd_error_level: float | None = None
    units: dict = field(default_factory=dict)

    def __post_init__(self):
        self.parameters = list(self.parameters)
        self.space = list(self.space)
        # Copies, so later changes to the caller's arrays cannot leave the
        # derived attributes stale. C order whatever the input's layout, so
        # the derived attributes (whose reductions depend on it in the last
        # bit) are the same for the same values, e.g. after a CSV round
        # trip.
        self.values = np.array(self.values, dtype=float)
        self.sensitivity = np.array(self.sensitivity, dtype=float, order='C')
        self.row_outputs = list(self.row_outputs)
        self.row_times = np.array(self.row_times, dtype=float)
        self.reliable = np.array(self.reliable, dtype=bool)
        self.units = dict(self.units)
        self.rank_rtol = _check_rank_rtol(self.rank_rtol)
        if self.relative_floor is not None:
            self.relative_floor = _check_relative_floor(self.relative_floor)
        if self.fd_error_level is not None:
            self.fd_error_level = _check_fd_error_level(self.fd_error_level)
        n_p = len(self.parameters)
        if len(set(self.parameters)) != n_p:
            raise ValueError(f'Duplicate parameters in {self.parameters}.')
        bad_space = [s for s in self.space if s not in _SPACES]
        if bad_space or len(self.space) != n_p:
            raise ValueError(f'space must hold one of {_SPACES} per '
                             'parameter.')
        n_rows = len(self.row_outputs)
        if (self.sensitivity.shape != (n_rows, n_p)
                or self.values.shape != (n_p,)
                or self.reliable.shape != (n_p,)
                or self.row_times.shape != (n_rows,)):
            raise ValueError('Inconsistent shapes: sensitivity must be '
                             '(n_rows, n_p) with one row output/time per row '
                             'and one value/reliable flag per parameter.')
        if not np.all(np.isfinite(self.values)):
            raise ValueError(f'values must be finite; got '
                             f'{self.values.tolist()}.')
        bad = [p for p, v, sp in zip(self.parameters, self.values, self.space)
               if sp == 'log' and not v > 0]
        if bad:
            raise ValueError(f"'log' columns need p > 0 (x = ln p); {bad} "
                             'are not.')
        svd = _svd_analysis(self.sensitivity, rank_rtol=self.rank_rtol)
        self._svd = svd
        delta, norms, normalized = _brun(self.sensitivity)
        self._norms, self._normalized = norms, normalized
        index = pd.Index(self.parameters, name='parameter', dtype=object)
        self.fim = pd.DataFrame(self.sensitivity.T @ self.sensitivity,
                                index=index, columns=list(self.parameters))
        self.singular_values = svd.singular_values.copy()
        self.condition_number = svd.condition_number
        self.numerical_rank = svd.numerical_rank
        names = [f'v{k + 1}' for k in range(n_p)]
        self.directions = pd.DataFrame(svd.vt.T.copy(), index=index,
                                       columns=names)
        self.null_directions = self.directions.iloc[
            :, self.numerical_rank:].copy()
        self.covariance = pd.DataFrame(svd.covariance.copy(), index=index,
                                       columns=list(self.parameters))
        self.se = pd.Series(svd.se.copy(), index=index, name='se')
        self.delta_msqr = pd.Series(delta, index=index, name='delta_msqr')
        # One column has no compensating combination: never doubtful.
        self.rank_test_doubtful = (
            False if n_p < 2 else None if self.fd_error_level is None
            else _rank_test_doubtful(self.fd_error_level, self.rank_rtol,
                                     svd.singular_values,
                                     svd.numerical_rank))

    def __repr__(self):
        doubt = ', rank_test=doubtful' if self.rank_test_doubtful else ''
        return (f'IdentifiabilityResult(parameters={self.parameters}, '
                f'numerical_rank={self.numerical_rank}/'
                f'{len(self.parameters)}, '
                f'condition_number={self.condition_number:.3g}, '
                f'unreliable={int((~self.reliable).sum())}{doubt})')

    # --- views ---------------------------------------------------------------
    @property
    def correlation(self):
        """Correlation matrix from :attr:`covariance` (``pandas.DataFrame``);
        ``NaN`` in the rows and columns of parameters with infinite or zero
        :attr:`se`."""
        se = self.se.to_numpy()
        ok = np.isfinite(se) & (se > 0)
        with np.errstate(all='ignore'):
            scale = np.where(ok, se, np.nan)
            corr = self.covariance.to_numpy() / np.outer(scale, scale)
        return pd.DataFrame(corr, index=self.covariance.index,
                            columns=list(self.parameters))

    @property
    def table(self):
        """One row per parameter (``pandas.DataFrame``, index ``parameter``).

        Columns: ``value`` (nominal ``p``), ``space``, ``delta_msqr``,
        ``se`` (Cramer-Rao standard error in the parameter's space),
        ``se_natural`` (first-order SE of ``p``: ``p * se`` for ``'log'``),
        ``rel_se`` (``se_natural / |p|``; ``NaN`` at ``p = 0``),
        ``null_loading`` (norm of the projection onto the null space),
        ``null_direction`` (``null_loading`` above ``1e-3``: the parameter
        moves along a direction the data cannot see, so ``se`` is
        infinite; finite-difference noise leaves small loadings on the
        others), ``identifiable`` and ``reliable``.

        ``identifiable`` is the practical verdict: not on a null direction
        *and* ``rel_se`` at most 1 (i.e. a 95% interval of ``ln p``
        narrower than about +-2, a factor of about 7 either way). A
        parameter off every null direction can still be too weakly
        determined to estimate -- e.g. ``Vmax`` and ``Km`` of
        Michaelis-Menten kinetics in its linear regime, each with a finite
        ``se`` above 10 -- and then reads ``False``. A ``'linear'``
        parameter at ``p = 0`` has no relative scale (``rel_se`` is
        ``NaN``) and passes the precision check. Its column, like that of
        any ``'linear'`` parameter, is per model unit of ``p``, and that
        unit sways the ``se``, ``null_direction`` and ``identifiable``
        entries of the other parameters too (see :func:`identifiability`).
        """
        se = self.se.to_numpy()
        null_loading = self._svd.null_loading.copy()
        natural, rel = _rel_se(self.values, self.space, se)
        return pd.DataFrame({
            'value': self.values.copy(),
            'space': pd.Series(self.space, dtype=object).to_numpy(),
            'delta_msqr': self.delta_msqr.to_numpy(),
            'se': se.copy(), 'se_natural': natural, 'rel_se': rel,
            'null_loading': null_loading,
            'null_direction': null_loading > NULL_LOADING,
            'identifiable': _identifiable(self.values, self.space, se),
            'reliable': self.reliable.copy()}, index=self.se.index.copy())

    def _indices(self, names):
        return [self.parameters.index(n) for n in names]

    def collinearity(self, subset_size=2, max_subsets=20000):
        """Collinearity index of parameter subsets (Brun et al. 2001).

        ``gamma_K = 1 / sqrt(min eig(S~_K^T S~_K))`` with the columns of
        ``S`` scaled to unit length: 1 for orthogonal sensitivities, large
        when a change of one parameter can be compensated by the others
        (``gamma > 10-20`` is usually taken as poorly identifiable;
        ``inf`` for an exactly singular subset or one holding a parameter
        without any effect).

        Parameters
        ----------
        subset_size : int
            Size ``k`` of the subsets (``1 <= k <= n_p``; default 2).
        max_subsets : int
            Evaluate all ``C(n_p, k)`` subsets if there are at most this many
            (default 20000); otherwise search greedily (see Notes).

        Returns
        -------
        pandas.DataFrame
            Columns ``parameters`` (tuple of str, in parameter order),
            ``gamma`` and ``rho`` (Brun's determinant measure
            ``det(S_K^T S_K / n_rows)^(1/2k)``), sorted by decreasing
            ``gamma`` (ties in enumeration order). ``attrs['search']`` is
            ``'exhaustive'`` or ``'greedy'``.

        Raises
        ------
        TypeError, ValueError
            For a non-integer or out-of-range ``subset_size`` or
            ``max_subsets``.

        Notes
        -----
        The greedy search grows one subset from each parameter by repeatedly
        adding the parameter that makes it most collinear (parameters
        without any effect, which make every subset infinitely collinear,
        only when no other is left), so it lists at most ``n_p`` subsets; it
        is not exhaustive.
        """
        n_p = len(self.parameters)
        k = _count('subset_size', subset_size, 1)
        if k > n_p:
            raise ValueError(f'subset_size must be at most the number of '
                             f'parameters ({n_p}); got {k}.')
        cap = _count('max_subsets', max_subsets, 1)
        zero = ~(self._norms > 0)
        if math.comb(n_p, k) <= cap:
            subsets = np.array(list(itertools.combinations(range(n_p), k)),
                               dtype=int).reshape(-1, k)
            search = 'exhaustive'
        else:
            subsets = np.array(_greedy_collinear(self._normalized, zero, k),
                               dtype=int).reshape(-1, k)
            search = 'greedy'
        gamma = _gamma(self._normalized, zero, subsets)
        rho = _rho(self.fim.to_numpy(), len(self.row_outputs), subsets)
        order = np.argsort(-gamma, kind='stable')
        names = [tuple(self.parameters[i] for i in row) for row in subsets]
        out = pd.DataFrame({
            'parameters': pd.Series([names[i] for i in order], dtype=object),
            'gamma': gamma[order], 'rho': rho[order]})
        out.attrs['search'] = search
        return out

    def gamma(self, names):
        """Collinearity index ``gamma_K`` of one subset.

        Parameters
        ----------
        names : sequence of str
            Parameter ids (at least one).

        Returns
        -------
        float

        Raises
        ------
        ValueError
            For unknown, repeated or no names.
        """
        names = _id_list('names', names)
        unknown = [n for n in names if n not in self.parameters]
        if unknown:
            raise ValueError(f'Unknown parameters {unknown}; parameters: '
                             f'{self.parameters}.')
        idx = np.array([sorted(self._indices(names))], dtype=int)
        return float(_gamma(self._normalized, ~(self._norms > 0), idx)[0])

    def best_subset(self, gamma_max=10.0):
        """A large identifiable subset: greedy by importance under a
        collinearity cap.

        Parameters are taken in order of decreasing :attr:`delta_msqr`
        (ties in parameter order), and each is kept if the kept set plus it
        meets both criteria:

        1. its collinearity index (:meth:`gamma`) is at most ``gamma_max``;
        2. analyzed on its own (its columns of :attr:`sensitivity`, the
           other parameters fixed), every member reads ``identifiable``
           as in :attr:`table`: a finite Cramer-Rao standard error (no
           loading on a null direction, so the set has full numerical
           rank) and a relative standard error of at most 1.

        Parameters
        ----------
        gamma_max : float
            Collinearity cap (``>= 1``; default 10, a common threshold).
            ``inf`` removes the cap; the subset is still identifiable by
            criterion 2.

        Returns
        -------
        list of str
            The kept parameters, most important first.

        Raises
        ------
        TypeError, ValueError
            If ``gamma_max`` is not a real number ``>= 1``.

        Notes
        -----
        Of two parameters that only act through their product (or any
        other compensating combination) only the first is kept, and a
        parameter without any effect, or whose column is at the
        finite-difference noise level, is never kept, however orthogonal
        its (noise) direction: its relative standard error is huge or
        infinite. Criterion 2 fixes the parameters outside the subset at
        their nominal values, so a kept parameter can read
        ``identifiable=False`` in :attr:`table`, where every parameter is
        free (e.g. ``Vmax`` of Michaelis-Menten kinetics in its linear
        regime, which is well determined once ``Km`` is fixed). The search
        is greedy, so it is fast for many parameters but not guaranteed to
        find the largest admissible subset. Unreliable columns (see
        :attr:`reliable`) take part like the others; check them first.
        """
        gamma_max = _real('gamma_max', gamma_max)
        if not gamma_max >= 1:
            raise ValueError(f'gamma_max must be >= 1; got {gamma_max!r}.')
        zero = ~(self._norms > 0)
        chosen = []
        for j in np.argsort(-self.delta_msqr.to_numpy(), kind='stable'):
            idx = sorted([*chosen, int(j)])
            gamma = _gamma(self._normalized, zero, np.array([idx]))[0]
            if gamma > gamma_max:
                continue
            svd = _svd_analysis(self.sensitivity[:, idx],
                                rank_rtol=self.rank_rtol)
            if _identifiable(self.values[idx], [self.space[i] for i in idx],
                             svd.se).all():
                chosen.append(int(j))
        return [self.parameters[j] for j in chosen]

    # --- plots ---------------------------------------------------------------
    def plot(self, kind='delta_msqr', ax=None):
        """Plot parameter importance or the singular values.

        Parameters
        ----------
        kind : {'delta_msqr', 'singular_values'}
            ``'delta_msqr'``: horizontal bars of :attr:`delta_msqr`, most
            important on top (unreliable parameters marked ``*``).
            ``'singular_values'``: :attr:`singular_values` on a log axis with
            the rank threshold ``rank_rtol * s_max`` dashed (zeros are not
            drawn). Default ``'delta_msqr'``.
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure otherwise.

        Returns
        -------
        (matplotlib.figure.Figure, matplotlib.axes.Axes)

        Raises
        ------
        ValueError
            For an unknown ``kind``.
        """
        if kind not in _PLOTS:
            raise ValueError(f'kind must be one of {_PLOTS}; got {kind!r}.')
        _set_font()
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        if kind == 'delta_msqr':
            delta = self.delta_msqr.to_numpy()
            order = np.argsort(-delta, kind='stable')[::-1]
            labels = [p if ok else f'{p}*' for p, ok
                      in zip(self.parameters, self.reliable)]
            ax.barh(range(len(order)), delta[order], color='steelblue')
            ax.set_yticks(range(len(order)))
            ax.set_yticklabels([labels[j] for j in order])
            ax.set_xlabel('delta_msqr')
        else:
            s = self.singular_values
            pos = s > 0
            x = np.arange(1, len(s) + 1)
            ax.semilogy(x[pos], s[pos], 'o-', color='steelblue')
            if pos.any():
                ax.axhline(self.rank_rtol * s[0], color='gray',
                           linestyle='--', linewidth=0.8)
            ax.set_xticks(x)
            ax.set_xlabel('index')
            ax.set_ylabel('singular value')
        _style_axes(ax)
        ax.set_title(f'rank {self.numerical_rank}/{len(self.parameters)}, '
                     'condition '
                     f'number {self.condition_number:.3g}')
        return fig, ax

    # --- IO ------------------------------------------------------------------
    def to_csv(self, path):
        """Write the result to a long-format CSV.

        The file has the columns ``quantity, key, parameter, index, value``:
        ``meta`` rows (``key``: field name), ``parameter_value`` and
        ``reliable`` (per ``parameter``), ``row`` (``index``: row,
        ``key``: output, ``value``: time) and ``sensitivity`` (``index``:
        row, ``parameter``). Values are JSON, so :meth:`from_csv` restores
        the matrix exactly and re-derives every other attribute from it.

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
            If ``sigma`` or the units are not JSON-safe.
        """
        sigma = None if self.sigma is None else _sigma_to_json(self.sigma)
        meta = {'nskinetics': _CSV_FORMAT, 'parameters': self.parameters,
                'space': self.space, 'sigma': sigma,
                'rank_rtol': self.rank_rtol,
                'relative_floor': self.relative_floor,
                'fd_error_level': self.fd_error_level, 'units': self.units,
                'n_rows': len(self.row_outputs)}
        dumps = json.dumps
        rows = [('meta', key, '', '', dumps(value, allow_nan=False))
                for key, value in meta.items()]
        for j, par in enumerate(self.parameters):
            rows.append(('parameter_value', '', par, '',
                         dumps(float(self.values[j]))))
            rows.append(('reliable', '', par, '',
                         dumps(bool(self.reliable[j]))))
        for i, (out, t) in enumerate(zip(self.row_outputs, self.row_times)):
            rows.append(('row', out, '', i, dumps(float(t))))
        for j, par in enumerate(self.parameters):
            rows += [('sensitivity', '', par, i, dumps(float(v)))
                     for i, v in enumerate(self.sensitivity[:, j])]
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
        IdentifiabilityResult
            With the same matrix, names, values and settings (exactly), and
            the derived attributes recomputed from them.

        Raises
        ------
        ValueError
            If the file was not written by :meth:`to_csv`, or is incomplete
            (a missing ``meta`` key, or a parameter value, flag, row or
            matrix entry missing or duplicated, e.g. a truncated or edited
            file).
        """
        meta, groups = _read_long_csv(path, _CSV_FORMAT, required=(
            'parameters', 'space', 'sigma', 'rank_rtol', 'relative_floor',
            'fd_error_level', 'units', 'n_rows'))
        params, n_rows = meta['parameters'], meta['n_rows']
        n_p = len(params)
        par_idx = {p: j for j, p in enumerate(params)}

        def scatter(quantity, shape, *axes, **kwargs):
            return _scatter(path, _CSV_FORMAT, groups, quantity, shape, axes,
                            **kwargs)

        values = scatter('parameter_value', (n_p,), ('parameter', par_idx))
        reliable = scatter('reliable', (n_p,), ('parameter', par_idx),
                           parse=bool, dtype=bool)
        row_times = scatter('row', (n_rows,), ('index', None))
        row_outputs = [''] * n_rows      # the rows cover 0..n_rows-1 once
        g = groups.get('row', _EMPTY)
        for i, out in zip(g['index'], g['key']):
            row_outputs[int(i)] = out
        sensitivity = scatter('sensitivity', (n_rows, n_p), ('index', None),
                              ('parameter', par_idx))
        return cls(parameters=params, values=values, space=meta['space'],
                   sensitivity=sensitivity, row_outputs=row_outputs,
                   row_times=row_times, reliable=reliable,
                   sigma=_sigma_from_json(meta['sigma']),
                   rank_rtol=meta['rank_rtol'],
                   relative_floor=meta['relative_floor'],
                   fd_error_level=meta['fd_error_level'], units=meta['units'])


# --- public API --------------------------------------------------------------
def _fd_error_warning(sens, flags, level, amp, result):
    """The rank warning (``''`` if none), issued when
    :func:`_rank_test_doubtful` holds for ``level``, the
    :func:`_fd.error_level` of the selected columns' ``flags``, scaled by
    ``amp``, the realized weight amplification of a lowered
    ``relative_floor`` (1 otherwise), and the rank test of ``result``;
    with the remedies that apply (judged on the unscaled components)."""
    s, rank = result.singular_values, result.numerical_rank
    rank_rtol = result.rank_rtol
    total = level.total
    effective = _amplified(total, amp)
    if not _rank_test_doubtful(effective, rank_rtol, s, rank):
        return ''
    ratio = s[rank - 1] / s[0]
    scaled = (f', scaled by {amp:.2g} to {effective:.2g} (relative_floor='
              f'{result.relative_floor:g}, below the default '
              f'{RELATIVE_FLOOR:g}, lets rows near 0 weigh that much more),'
              if amp > 1 else '')
    parts = [
        'The estimated finite-difference error level of these '
        'sensitivities (integration noise plus truncation error, about '
        f'{total:.2g} of the largest singular value){scaled} exceeds '
        f'rank_rtol={rank_rtol:g}, and the smallest singular value the '
        f'rank test kept is within {FD_ERROR_MARGIN:g} x of it '
        f'(singular_values[{rank - 1}] / singular_values[0] = '
        f'{ratio:.1g}): the rank test cannot tell an exact '
        'non-identifiability lifted by that error, which reads as a finite '
        'se, from a genuine weak direction here.']
    if amp > 1:
        parts.append(f'Raise relative_floor to the default '
                     f'{RELATIVE_FLOOR:g}, which bounds the weight of those '
                     'rows.')
    half = rank_rtol / 2
    truncation = max(level.one_sided, level.central)
    if level.noise > half:
        # rtol below `needed` brings noise + truncation under rank_rtol
        # (the noise alone if the truncation is above it); rounded down.
        atol, rtol = sens.tolerances
        room = rank_rtol - truncation if truncation < rank_rtol else rank_rtol
        needed = room * level.h_min / _fd.FD_NOISE
        scale = 10.0 ** math.floor(math.log10(needed))
        needed = math.floor(needed / scale) * scale
        loose = (f'; atol={atol:g} must come down too, to at most '
                 f'{_fd.FD_TOLERANCES[0]:g} as by default'
                 if atol > _fd.FD_TOLERANCES[0] else '')
        parts.append(
            f'Integration noise gives about {level.noise:.1g} '
            f'({_fd.FD_NOISE:g} x the integrator rtol {rtol:g} / the '
            f'smallest relative step {level.h_min:g} in sens.flags): to '
            'lower it, recompute the sensitivities with '
            f'tolerances=({min(atol, _fd.FD_TOLERANCES[0]):g}, rtol), rtol '
            f'below {needed:.1g}{loose}.')
    if level.one_sided > half:
        fix = ("recompute them with method='central' (a smaller step does "
               'not help: the noise grows as 1/step)'
               if sens.method == 'forward' else
               "they fell back to one side because a run on the other side "
               'failed (sens.flags); leave them out with params= to check '
               'the verdict')
        parts.append(
            f'The one-sided differences of {level.one_sided_names} '
            "(sens.flags['scheme']) are first-order accurate, with a "
            f'truncation error of about {level.one_sided:.1g}: {fix}.')
    if level.central > half:
        parts += _central_advice(sens, flags, half, rank_rtol)
    parts.append(
        'To check the verdict, recompute at rel_step / 3 with the same '
        'tolerances and compare that ratio (1 / condition_number at full '
        'rank): it falls about 9x for a null direction lifted by '
        'central-difference truncation error (about 3x for one-sided '
        'differences), rises about 3x for one lifted by integration noise, '
        'and barely changes for a genuine weak direction (see the Notes of '
        'identifiability).')
    return ' '.join(parts)


def _central_advice(sens, flags, half, rank_rtol):
    """Remedies for the reliable central differences of ``flags`` whose
    truncation error exceeds ``half``: a smaller ``abs_step`` for linear
    columns whose step it set (wider than ``rel_step``, or at ``p > 0``),
    ``rel_step`` at most 1e-3 for the others, and none for steps already
    at most 1e-3 (only below the default ``rank_rtol``)."""
    groups = {'abs': [], 'rel': [], 'small': []}
    for name, row in flags.iterrows():
        if not row['reliable'] or row['scheme'] in ('forward', 'backward'):
            continue
        h = _fd.relative_step(row['value'], row['log'], row['step'])
        if not (h > 0 and _fd.TRUNCATION_CENTRAL * h ** 2 > half):
            continue
        by_abs = not row['log'] and (row['value'] > 0
                                     or h > sens.rel_step * (1 + 1e-9))
        key = 'abs' if by_abs else 'rel' if h > 1e-3 else 'small'
        groups[key].append((name, h))
    parts = []
    for key, rows in groups.items():
        if not rows:
            continue
        names = [n for n, _ in rows]
        h = max(h for _, h in rows)
        error = _fd.TRUNCATION_CENTRAL * h ** 2
        if key == 'abs':
            parts.append(
                f'The central differences of {names} use an abs_step of up '
                f'to {h:g} of the value, with a truncation error of about '
                f'{error:.1g}: recompute the sensitivities with a smaller '
                'abs_step for them (at most about 1e-3 x |value|).')
        elif key == 'rel':
            parts.append(
                f'Central differences at a relative step of {h:g} '
                f'({names}) have a truncation error of about {error:.1g}: '
                'recompute the sensitivities with rel_step at most 1e-3.')
        else:
            parts.append(
                f'The central differences of {names}, at relative steps of '
                f'at most {h:g}, have a truncation error of about '
                f'{error:.1g}, too large for rank_rtol={rank_rtol:g}.')
    return parts


def _time_indices(grid, times):
    """``(indices, order)``: the indices of ``times`` on the output ``grid``
    in increasing order (all if ``None``), and the permutation that sorts
    ``times`` into that order (``times[order]`` is sorted)."""
    if times is None:
        return np.arange(len(grid)), np.arange(len(grid))
    if isinstance(times, (str, bytes, Mapping)):
        raise TypeError(f'times must be a sequence of numbers; got '
                        f'{type(times).__name__}.')
    try:
        t = np.asarray(times, dtype=float)
    except (TypeError, ValueError) as e:
        raise TypeError(f'times must be a sequence of numbers: {e}') from e
    if t.ndim == 0:
        t = t[None]
    if t.ndim != 1 or t.size == 0 or not np.all(np.isfinite(t)):
        raise ValueError('times must be a non-empty 1-D sequence of finite '
                         'numbers.')
    span = max(1.0, float(np.max(np.abs(grid)))) if len(grid) else 1.0
    idx = []
    for value in t:
        j = int(np.argmin(np.abs(grid - value)))
        if abs(grid[j] - value) > 1e-9 * span:
            raise ValueError(f'time {value!r} is not an output time of the '
                             'sensitivity run; put it on the condition\'s '
                             'output grid (times=...).')
        idx.append(j)
    if len(set(idx)) != len(idx):
        raise ValueError('times lists an output time more than once.')
    order = np.argsort(idx, kind='stable')
    return np.asarray(idx, dtype=int)[order], order


def _positive(arg, value, n, order):
    """A positive, finite scalar or length-``n`` 1-D array as an array of
    length ``n`` and its JSON-safe record (float or list). An array is given
    in the order of the user's ``times`` and returned (and recorded)
    permuted by ``order``, i.e. in increasing time like the rows."""
    if isinstance(value, (bool, np.bool_, str, bytes, Mapping)):
        raise TypeError(f'{arg} must be a positive number or a sequence of '
                        f'them; got {value!r}.')
    try:
        arr = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as e:
        raise TypeError(f'{arg} must be a positive number or a sequence of '
                        f'them: {e}') from e
    if arr.ndim == 0:
        record = float(arr)
        arr = np.full(n, record)
    elif arr.ndim == 1 and arr.size == n:
        arr = arr[order]
        record = [float(v) for v in arr]
    else:
        raise ValueError(f'{arg} must be a number or have one value per '
                         f'selected time ({n}); got shape {arr.shape}.')
    if not (np.all(np.isfinite(arr)) and np.all(arr > 0)):
        raise ValueError(f'{arg} must be positive and finite.')
    return arr, record


def _relative_value(arg, value):
    """A relative sigma as a positive, finite float."""
    value = _real(arg, value)
    if not (math.isfinite(value) and value > 0):
        raise ValueError(f'{arg} must be positive and finite; got '
                         f'{value!r}.')
    return value


def _is_tagged(value):
    """Whether ``value`` is a tagged sigma entry such as
    ``('relative', rel)`` (a tuple or list starting with a string)."""
    return (isinstance(value, (tuple, list)) and len(value) > 0
            and isinstance(value[0], str))


def _sigma_plan(sens, sigma, n_t, order):
    """``(plan, record)``: ``plan`` lists ``(output, absolute sigma array or
    None, relative sigma or None)``; ``record`` is the JSON-safe form."""
    outputs = sens.outputs
    if isinstance(sigma, Mapping):
        if not sigma:
            raise ValueError('sigma must name at least one output.')
        plan, record = [], {}
        for out, value in sigma.items():
            if out not in outputs:
                raise ValueError(f'sigma names {out!r}, which is not an '
                                 f'output of the sensitivity run; outputs: '
                                 f'{outputs}.')
            if _is_tagged(value):
                # A relative entry among absolute ones: {out: ('relative',
                # rel)}.
                if len(value) != 2 or value[0] != 'relative':
                    raise TypeError(
                        f'sigma[{out!r}] must be a positive number, one per '
                        "selected time, or ('relative', rel); got "
                        f'{value!r}.')
                rel = _relative_value(f'relative sigma of {out!r}', value[1])
                record[out] = ('relative', rel)
                plan.append((out, None, rel))
                continue
            arr, record[out] = _positive(f'sigma[{out!r}]', value, n_t,
                                         order)
            plan.append((out, arr, None))
        return plan, record
    if (isinstance(sigma, (tuple, list)) and len(sigma) == 2
            and isinstance(sigma[0], str) and sigma[0] == 'relative'):
        rel = sigma[1]
        if isinstance(rel, Mapping):
            raise TypeError(
                "rel in sigma=('relative', rel) must be one number for every "
                'output; per-output relative sigmas are written as a '
                "mapping, e.g. {'[S]': ('relative', 0.1), '[P]': "
                f"('relative', 0.05)}}; got {sigma!r}.")
        # One value for every output: name the argument, not '[S]'.
        rel = _relative_value("rel in sigma=('relative', rel)", rel)
        return [(o, None, rel) for o in outputs], ('relative', rel)
    raise TypeError("sigma must be a mapping {output: sigma or ('relative', "
                    "rel)} or ('relative', rel); got " f'{sigma!r}.')


def identifiability(sens, sigma, *, times=None, params=None,
                    rank_rtol=RANK_RTOL, relative_floor=RELATIVE_FLOOR):
    """Practical identifiability of parameters from local sensitivities.

    Builds the weighted, log-scaled sensitivity matrix
    ``S_ij = (p_j dy_i / dp_j) / sigma_i`` -- one row per measured output
    and time -- and from it the Fisher information ``FIM = S^T S`` (that is
    ``S^T W S`` with ``W = diag(1 / sigma^2)`` for the unweighted
    log-scaled sensitivities), its singular values, condition number,
    numerical rank and null-space directions, the Cramer-Rao standard
    errors and Brun et al. (2001) importance ``delta_msqr``; see
    :class:`IdentifiabilityResult` for the collinearity indices and the
    subset selection.

    Parameters
    ----------
    sens : LocalSensitivityResult
        From :func:`local_sensitivity`; its output time courses are used
        (not its metrics).
    sigma : mapping or tuple
        Measurement noise (standard deviation) of the outputs that are
        measured: ``{output: sigma}`` with absolute values in the output's
        units (a number, or one value per selected time, in the order of
        ``times``; rows, and the recorded sigma, are sorted by time), or
        ``('relative', rel)`` for a relative error ``rel`` on every output,
        bounded below: the row of output ``y`` at time ``t_i`` has
        ``sigma_i = rel * max(|y_i|, relative_floor * max_t |y|)``, the
        maximum taken over the output's whole time grid in ``sens`` (see
        Notes). A mapping can mix both, and gives outputs their own
        ``rel``: ``{'[S]': 0.5, '[P]': ('relative', 0.05)}``. Outputs not
        named are not measured.
    times : sequence of float, optional
        Measurement times (model time units); each must be an output time
        of ``sens``. Defaults to all output times (including ``t0``, where
        only initial values have an effect).
    params : sequence of str, optional
        Subset of ``sens.parameters`` (default: all).
    rank_rtol : float
        Singular values ``<= rank_rtol * s_max`` count as zero (default
        ``1e-6``; see Notes).
    relative_floor : float
        Detection limit of the relative noise model, as a fraction of each
        output's peak: ``|y|`` counts as at least
        ``relative_floor * max_t |y|``, so
        ``sigma >= rel * relative_floor * max_t |y|``; in ``(0, 1]``
        (default ``1e-3``). Not used with an absolute ``sigma``.

    Returns
    -------
    IdentifiabilityResult
        Its ``relative_floor`` is recorded when some output with a relative
        sigma contributes rows (``None`` when all rows are absolute), and
        its ``fd_error_level`` and ``rank_test_doubtful`` record whether
        the finite-difference error estimate casts doubt on the rank test;
        no doubt is no proof (see Notes).

    Raises
    ------
    TypeError, ValueError
        For an invalid ``sigma``, ``times`` (not on the output grid),
        ``params``, ``rank_rtol`` or ``relative_floor`` (all checked first);
        if ``sens`` has no output time courses (``outputs=[]``); if under a
        relative noise model every measured output is 0 at every output
        time, or an output's nominal values are not finite or its sigma is
        not finite and positive (e.g. an overflowing ``rel``); if a sigma
        is so small that the weighted sensitivities overflow; or if a
        selected parameter has a non-finite sensitivity (a failed
        perturbation, see ``sens.flags['error']``): exclude it with
        ``params``.

    Warns
    -----
    UserWarning
        One warning per condition that applies:

        - a selected parameter's finite differences were unreliable (an
          event count changed, see ``sens.flags``);
        - a selected parameter is at ``p <= 0`` (a ``'linear'`` column,
          see Notes);
        - under a relative noise model, an output is 0 at every output
          time (it is left out);
        - for two or more parameters, the finite-difference error level is
          unknown (an ``atol`` above 1e-8, or a step so small that the
          estimate overflows; ``fd_error_level`` is then ``None``);
        - for two or more parameters, the rank test is doubtful
          (``rank_test_doubtful=True``, see Notes); the warning names the
          remedies and ends with the step-convergence check.

    Notes
    -----
    Columns are in ``ln p`` for ``p > 0``. A parameter at ``p <= 0`` has a
    ``'linear'`` column, ``dy/dp / sigma`` per model unit of ``p``, which
    makes the other parameters' verdicts depend on that unit: judge them
    with ``params=`` excluding it. The analysis is local, at the nominal
    values and for the condition of ``sens``.

    The relative noise model is bounded below,
    ``sigma = rel * max(|y|, relative_floor * max_t |y|)``, like a
    detection limit at a small fraction of the output's peak over the
    whole grid of ``sens``: a pure ``rel * |y|`` gives unbounded weight to
    rows near 0, whose finite-difference errors do not shrink with ``y``.
    Keep ``relative_floor`` at the default: smaller floors let that noise
    back in, and the rank test scales its error estimate by the realized
    weight amplification of the rows (at most ``1e-3 / relative_floor``).

    An exact non-identifiability shows up as a singular value at the
    finite-difference error level, not at 0: about ``3 rtol / h`` of
    integration noise plus ``0.3 h^2`` (central) or ``h / 2`` (one-sided)
    of truncation error where the compensation is not linear in ``ln p``
    (``(k1 + k2) * S``; ``k1 * k2`` cancels). Keep central differences,
    ``rel_step=1e-3`` and the default tolerances (an ``atol`` above 1e-8
    is not modelled). The rank warning and ``rank_test_doubtful`` flag an
    estimate above ``rank_rtol`` with the smallest retained singular value
    within 100 times it. They are a screen, not a proof: at a depletion
    front (Monod growth, Michaelis-Menten kinetics with a small ``Km``)
    the truncation error is far larger, and a null can read full rank and
    ``identifiable`` with no warning.

    The decisive test is step convergence: recompute
    :func:`local_sensitivity` at ``rel_step / 3`` with the same tolerances
    and compare ``singular_values[numerical_rank - 1] /
    singular_values[0]`` (``1 / condition_number`` at full rank). A null
    direction lifted by truncation error falls about 9 times for central
    differences (as ``h^2``) and 3 times for one-sided ones (as ``h``),
    one lifted by integration noise grows about 3 times (as ``1 / h``),
    and a genuine weak direction barely changes.
    :ref:`identifiability-noise` explains these rules, with an example.

    A parameter is on a null direction (``se = inf``,
    ``table['null_direction']``) when its loading on the null space exceeds
    ``1e-3``; ``table['identifiable']`` also requires a relative standard
    error of at most 1. The Cramer-Rao bound is a lower bound on the
    standard error of any unbiased estimator, valid for small, Gaussian
    errors.

    Examples
    --------
    >>> ident = identifiability(ls, sigma={'[P]': 0.5})   # doctest: +SKIP
    >>> ident.numerical_rank, ident.null_directions        # doctest: +SKIP
    >>> ident.collinearity(subset_size=2).head()           # doctest: +SKIP
    >>> ident.best_subset(gamma_max=10)                    # doctest: +SKIP
    """
    if not isinstance(sens, LocalSensitivityResult):
        raise TypeError(f'sens must be a LocalSensitivityResult; got '
                        f'{type(sens).__name__}.')
    if not sens.outputs:
        raise ValueError(
            'identifiability needs output time courses, and this '
            'LocalSensitivityResult has none (it was computed with '
            'outputs=[]); metric sensitivities are not used. Rerun '
            'local_sensitivity with outputs=[...].')
    names = (list(sens.parameters) if params is None
             else _id_list('params', params))
    unknown = [p for p in names if p not in sens.parameters]
    if unknown:
        raise ValueError(f'params {unknown} are not parameters of the '
                         f'sensitivity run; parameters: {sens.parameters}.')
    rank_rtol = _check_rank_rtol(rank_rtol)
    relative_floor = _check_relative_floor(relative_floor)
    cols = [sens.parameters.index(p) for p in names]
    t_idx, order = _time_indices(sens.time, times)
    plan, record = _sigma_plan(sens, sigma, len(t_idx), order)
    values = sens.values[cols]
    space = ['log' if v > 0 else 'linear' for v in values]
    scale = np.where(values > 0, values, 1.0)
    # Whether some relative block is built (an always-zero output is left
    # out): only then does the floor shape the matrix and get recorded.
    relative_rows = False
    # The realized weight amplification of a relative_floor below the
    # default: max over the relative rows of sigma_default / sigma (>= 1),
    # by which their finite-difference error is amplified too.
    amp = 1.0
    blocks, row_outputs, zero = [], [], []
    for out, absolute, relative in plan:
        k = sens.outputs.index(out)
        if absolute is None:
            # Bounded relative noise: never below relative_floor of the
            # output's peak over the whole grid (a detection limit).
            with np.errstate(all='ignore'):
                peak = float(np.max(np.abs(sens.nominal[:, k])))
            if not math.isfinite(peak):
                raise ValueError(f'The nominal values of {out!r} are not '
                                 'finite, so its relative sigma is '
                                 'undefined.')
            if peak == 0:
                zero.append(out)
                continue
            sig = _relative_sigma(sens.nominal[t_idx, k], relative, peak,
                                  relative_floor)
            if not np.all(np.isfinite(sig) & (sig > 0)):
                raise ValueError(
                    f'The relative sigma of {out!r} is not finite and '
                    f'positive at the selected times (rel={relative!r}, '
                    f'relative_floor={relative_floor!r}).')
            relative_rows = True
            with np.errstate(all='ignore'):       # sig > 0: no 0 / 0
                y = np.abs(sens.nominal[t_idx, k])
                amp = max(amp, float(np.max(
                    np.maximum(y, RELATIVE_FLOOR * peak)
                    / np.maximum(y, relative_floor * peak))))
        else:
            sig = absolute
        with np.errstate(all='ignore'):
            deriv = sens.raw[t_idx][:, k, :][:, cols] * scale
            block = deriv / sig[:, None]
        if not np.all(np.isfinite(block) | ~np.isfinite(deriv)):
            raise ValueError(f'The sigma of {out!r} is too small: the '
                             'weighted sensitivities overflow.')
        blocks.append(block)
        row_outputs += [out] * len(t_idx)
    if zero:
        why = (f'{zero} are 0 at every output time of the sensitivity run, '
               'so they have no relative sigma')
        advice = ('To measure them, give them an absolute sigma in a '
                  "mapping, e.g. {output: sigma, other: ('relative', rel)}.")
        if not blocks:
            raise ValueError(f'No measured output: under the relative noise '
                             f'model, {why}. {advice}')
        warnings.warn(f'Under the relative noise model, {why}; they are left '
                      f'out. {advice}', UserWarning, stacklevel=2)
    matrix = np.vstack(blocks)
    row_times = np.tile(sens.time[t_idx], len(blocks))
    bad = [p for j, p in enumerate(names)
           if not np.all(np.isfinite(matrix[:, j]))]
    if bad:
        raise ValueError(f'Non-finite sensitivities for {bad} (see '
                         "sens.flags['error']); exclude them with params=.")
    flags = sens.flags.iloc[cols]
    reliable = flags['reliable'].to_numpy(dtype=bool)
    unreliable = [p for p, ok in zip(names, reliable) if not ok]
    if unreliable:
        warnings.warn(f'The finite-difference sensitivities of {unreliable} '
                      'are unreliable (an event count changed; see '
                      'sens.flags); their identifiability measures are not '
                      'meaningful.', UserWarning, stacklevel=2)
    linear = [p for p, sp in zip(names, space) if sp == 'linear']
    if linear:
        weight = (' Under a relative sigma, the t0 row of an output that '
                  'starts at 0 is weighted 1 / (rel * relative_floor * '
                  'max_t |y|).' if relative_rows else '')
        warnings.warn(
            f'{linear} are at p <= 0, so their columns are dy/dp / sigma '
            'per model unit of p, not per ln p: the singular values, '
            'condition number and numerical rank depend on that unit, and '
            "with them the other parameters' se, null directions and "
            'identifiable verdicts, delta_msqr and the order of '
            f'best_subset.{weight} Judge the other parameters with params= '
            'excluding them.', UserWarning, stacklevel=2)
    if (len(names) > 1 and sens.tolerances is not None
            and sens.tolerances[0] > _fd.ATOL_MAX):
        atol = sens.tolerances[0]
        warnings.warn(
            f'The sensitivities were computed with atol={atol:g}, but the '
            'finite-difference error level behind the rank test assumes '
            f'atol <= {_fd.ATOL_MAX:g} (the default is '
            f'{_fd.FD_TOLERANCES[0]:g}, scaled to the outputs): rows of '
            'depleted outputs or of outputs that start at 0 can carry '
            'noise that lifts an exact non-identifiability to full rank '
            'with a finite se. Recompute them with '
            f'tolerances=({_fd.FD_TOLERANCES[0]:g}, rtol).', UserWarning,
            stacklevel=2)
    level = _fd.error_level(flags, sens.tolerances)
    # A relative step so small that the noise estimate overflows leaves the
    # level unknown (the step itself was warned about by local_sensitivity).
    finite = math.isfinite(level.total)
    if len(names) > 1 and not finite:
        warnings.warn(
            'The finite-difference error level of these sensitivities could '
            f'not be estimated (a relative step of {level.h_min:g} in '
            'sens.flags): fd_error_level and rank_test_doubtful are None.',
            UserWarning, stacklevel=2)
    # The level models rtol only: above ATOL_MAX it is unknown, so the
    # result records no level (and no verdict on the rank test).
    modelled = finite and (sens.tolerances is None
                           or not sens.tolerances[0] > _fd.ATOL_MAX)
    result = IdentifiabilityResult(
        parameters=names, values=values, space=space, sensitivity=matrix,
        row_outputs=row_outputs, row_times=row_times, reliable=reliable,
        sigma=record, rank_rtol=rank_rtol,
        relative_floor=relative_floor if relative_rows else None,
        fd_error_level=_amplified(level.total, amp) if modelled else None,
        units=dict(sens.units))
    # Also with an unmodelled atol, whose own warning comes first.
    message = (_fd_error_warning(sens, flags, level, amp, result)
               if finite else '')
    if message:
        warnings.warn(message, UserWarning, stacklevel=2)
    return result
