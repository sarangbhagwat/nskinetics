# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Parameter ranges and seeded samplers (shared, mostly private).

:class:`ParamRange` is public (re-exported from :mod:`nskinetics.analysis`).
The helpers are the one sampling layer behind screening (Morris/Sobol
designs), calibration multistarts and Laplace posterior draws:

- :func:`_seed_sequence` / :func:`_rng`: every random stream comes from an
  explicit seed through :class:`numpy.random.SeedSequence`;
  :func:`_child_seed` derives independent streams from one seed without
  mutating it and :func:`_seed_record` gives the int that reproduces it;
- :func:`_unit_sample`: scrambled Sobol' or Latin-hypercube points in the
  unit cube (``scipy.stats.qmc``);
- :func:`_as_ranges` / :func:`_scale` / :func:`_sample`: map unit-cube
  points onto parameter ranges (uniform in linear or log space);
- :func:`_gaussian`: multivariate normal draws (for Laplace sampling).
"""
import math
import numbers
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import qmc

from ._common import _finite

__all__ = ('ParamRange',)

_METHODS = ('sobol', 'lhs')


@dataclass(frozen=True)
class ParamRange:
    """A parameter range, sampled uniformly in linear or log space.

    Parameters
    ----------
    lo, hi : float
        Lower and upper bounds, in the parameter's own (model) units;
        ``lo < hi``.
    log : bool
        Sample uniformly in ``ln(value)`` (requires ``lo > 0``) instead of
        in ``value``. Defaults to ``False``.

    Raises
    ------
    TypeError
        For non-numeric bounds or a non-bool ``log``.
    ValueError
        For non-finite bounds, ``lo >= hi``, or ``log=True`` with
        ``lo <= 0``.

    Examples
    --------
    >>> r = ParamRange(0.05, 5, log=True)
    >>> round(r.from_unit(0.5), 12)      # geometric midpoint
    0.5
    """
    lo: float
    hi: float
    log: bool = False

    def __post_init__(self):
        lo, hi = _finite('lo', self.lo), _finite('hi', self.hi)
        if not isinstance(self.log, (bool, np.bool_)):
            raise TypeError(f'log must be a bool; got {self.log!r}.')
        if not lo < hi:
            raise ValueError(f'ParamRange needs lo < hi; got lo={lo!r}, '
                             f'hi={hi!r}.')
        if self.log and lo <= 0:
            raise ValueError(f'A log-scaled ParamRange needs lo > 0; got '
                             f'lo={lo!r}.')
        object.__setattr__(self, 'lo', lo)
        object.__setattr__(self, 'hi', hi)
        object.__setattr__(self, 'log', bool(self.log))

    @property
    def bounds(self):
        """``(lo, hi)`` in the parameter's units (tuple of float)."""
        return (self.lo, self.hi)

    @property
    def transformed_bounds(self):
        """``(ln lo, ln hi)`` for a log range, else ``(lo, hi)``."""
        if self.log:
            return (math.log(self.lo), math.log(self.hi))
        return (self.lo, self.hi)

    def to_transformed(self, value):
        """Map values to the sampling space.

        Parameters
        ----------
        value : float or array_like
            Values in the parameter's own (model) units.

        Returns
        -------
        float or numpy.ndarray
            ``ln(value)`` for a log range (``NaN`` for ``value < 0``,
            ``-inf`` for 0), else ``value``; a float for a scalar input.
        """
        x = np.asarray(value, dtype=float)
        if self.log:
            with np.errstate(all='ignore'):
                x = np.log(x)
        return float(x) if x.ndim == 0 else x

    def from_transformed(self, z):
        """Inverse of :meth:`to_transformed`.

        Parameters
        ----------
        z : float or array_like
            Values in the sampling space (``ln`` of the model units for a
            log range).

        Returns
        -------
        float or numpy.ndarray
            ``exp(z)`` for a log range, else ``z`` (model units); a float
            for a scalar input. Not clipped to the range.
        """
        z = np.asarray(z, dtype=float)
        if self.log:
            with np.errstate(all='ignore'):
                z = np.exp(z)
        return float(z) if z.ndim == 0 else z

    def from_unit(self, u):
        """Map unit-interval values ``u`` in ``[0, 1]`` onto the range.

        Parameters
        ----------
        u : float or array_like
            Positions in ``[0, 1]`` (0 maps to ``lo``, 1 to ``hi``).

        Returns
        -------
        float or numpy.ndarray
            ``lo + u (hi - lo)``, or ``exp(ln lo + u (ln hi - ln lo))`` for a
            log range, clipped to ``[lo, hi]`` against rounding.
        """
        a, b = self.transformed_bounds
        u = np.asarray(u, dtype=float)
        with np.errstate(all='ignore'):
            x = np.clip(self.from_transformed(a + u * (b - a)),
                        self.lo, self.hi)
            # exact end points (exp(ln lo) can be lo + 1 ulp)
            x = np.where(u <= 0, self.lo, np.where(u >= 1, self.hi, x))
        return float(x) if np.ndim(x) == 0 else x

    def to_unit(self, value):
        """Inverse of :meth:`from_unit`.

        Parameters
        ----------
        value : float or array_like
            Values in the parameter's own (model) units.

        Returns
        -------
        float or numpy.ndarray
            Positions ``(x - a) / (b - a)`` with ``x``, ``a``, ``b`` the
            value and bounds in the sampling space (``ln`` for a log
            range): 0 at ``lo``, 1 at ``hi``; values outside the range map
            outside ``[0, 1]``.
        """
        a, b = self.transformed_bounds
        z = np.asarray(self.to_transformed(value), dtype=float)
        with np.errstate(all='ignore'):
            u = (z - a) / (b - a)
        return float(u) if u.ndim == 0 else u

    def contains(self, value):
        """Whether values lie in the closed range.

        Parameters
        ----------
        value : float or array_like
            Values in the parameter's own (model) units.

        Returns
        -------
        bool or numpy.ndarray
            ``lo <= value <= hi`` (``False`` for ``NaN``); a bool for a
            scalar input.
        """
        v = np.asarray(value, dtype=float)
        out = (v >= self.lo) & (v <= self.hi)
        return bool(out) if out.ndim == 0 else out


def _as_ranges(bounds):
    """``{name: ParamRange}`` from ``{name: ParamRange | (lo, hi)}``.

    Raises
    ------
    TypeError, ValueError
        For a non-mapping, an empty mapping, a non-string name or an invalid
        range.
    """
    if not isinstance(bounds, Mapping):
        raise TypeError('bounds must be a mapping {name: (lo, hi) or '
                        f'ParamRange}}; got {type(bounds).__name__}.')
    if not bounds:
        raise ValueError('bounds must name at least one parameter.')
    out = {}
    for name, spec in bounds.items():
        if not isinstance(name, str) or not name:
            raise TypeError(f'bounds keys must be non-empty strings; got '
                            f'{name!r}.')
        if isinstance(spec, ParamRange):
            out[name] = spec
        elif (isinstance(spec, (tuple, list)) and len(spec) == 2):
            try:
                out[name] = ParamRange(*spec)
            except (TypeError, ValueError) as e:
                raise type(e)(f'bounds[{name!r}]: {e}') from e
        else:
            raise TypeError(f'bounds[{name!r}] must be a (lo, hi) pair or a '
                            f'ParamRange; got {spec!r}.')
    return out


def _seed_sequence(seed):
    """A :class:`numpy.random.SeedSequence` from ``seed``.

    Parameters
    ----------
    seed : None, int or numpy.random.SeedSequence
        ``None`` draws fresh OS entropy (not reproducible); a non-negative
        int or a SeedSequence is reproducible.

    Raises
    ------
    TypeError, ValueError
        For any other type, or a negative int.
    """
    if seed is None:
        return np.random.SeedSequence()
    if isinstance(seed, np.random.SeedSequence):
        return seed
    if isinstance(seed, bool) or not isinstance(seed, numbers.Integral):
        raise TypeError(f'seed must be None, a non-negative int or a '
                        f'numpy SeedSequence; got {seed!r}.')
    if seed < 0:
        raise ValueError(f'seed must be non-negative; got {seed!r}.')
    return np.random.SeedSequence(int(seed))


def _rng(seed):
    """A :class:`numpy.random.Generator` seeded through :func:`_seed_sequence`."""
    return np.random.default_rng(_seed_sequence(seed))


def _child_seed(ss, i):
    """The ``i``-th independent child stream of the SeedSequence ``ss``.

    Same as the ``i``-th child of :meth:`numpy.random.SeedSequence.spawn`
    on a fresh sequence, but without mutating ``ss`` (``spawn`` advances its
    counter), so one user seed gives the same streams on every call, e.g.
    a design (child 0) and its bootstrap (child 1), or calibration starts
    and posterior draws.
    """
    return np.random.SeedSequence(ss.entropy, spawn_key=(*ss.spawn_key, i),
                                  pool_size=ss.pool_size)


def _seed_record(ss):
    """The int that reproduces ``ss`` when passed back as ``seed=`` (its
    entropy; for ``seed=None`` the fresh OS entropy drawn), or ``None`` for
    a child sequence or one built from several ints."""
    if ss.spawn_key or not isinstance(ss.entropy, numbers.Integral):
        return None
    return int(ss.entropy)


def _unit_sample(n, d, *, method='sobol', seed=None, scramble=True):
    """``n`` points in the ``d``-dimensional unit cube.

    Parameters
    ----------
    n, d : int
        Number of points (>= 1) and dimensions (>= 1).
    method : {'sobol', 'lhs'}
        Scrambled Sobol' sequence (``scipy.stats.qmc.Sobol``; the first
        ``n`` points of a ``2**ceil(log2 n)`` block, so a power of two gives
        the balanced set) or a Latin hypercube
        (``scipy.stats.qmc.LatinHypercube``).
    seed : None, int or SeedSequence
        See :func:`_seed_sequence`; the same seed gives the same points.
    scramble : bool
        Scramble the Sobol' sequence (default ``True``); ignored for 'lhs'.

    Returns
    -------
    numpy.ndarray
        Shape ``(n, d)``, values in ``[0, 1)``.
    """
    for arg, value in (('n', n), ('d', d)):
        if isinstance(value, bool) or not isinstance(value, numbers.Integral):
            raise TypeError(f'{arg} must be an integer; got {value!r}.')
        if value < 1:
            raise ValueError(f'{arg} must be >= 1; got {value!r}.')
    if method not in _METHODS:
        raise ValueError(f'method must be one of {_METHODS}; got {method!r}.')
    rng = _rng(seed)
    if method == 'sobol':
        engine = qmc.Sobol(int(d), scramble=bool(scramble), rng=rng)
        m = max(0, math.ceil(math.log2(int(n))))
        return engine.random_base2(m)[:int(n)]
    return qmc.LatinHypercube(int(d), rng=rng).random(int(n))


def _scale(unit, ranges):
    """Map unit-cube points (``(n, d)``) column-wise onto ``ranges`` (a
    sequence of :class:`ParamRange`); returns a new ``(n, d)`` array."""
    unit = np.asarray(unit, dtype=float)
    ranges = list(ranges)
    if unit.ndim != 2 or unit.shape[1] != len(ranges):
        raise ValueError(f'unit must have shape (n, {len(ranges)}); got '
                         f'{unit.shape}.')
    out = np.empty_like(unit)
    for j, r in enumerate(ranges):
        out[:, j] = r.from_unit(unit[:, j])
    return out


def _sample(bounds, n, *, method='sobol', seed=None):
    """``n`` samples of the parameters in ``bounds`` as a DataFrame.

    Parameters
    ----------
    bounds : mapping
        ``{name: (lo, hi) or ParamRange}`` (see :func:`_as_ranges`).
    n : int
        Number of samples.
    method, seed
        See :func:`_unit_sample`.

    Returns
    -------
    pandas.DataFrame
        One column per parameter (in the order of ``bounds``), ``n`` rows,
        values in the parameters' own units.
    """
    ranges = _as_ranges(bounds)
    unit = _unit_sample(n, len(ranges), method=method, seed=seed)
    return pd.DataFrame(_scale(unit, ranges.values()), columns=list(ranges))


def _gaussian(mean, cov, n, *, seed=None, lo=None, hi=None):
    """``n`` draws from a multivariate normal ``N(mean, cov)``.

    The covariance may be singular (positive semi-definite): draws use its
    symmetric eigendecomposition with negative round-off eigenvalues set to
    zero, so they stay in the covariance's range. Optional ``lo`` / ``hi``
    clip each column (e.g. to calibration bounds in transformed space).

    Parameters
    ----------
    mean : array_like, shape (d,)
    cov : array_like, shape (d, d)
        Symmetric; finite.
    n : int
        Number of draws (>= 1).
    seed : None, int or SeedSequence
        See :func:`_seed_sequence`.
    lo, hi : array_like, optional
        Per-column clip bounds.

    Returns
    -------
    numpy.ndarray
        Shape ``(n, d)``.
    """
    mean = np.asarray(mean, dtype=float)
    cov = np.asarray(cov, dtype=float)
    if mean.ndim != 1 or cov.shape != (mean.size, mean.size):
        raise ValueError(f'cov must have shape ({mean.size}, {mean.size}) to '
                         f'match mean; got {cov.shape}.')
    if not (np.all(np.isfinite(mean)) and np.all(np.isfinite(cov))):
        raise ValueError('mean and cov must be finite.')
    if isinstance(n, bool) or not isinstance(n, numbers.Integral) or n < 1:
        raise ValueError(f'n must be a positive integer; got {n!r}.')
    with np.errstate(all='ignore'):
        w, v = np.linalg.eigh(0.5 * (cov + cov.T))
        factor = v * np.sqrt(np.clip(w, 0.0, None))
        z = _rng(seed).standard_normal((int(n), mean.size))
        draws = mean + z @ factor.T
    if lo is not None or hi is not None:
        draws = np.clip(draws, -np.inf if lo is None else lo,
                        np.inf if hi is None else hi)
    return draws
