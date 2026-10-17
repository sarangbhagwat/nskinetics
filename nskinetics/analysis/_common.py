# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Shared private helpers of the analysis package.

Argument validators and the long-format CSV reader used by every analysis
result (local sensitivity, identifiability, screening and calibration), so
each check behaves -- and fails -- the same way everywhere. Imports only the
standard library, numpy and pandas, so any analysis module can use it.
"""
import json
import math
import numbers
from collections.abc import Mapping

import numpy as np
import pandas as pd

__all__ = ()

#: Columns of the long-format CSV every analysis result writes.
_CSV_COLUMNS = ('quantity', 'key', 'parameter', 'index', 'value')
#: An empty group of that CSV (for quantities a file does not contain).
_EMPTY = pd.DataFrame(columns=list(_CSV_COLUMNS), dtype=str)


def _real(arg, value):
    """``value`` as a float, which may be infinite or ``NaN`` (for caps such
    as ``gamma_max=inf``); ``TypeError`` for a non-real or a bool."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(
            value, numbers.Real):
        raise TypeError(f'{arg} must be a real number; got {value!r}.')
    return float(value)


def _finite(arg, value):
    """``value`` as a finite float.

    Raises ``TypeError`` for anything that is not a real number (bools
    included) and ``ValueError`` for ``NaN`` or an infinity.
    """
    if isinstance(value, (bool, np.bool_)) or not isinstance(
            value, numbers.Real):
        raise TypeError(f'{arg} must be a real number; got {value!r}.')
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f'{arg} must be finite; got {value!r}.')
    return value


def _count(arg, value, minimum):
    """``value`` as an int ``>= minimum``.

    Raises ``TypeError`` for a non-integer (bools included) and
    ``ValueError`` below ``minimum``.
    """
    if isinstance(value, (bool, np.bool_)) or not isinstance(
            value, numbers.Integral):
        raise TypeError(f'{arg} must be an integer; got {value!r}.')
    if value < minimum:
        raise ValueError(f'{arg} must be >= {minimum}; got {value!r}.')
    return int(value)


def _id_list(arg, ids, *, allow_empty=False):
    """``ids`` as a list of unique non-empty strings."""
    if isinstance(ids, str):
        ids = [ids]
    if isinstance(ids, (bytes, Mapping)) or not hasattr(ids, '__iter__'):
        raise TypeError(f'{arg} must be a sequence of ids; got '
                        f'{type(ids).__name__}.')
    ids = list(ids)
    for i in ids:
        if not isinstance(i, str) or not i:
            raise TypeError(f'{arg} must hold non-empty strings; got {i!r}.')
    dup = sorted({i for i in ids if ids.count(i) > 1})
    if dup:
        raise ValueError(f'{arg} lists {dup} more than once.')
    if not ids and not allow_empty:
        raise ValueError(f'{arg} must not be empty.')
    return ids


def _check_ids(model, arg, ids, *, settable=True):
    """``model._validate_ids(ids)``, naming the argument the ids came from
    in the ``ValueError`` (e.g. ``"params: unknown model id(s) ['q']"``)."""
    try:
        model._validate_ids(list(ids), settable=settable)
    except ValueError as e:
        raise ValueError(f'{arg}: {e}') from e


def _incomplete(path, kind, detail):
    """The ``ValueError`` for a truncated or edited long CSV."""
    return ValueError(f'{path!s}: incomplete {kind} CSV ({detail}).')


def _read_long_csv(path, kind, required=()):
    """``(meta, {quantity: rows})`` of a long CSV written by ``kind.to_csv``.

    Raises ``ValueError`` if the file does not have the long-format header,
    was written for another kind of result, or lacks a ``meta`` key listed
    in ``required``.
    """
    with open(path, encoding='utf-8', newline='') as f:
        header = f.readline().rstrip('\r\n')
    if header != ','.join(_CSV_COLUMNS):
        raise ValueError(f'{path!s} was not written by {kind}.to_csv().')
    df = pd.read_csv(path, dtype=str, keep_default_na=False, encoding='utf-8')
    groups = {q: g for q, g in df.groupby('quantity', sort=False)}
    meta_rows = groups.get('meta', _EMPTY)
    meta = {key: json.loads(value) for key, value in
            zip(meta_rows['key'], meta_rows['value'])}
    if meta.get('nskinetics') != kind:
        raise ValueError(f'{path!s} was not written by {kind}.to_csv().')
    missing = [key for key in required if key not in meta]
    if missing:
        raise _incomplete(path, kind, f'missing meta keys {missing}')
    return meta, groups


def _scatter(path, kind, groups, quantity, shape, axes, parse=float,
             dtype=float):
    """The array of ``shape`` held by the rows of ``quantity`` in a long CSV
    (``groups`` from :func:`_read_long_csv`), one row per entry.

    ``axes`` gives one ``(column, mapping)`` per dimension: the column that
    holds the entry's index along it, read as ``int(value)`` if ``mapping``
    is ``None`` and as ``mapping[value]`` (a name) otherwise. Each ``value``
    is JSON, passed through ``parse``. Raises ``ValueError`` (via
    :func:`_incomplete`) unless the rows cover every entry exactly once, so
    a truncated or edited file never leaves an entry uninitialized; an
    empty ``shape`` accepts a missing group.
    """
    g = groups.get(quantity, _EMPTY)
    n = math.prod(shape)
    if len(g) != n:
        raise _incomplete(path, kind, f'{quantity}: expected {n} rows, found '
                                      f'{len(g)}')
    out = np.empty(shape, dtype=dtype)
    if n == 0:
        return out
    try:
        idx = tuple(g[column].astype(int).to_numpy() if mapping is None
                    else np.array([mapping[v] for v in g[column]], dtype=int)
                    for column, mapping in axes)
    except (KeyError, ValueError) as e:
        raise _incomplete(path, kind, f'{quantity}: unknown index {e}'
                          ) from e
    if any(np.any((i < 0) | (i >= size)) for i, size in zip(idx, shape)):
        raise _incomplete(path, kind, f'{quantity}: index out of range')
    flat = np.ravel_multi_index(idx, shape)
    if np.unique(flat).size != n:
        raise _incomplete(path, kind, f'{quantity}: duplicate entries')
    values = [parse(json.loads(v)) for v in g['value']]
    target = out.reshape(-1)                    # a view: out is contiguous
    if out.dtype == object:
        for i, value in zip(flat, values):     # keeps sequences whole
            target[i] = value
    else:
        target[flat] = values
    return out
