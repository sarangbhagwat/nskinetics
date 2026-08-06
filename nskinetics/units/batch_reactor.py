# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import numpy as np

from ..utils import get_index_nearest_element_from_sorted_array

__all__ = ('NSKBatchReactor', 'AerationSpec', 'SpikeReduceRetry',
           'select_tau_index')


def select_tau_index(results, col_names, tau, policy, n_decimals=2):
    """Select the results-array row index for the reported ``tau``.

    Parameters
    ----------
    results : numpy.ndarray
        2-D array of simulation results (rows = timesteps).
    col_names : list of str
        Column names for ``results``.
    tau : float or None
        Reaction time; used only when ``policy is None``.
    policy : None or tuple
        ``None`` -> nearest-time row. ``('max'|'min', var)`` -> first row where
        ``var`` (rounded to ``n_decimals``) equals its rounded max/min.
        ``('equals', var, value)`` -> first row where ``var`` (rounded) equals
        ``value`` (rounded).
    n_decimals : int
        Rounding used by the ``max``/``min``/``equals`` policies.

    Returns
    -------
    (index, success) : tuple[int, bool]
        ``success`` is ``False`` (and ``index`` is ``-1``) only when an
        ``'equals'`` policy finds no match.
    """
    if policy is None:
        idx = get_index_nearest_element_from_sorted_array(
            results[:, col_names.index('time')], tau)
        return idx, True

    if policy[0] in ('max', 'min'):
        col = col_names.index(policy[1])
        column = results[:, col]
        target = getattr(column, policy[0])()
        idx = np.where(np.round(column, n_decimals) == np.round(target, n_decimals))[0][0]
        return idx, True

    if policy[0] == 'equals':
        col = col_names.index(policy[1])
        column = results[:, col]
        matches = np.where(np.round(column, n_decimals) == np.round(policy[2], n_decimals))[0]
        if len(matches):
            return int(matches[0]), True
        return -1, False

    raise ValueError(f'Unknown tau_update_policy {policy!r}.')


class AerationSpec:  # replaced in Task 8
    def __init__(self, *a, **k):
        raise NotImplementedError


class SpikeReduceRetry:  # replaced in Task 8
    def __init__(self, *a, **k):
        raise NotImplementedError


class NSKBatchReactor:  # replaced in Task 9
    def __init__(self, *a, **k):
        raise NotImplementedError
