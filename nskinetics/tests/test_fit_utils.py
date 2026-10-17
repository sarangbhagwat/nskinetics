# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests of the repaired hybrid global optimizer (nskinetics.utils._fit)."""
import inspect

import numpy as np

from nskinetics.utils import _fit


def _rosenbrock(x):
    return float((1 - x[0]) ** 2 + 100 * (x[1] - x[0] ** 2) ** 2)


def test_default_parallel_is_serial():
    sig = inspect.signature(_fit.hybrid_global_optimize)
    assert sig.parameters['parallel'].default == 'none'


def test_global_opt_result_is_exported():
    from nskinetics.utils import GlobalOptResult, hybrid_global_optimize
    assert GlobalOptResult is _fit.GlobalOptResult
    assert 'GlobalOptResult' in _fit.__all__
    assert hybrid_global_optimize is _fit.hybrid_global_optimize


def test_nan_objective_is_not_a_success():
    calls = []

    def nan_objective(x):
        calls.append(1)
        return float('nan')

    res = _fit.hybrid_global_optimize(
        nan_objective, [(-1.0, 1.0), (-1.0, 1.0)], de_popsize=5,
        de_maxiter=3, do_basinhop=False, polish_top_k=1, local_maxiter=5,
        random_state=0)
    assert calls
    assert res.success is False
    assert 'did not converge' in res.message
    assert '0/' in res.message           # no polish converged


def test_rosenbrock_2d_small_run():
    res = _fit.hybrid_global_optimize(
        _rosenbrock, [(-2.0, 2.0), (-1.0, 3.0)], de_popsize=8,
        de_maxiter=40, do_basinhop=False, polish_top_k=2, random_state=1)
    assert isinstance(res, _fit.GlobalOptResult)
    assert res.success is True
    assert np.allclose(res.x, [1.0, 1.0], atol=1e-3)
    assert res.fun < 1e-6
    assert res.meta['parallel'] == 'none'
    assert res.nfev > 0
    # the same seed gives the same result
    again = _fit.hybrid_global_optimize(
        _rosenbrock, [(-2.0, 2.0), (-1.0, 3.0)], de_popsize=8,
        de_maxiter=40, do_basinhop=False, polish_top_k=2, random_state=1)
    assert np.array_equal(again.x, res.x)
