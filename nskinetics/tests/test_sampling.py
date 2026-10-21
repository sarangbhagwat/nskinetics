# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""ParamRange and the seeded samplers (F2 item 6)."""
import math

import numpy as np
import pytest

from nskinetics.analysis import ParamRange
from nskinetics.analysis import _sampling


def test_param_range_linear_and_log_mapping():
    lin = ParamRange(1, 20)
    assert lin.bounds == (1.0, 20.0) and lin.log is False
    assert lin.from_unit(0.0) == 1.0 and lin.from_unit(1.0) == 20.0
    assert lin.from_unit(0.5) == 10.5
    assert np.allclose(lin.to_unit(lin.from_unit([0.1, 0.7])), [0.1, 0.7])
    log = ParamRange(0.05, 5, log=True)
    assert log.transformed_bounds == (math.log(0.05), math.log(5))
    assert log.from_unit(0.5) == pytest.approx(0.5, rel=1e-12)  # sqrt(lo hi)
    assert log.from_unit(0.0) == 0.05 and log.from_unit(1.0) == 5.0
    assert log.to_unit(0.5) == pytest.approx(0.5, rel=1e-12)
    u = np.linspace(0, 1, 11)
    x = log.from_unit(u)
    assert np.allclose(np.diff(np.log(x)), math.log(100) / 10)
    assert log.contains(0.05) and not log.contains(5.01)
    assert list(log.contains([0.01, 1.0])) == [False, True]


@pytest.mark.parametrize('args, error', [
    ((1, 1), ValueError),
    ((2, 1), ValueError),
    ((0, 1, True), ValueError),
    ((-1, 1, True), ValueError),
    ((0, float('inf')), ValueError),
    (('0', 1), TypeError),
    ((0, 1, 'yes'), TypeError),
])
def test_param_range_validation(args, error):
    with pytest.raises(error):
        ParamRange(*args)


def test_as_ranges():
    ranges = _sampling._as_ranges({'k': ParamRange(0.05, 5, log=True),
                                   'threshold': (1, 20)})
    assert list(ranges) == ['k', 'threshold']
    assert ranges['threshold'] == ParamRange(1.0, 20.0)
    with pytest.raises(ValueError, match="bounds\\['k'\\]"):
        _sampling._as_ranges({'k': (2, 1)})
    with pytest.raises(TypeError):
        _sampling._as_ranges([('k', (0, 1))])
    with pytest.raises(ValueError):
        _sampling._as_ranges({})
    with pytest.raises(TypeError):
        _sampling._as_ranges({'k': 3})


@pytest.mark.parametrize('method', ['sobol', 'lhs'])
def test_unit_sample_is_seeded(method):
    a = _sampling._unit_sample(16, 3, method=method, seed=0)
    b = _sampling._unit_sample(16, 3, method=method, seed=0)
    c = _sampling._unit_sample(16, 3, method=method, seed=1)
    assert a.shape == (16, 3) and np.array_equal(a, b)
    assert not np.array_equal(a, c)
    assert np.all((a >= 0) & (a < 1))
    ss = np.random.SeedSequence(0)
    assert np.array_equal(_sampling._unit_sample(16, 3, method=method,
                                                 seed=ss), a)


def test_latin_hypercube_and_sobol_balance():
    n = 8
    lhs = _sampling._unit_sample(n, 2, method='lhs', seed=3)
    for j in range(2):      # one point per stratum in every dimension
        assert sorted(np.floor(lhs[:, j] * n).astype(int)) == list(range(n))
    sob = _sampling._unit_sample(n, 2, method='sobol', seed=3)
    for j in range(2):      # a power-of-two Sobol' block is balanced too
        assert sorted(np.floor(sob[:, j] * n).astype(int)) == list(range(n))
    # a non-power-of-two n is the head of the next power-of-two block
    head = _sampling._unit_sample(5, 2, method='sobol', seed=3)
    assert np.array_equal(head, sob[:5])


def test_sample_maps_onto_ranges():
    df = _sampling._sample({'k': ParamRange(0.05, 5, log=True),
                            'threshold': (1, 20)}, 64, seed=2)
    assert list(df.columns) == ['k', 'threshold'] and len(df) == 64
    assert df['k'].between(0.05, 5).all()
    assert df['threshold'].between(1, 20).all()
    # uniform in ln k: about half the points below the geometric midpoint
    assert abs((df['k'] < 0.5).mean() - 0.5) <= 1 / 64
    unit = _sampling._unit_sample(64, 2, seed=2)
    assert np.array_equal(df.to_numpy(), _sampling._scale(
        unit, [ParamRange(0.05, 5, log=True), ParamRange(1, 20)]))


def test_gaussian_draws():
    mean = np.array([1.0, -2.0])
    cov = np.array([[0.04, 0.01], [0.01, 0.09]])
    draws = _sampling._gaussian(mean, cov, 20000, seed=0)
    assert draws.shape == (20000, 2)
    assert np.allclose(draws.mean(axis=0), mean, atol=0.01)
    assert np.allclose(np.cov(draws.T), cov, atol=0.005)
    assert np.array_equal(draws, _sampling._gaussian(mean, cov, 20000,
                                                     seed=0))
    # a singular covariance keeps draws on its line x0 = x1
    sing = _sampling._gaussian([0.0, 0.0], [[1.0, 1.0], [1.0, 1.0]], 50,
                               seed=1)
    assert np.allclose(sing[:, 0], sing[:, 1])
    clipped = _sampling._gaussian(mean, cov, 500, seed=0, lo=[0.9, -2.1],
                                  hi=[1.1, -1.9])
    assert clipped[:, 0].min() >= 0.9 and clipped[:, 1].max() <= -1.9
    with pytest.raises(ValueError):
        _sampling._gaussian(mean, np.eye(3), 10)


@pytest.mark.parametrize('seed, error', [(-1, ValueError), (1.5, TypeError),
                                         (True, TypeError), ('0', TypeError)])
def test_seed_validation(seed, error):
    with pytest.raises(error):
        _sampling._seed_sequence(seed)


def test_child_seeds_are_independent_and_do_not_mutate_the_parent():
    ss = _sampling._seed_sequence(5)
    children = [_sampling._child_seed(ss, i) for i in range(2)]
    # The same streams as spawn() on a fresh sequence, without advancing ss.
    spawned = np.random.SeedSequence(5).spawn(2)
    for mine, theirs in zip(children, spawned):
        assert np.array_equal(mine.generate_state(4), theirs.generate_state(4))
    assert ss.n_children_spawned == 0
    a, b = (np.random.default_rng(c).random(3) for c in children)
    assert not np.array_equal(a, b)
    assert _sampling._seed_record(ss) == 5
    assert _sampling._seed_record(children[0]) is None
    fresh = _sampling._seed_sequence(None)
    again = _sampling._seed_sequence(_sampling._seed_record(fresh))
    assert np.array_equal(fresh.generate_state(4), again.generate_state(4))
