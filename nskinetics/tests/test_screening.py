# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Kinetics-level global screening (F2 item 5; spec tests g, h, i, j, k):
pure-numpy Morris/Sobol' cores on analytic functions, then morris/sobol on
inline toys through run_many."""
import math
import warnings

import matplotlib
import numpy as np
import pandas as pd
import pytest

import nskinetics as nsk
from nskinetics.analysis import (At, CallableMetric, Final, MorrisResult,
                                 ParamRange, SobolResult, morris, sobol)
from nskinetics.analysis._sampling import _scale
from nskinetics.analysis.screening import (_morris_design, _morris_indices,
                                           _sobol_design, _sobol_indices,
                                           _sobol_split)

matplotlib.use('Agg')

UNITS = {'time': 'h', 'conc': 'g/L'}

# S => P with an inert parameter z: P(T) = P0 + S0 (1 - e^(-kT)).
INERT_TOY = """
model inert()
  S = 10; P = 0; k = 0.5; z = 1;
  J: S => P; k*S;
end
"""

# dS/dt = k S^2 with S(0) = 1 blows up at t = 1/k: runs with k > 1/T fail.
BLOW_TOY = """
model blow()
  species S; S = 1; k = 0.2; z = 1;
  J: => S; k*S^2;
end
"""


def _model_state(km):
    """Everything A7 says an analysis must leave unchanged."""
    integrator = km._te.getIntegrator()
    return (km._user_parameter_values(),
            {v: km._te[v] for v in km.baseline_variables},
            list(km.events), km._events_compiled,
            {k: integrator.getValue(k) for k in integrator.getSettings()},
            km._te.getSBML())


@pytest.fixture(scope='module')
def inert_model():
    return nsk.KineticModel.from_antimony(INERT_TOY, units=UNITS)


# The base condition sets its own k (replaced by the samples) and P0 = 1
# (kept in every sample).
BASE = nsk.Condition(t_end=5, n_steps=11, parameters={'k': 0.3},
                     initial={'[P]': 1.0}, label='base')
BOUNDS = {'k': ParamRange(0.1, 1.0, log=True), 'z': (0.0, 5.0),
          '[S]': (5.0, 15.0)}


# --- (g) pure cores on analytic functions ------------------------------------
def _ishigami(x, a=7.0, b=0.1):
    return (np.sin(x[:, 0]) + a * np.sin(x[:, 1]) ** 2
            + b * x[:, 2] ** 4 * np.sin(x[:, 0]))


def test_sobol_core_on_ishigami():
    n, d = 2 ** 12, 3
    unit = _sobol_design(d, n, seed=0)
    assert unit.shape == (n * (d + 2), d)
    x = -np.pi + 2 * np.pi * unit
    fa, fb, fab = _sobol_split(_ishigami(x), n, d)
    out = _sobol_indices(fa, fb, fab, n_boot=50, seed=1)
    assert out['S1'].shape == (3, 1) and list(out['n']) == [n]
    assert np.allclose(out['S1'][:, 0], [0.314, 0.442, 0.0], atol=0.05)
    assert np.allclose(out['ST'][:, 0], [0.558, 0.442, 0.244], atol=0.05)
    conf = out['ST_conf'][:, 0]
    assert np.all(conf > 0) and np.all(conf < 0.05)


def test_sobol_design_is_the_saltelli_layout():
    n, d = 8, 3
    unit = _sobol_design(d, n, seed=3)
    a, b, ab = unit[:n], unit[n:2 * n], unit[2 * n:].reshape(d, n, d)
    for i in range(d):
        others = [j for j in range(d) if j != i]
        assert np.array_equal(ab[i][:, i], b[:, i])
        assert np.array_equal(ab[i][:, others], a[:, others])
    assert not np.array_equal(a, b)
    assert np.array_equal(unit, _sobol_design(d, n, seed=3))
    for bad in (0, 3, 100):
        with pytest.raises(ValueError, match='power of 2|>= 2'):
            _sobol_design(d, bad)


def test_sobol_core_on_a_linear_function():
    coef = np.array([1.0, -2.0, 0.5, 0.0])
    ranges = [ParamRange(0, 1), ParamRange(-1, 3), ParamRange(10, 12),
              ParamRange(0, 9)]
    n, d = 2 ** 12, 4
    unit = _sobol_design(d, n, seed=0)
    y = _scale(unit, ranges) @ coef
    out = _sobol_indices(*_sobol_split(y, n, d), n_boot=0)
    var = coef ** 2 * np.array([r.hi - r.lo for r in ranges]) ** 2 / 12
    expected = var / var.sum()
    assert np.allclose(out['S1'][:, 0], expected, atol=0.02)
    assert np.allclose(out['ST'][:, 0], expected, atol=0.02)
    assert out['ST'][3, 0] == 0.0                  # no effect: exactly 0
    assert np.isnan(out['S1_conf']).all()          # n_boot=0


def test_morris_core_on_a_linear_function():
    coef = np.array([3.0, -1.5, 0.0])
    ranges = [ParamRange(0, 2), ParamRange(1, 5), ParamRange(-1, 1)]
    r, d = 12, 3
    unit = _morris_design(d, r, levels=4, seed=0)
    y = _scale(unit, ranges) @ coef
    out = _morris_indices(unit, y, r=r, n_boot=20, seed=1)
    span = np.array([rg.hi - rg.lo for rg in ranges])
    # mu* is |coefficient| x range; every effect is the same.
    assert np.allclose(out['mu_star'][:, 0], np.abs(coef) * span, rtol=1e-12)
    assert np.allclose(out['mu'][:, 0], coef * span, rtol=1e-12)
    assert np.allclose(out['sigma'][:, 0], 0.0, atol=1e-12)
    assert np.allclose(out['mu_star_conf'][:, 0], 0.0, atol=1e-12)
    assert out['ee'].shape == (r, d, 1) and np.all(out['n'] == r)


def test_morris_design_is_one_at_a_time_on_the_grid():
    r, d, levels = 5, 4, 4
    unit = _morris_design(d, r, levels, seed=2)
    assert unit.shape == (r * (d + 1), d)
    grid = np.arange(levels) / (levels - 1)
    assert np.all(np.min(np.abs(unit[..., None] - grid), axis=-1) < 1e-12)
    steps = np.diff(unit.reshape(r, d + 1, d), axis=1)
    moved = steps != 0
    assert np.all(moved.sum(axis=2) == 1) and np.all(moved.sum(axis=1) == 1)
    assert np.allclose(np.abs(steps[moved]), 2 / 3)
    assert np.array_equal(unit, _morris_design(d, r, levels, seed=2))
    assert not np.array_equal(unit, _morris_design(d, r, levels, seed=3))
    with pytest.raises(ValueError, match='even'):
        _morris_design(d, r, levels=5)
    with pytest.raises(ValueError, match='one coordinate'):
        _morris_indices(np.zeros((r * (d + 1), d)), np.zeros(r * (d + 1)),
                        r=r)


def test_cores_drop_non_finite_values():
    r, d = 6, 2
    unit = _morris_design(d, r, seed=0)
    y = unit @ np.array([1.0, 2.0])
    y[0] = np.nan                       # first point of trajectory 0
    out = _morris_indices(unit, y, r=r, n_boot=0)
    # Only the first step of trajectory 0 uses point 0.
    assert out['n'].sum() == r * d - 1
    assert np.allclose(out['mu_star'][:, 0], [1.0, 2.0])
    n = 16
    unit = _sobol_design(d, n, seed=0)
    y = unit @ np.array([1.0, 2.0])
    y[2 * n + 5] = np.inf               # AB_1, row 5
    out = _sobol_indices(*_sobol_split(y, n, d), n_boot=10, seed=0)
    assert list(out['n']) == [n - 1]
    assert np.all(np.isfinite(out['ST']))
    # The cores reject one bootstrap resample too (its spread is undefined).
    with pytest.raises(ValueError, match='n_boot'):
        _sobol_indices(*_sobol_split(y, n, d), n_boot=1)
    with pytest.raises(ValueError, match='n_boot'):
        _morris_indices(_morris_design(d, r, seed=0), np.zeros(r * (d + 1)),
                        r=r, n_boot=1)


# --- (h) the toys ------------------------------------------------------------
def test_morris_on_the_toy_finds_the_inert_parameter(inert_model):
    before = _model_state(inert_model)
    g = morris(inert_model, BASE, BOUNDS, {'P': Final('[P]')}, r=8, seed=0,
               n_boot=50)
    assert isinstance(g, MorrisResult)
    assert g.parameters == ['k', 'z', '[S]'] and g.metric_names == ['P']
    assert g.n_runs == 8 * 4 and g.n_failed == 0
    assert g.mu_star.loc['z', 'P'] == 0.0 and g.sigma.loc['z', 'P'] == 0.0
    assert g.mu_star.loc['k', 'P'] > 1 and g.mu_star.loc['[S]', 'P'] > 1
    assert list(g.rank('P')['parameter'])[-1] == 'z'
    # Each sample was run with its own k and [S] and the base P0 = 1.
    samples = g.samples
    assert list(samples.columns) == ['k', 'z', '[S]', 'status', 'error', 'P']
    expected = 1 + samples['[S]'] * (1 - np.exp(-5 * samples['k']))
    assert np.allclose(samples['P'], expected, rtol=1e-5)
    assert samples['k'].between(0.1, 1.0).all()
    assert set(np.round(np.log(samples['k'] / 0.1) / np.log(10) * 3, 9)) \
        <= {0.0, 1.0, 2.0, 3.0}         # the 4-level grid, in ln k
    # mu* of [S] is (1 - e^(-kT)) x range on average: between the extremes.
    assert 10 * (1 - math.exp(-0.5)) <= g.mu_star.loc['[S]', 'P'] \
        <= 10 * (1 - math.exp(-5))
    assert g.elementary_effects.shape == (8, 3, 1)
    assert _model_state(inert_model) == before


def test_sobol_on_the_toy(inert_model):
    s = sobol(inert_model, BASE, {'k': ParamRange(0.1, 1.0, log=True),
                                  'z': (0.0, 5.0)},
              [Final('[P]', name='P')], n=32, seed=0, n_boot=20)
    assert isinstance(s, SobolResult)
    assert s.n_runs == 32 * 4 and s.n_failed == 0
    assert s.ST.loc['z', 'P'] == 0.0
    assert s.ST.loc['k', 'P'] == pytest.approx(1.0, abs=0.1)
    assert s.S1.loc['k', 'P'] == pytest.approx(1.0, abs=0.1)
    assert s.ST.columns.name == 'metric' and s.ST.index.name == 'parameter'
    assert list(s.rank()['parameter']) == ['k', 'z']


# --- (i) seeds ---------------------------------------------------------------
def test_same_seed_gives_identical_indices(inert_model):
    metrics = {'P': Final('[P]'), 'S': Final('[S]')}
    kw = dict(r=4, n_boot=30)
    a = morris(inert_model, BASE, BOUNDS, metrics, seed=7, **kw)
    b = morris(inert_model, BASE, BOUNDS, metrics, seed=7, **kw)
    for name in MorrisResult._TABLES:
        pd.testing.assert_frame_equal(getattr(a, name), getattr(b, name),
                                      check_exact=True)
    assert np.array_equal(a.unit_design, b.unit_design) and a.seed == 7
    c = morris(inert_model, BASE, BOUNDS, metrics, seed=8, **kw)
    assert not np.array_equal(a.unit_design, c.unit_design)
    # A run without a seed records the entropy that reproduces it.
    d = morris(inert_model, BASE, BOUNDS, metrics, seed=None, **kw)
    e = morris(inert_model, BASE, BOUNDS, metrics, seed=d.seed, **kw)
    assert isinstance(d.seed, int)
    pd.testing.assert_frame_equal(d.mu_star_conf, e.mu_star_conf,
                                  check_exact=True)
    bounds = {'k': BOUNDS['k'], '[S]': BOUNDS['[S]']}
    s1 = sobol(inert_model, BASE, bounds, metrics, n=8, seed=3, n_boot=20)
    s2 = sobol(inert_model, BASE, bounds, metrics, n=8, seed=3, n_boot=20)
    for name in SobolResult._TABLES:
        pd.testing.assert_frame_equal(getattr(s1, name), getattr(s2, name),
                                      check_exact=True)
    s3 = sobol(inert_model, BASE, bounds, metrics, n=8, seed=4, n_boot=20)
    assert not s1.ST.equals(s3.ST)


# --- failures ----------------------------------------------------------------
def test_failed_runs_are_excluded_counted_and_warned():
    km = nsk.KineticModel.from_antimony(BLOW_TOY, units=UNITS)
    cond = nsk.Condition(t_end=2, n_steps=5)
    bounds = {'k': (0.1, 0.7), 'z': (0, 1)}
    with pytest.warns(UserWarning, match='runs failed'):
        s = sobol(km, cond, bounds, {'S': Final('[S]')}, n=16, seed=0,
                  n_boot=10)
    failed = s.samples['status'] != 'ok'
    assert s.n_failed == failed.sum() > 0.05 * s.n_runs
    # The runs with k beyond the blow-up threshold 1/T failed.
    k = s.design['k']
    assert (k[failed] > 0.5).all() and failed[k > 0.52].all()
    assert s.samples.loc[failed, 'error'].str.contains('Error').all()
    assert np.isnan(s.samples.loc[failed, 'S']).all()
    assert np.isfinite(s.ST['S']).all()
    with pytest.warns(UserWarning, match='runs failed'):
        g = morris(km, cond, bounds, {'S': Final('[S]')}, r=6, seed=0,
                   n_boot=10)
    assert g.n_failed > 0
    ee = g.elementary_effects
    assert np.isnan(ee).any() and np.isfinite(ee).any()
    # Below the threshold nothing fails and nothing is warned.
    with warnings.catch_warnings():
        warnings.simplefilter('error', UserWarning)
        ok = morris(km, cond, {'k': (0.1, 0.4)}, {'S': Final('[S]')}, r=3,
                    seed=0, n_boot=0)
    assert ok.n_failed == 0 and np.isnan(ok.mu_star_conf.iloc[0, 0])


def test_non_finite_metric_values_are_counted_and_warned(inert_model):
    # A metric that is NaN on part of the design (here when little substrate
    # is left) is excluded there, like a failed run, but the runs succeed.
    def low_s(res):
        s_end = res['[S]'][-1]
        return s_end if s_end > 1.0 else math.nan

    metric = CallableMetric(low_s, 'h', selections=('[S]',))
    bounds = {'k': ParamRange(0.1, 1.0, log=True), 'z': (0.0, 5.0)}
    cond = nsk.Condition(t_end=5, n_steps=11)
    with pytest.warns(UserWarning, match=r"not finite .*'h' on \d+ of 64"):
        s = sobol(inert_model, cond, bounds, {'h': metric}, n=16, seed=0,
                  n_boot=0)
    assert s.n_failed == 0 and s.n_runs == 64
    nan_rows = int(np.isnan(s.metric_values[:, 0]).sum())
    assert 0.05 * 64 < nan_rows < 64
    assert 0 < s.n_used['h'] < 16
    with pytest.warns(UserWarning, match="'h'"):
        g = morris(inert_model, cond, bounds, {'h': metric}, r=8, seed=0,
                   n_boot=0)
    assert g.n_failed == 0
    used = g.n_used
    assert list(used.columns) == ['h'] and (used['h'] < 8).any()
    assert (used['h'] <= 8).all() and used.dtypes['h'] == np.int64
    # a fully finite metric is not warned about and uses every row
    with warnings.catch_warnings():
        warnings.simplefilter('error', UserWarning)
        ok = sobol(inert_model, cond, bounds, {'P': Final('[P]')}, n=16,
                   seed=0, n_boot=0)
    assert ok.n_used['P'] == 16


def test_ids_read_by_initial_assignments_are_warned_about():
    km = nsk.KineticModel.from_antimony("""
    model ia()
      S_init = 10; S = S_init; P = 0; k = 0.5;
      J: S => P; k*S;
    end
    """, units=UNITS)
    cond = nsk.Condition(t_end=6, n_steps=7)
    with pytest.warns(UserWarning, match=r"init\(S_init\)") as record:
        g = morris(km, cond, {'S_init': (5.0, 20.0), 'k': (0.2, 0.8)},
                   {'P': Final('[P]')}, r=2, seed=0, n_boot=0)
    assert record[0].filename == __file__
    assert g.mu_star.loc['S_init', 'P'] == 0.0      # looks inert
    with warnings.catch_warnings():
        warnings.simplefilter('error', UserWarning)
        fixed = morris(km, cond, {'init(S_init)': (5.0, 20.0)},
                       {'P': Final('[P]')}, r=2, seed=0, n_boot=0)
    # P(6) = S_init (1 - e^-3): the full range moves it by 15 (1 - e^-3)
    assert fixed.mu_star.loc['init(S_init)', 'P'] == pytest.approx(
        15 * (1 - math.exp(-3)), rel=1e-6)


def test_jansen_core_without_variance_or_samples():
    n, d = 8, 2
    const = np.ones(n)
    out = _sobol_indices(const, const, np.ones((d, n)), n_boot=5, seed=0)
    assert np.isnan(out['S1']).all() and np.isnan(out['ST']).all()
    assert out['n'][0] == n
    nan = np.full(n, np.nan)
    out = _sobol_indices(nan, nan, np.full((d, n), np.nan), n_boot=5, seed=0)
    assert out['n'][0] == 0
    assert np.isnan(out['S1']).all() and np.isnan(out['ST']).all()


def test_a_constant_metric_has_nan_sobol_indices_and_is_warned(inert_model):
    # zero variance over the design: the indices are undefined
    const = CallableMetric(lambda res: 1.0, 'one')
    bounds = {'k': BOUNDS['k'], 'z': BOUNDS['z']}
    with pytest.warns(UserWarning, match=r"\['one'\] do not vary"):
        s = sobol(inert_model, BASE, bounds,
                  {'one': const, 'P': Final('[P]')}, n=8, seed=0, n_boot=5)
    for table in ('S1', 'ST', 'S1_conf', 'ST_conf'):
        assert getattr(s, table)['one'].isna().all()
    assert np.isfinite(s.ST['P']).all()


# --- validation --------------------------------------------------------------
@pytest.mark.parametrize('call, kwargs, error, match', [
    ('sobol', dict(n=100), ValueError, 'power of 2'),
    ('sobol', dict(n=1), ValueError, 'n must be'),
    ('morris', dict(levels=3), ValueError, 'even'),
    ('morris', dict(r=1), ValueError, 'r must be'),
    ('morris', dict(r=2.5), TypeError, 'r must be'),
    ('morris', dict(n_boot=-1), ValueError, 'n_boot'),
    # one resample has no spread: every half-width would be NaN
    ('morris', dict(n_boot=1), ValueError, 'n_boot must be 0 .* at least 2'),
    ('sobol', dict(n_boot=1), ValueError, 'n_boot must be 0 .* at least 2'),
    ('morris', dict(seed=-1), ValueError, 'seed'),
    ('morris', dict(workers=0), ValueError, 'workers'),
    ('morris', dict(bounds={'q': (0, 1)}), ValueError, r"^bounds: .*'q'"),
    ('morris', dict(bounds={'k': (1, 0)}), ValueError, 'lo < hi'),
    ('morris', dict(bounds={'k': 3}), TypeError, 'pair'),
    ('morris', dict(bounds={}), ValueError, 'at least one'),
    ('morris', dict(metrics=None), ValueError, 'metrics'),
    ('morris', dict(metrics={}), ValueError, 'at least one'),
    ('morris', dict(metrics={'k': Final('[P]')}), ValueError, 'clash'),
    ('morris', dict(metrics={'status': Final('[P]')}), ValueError, 'clash'),
    ('morris', dict(metrics={'P': Final('[Q]')}), ValueError,
     r'^metrics: .*\[Q\]'),
    ('morris', dict(condition=BASE.replace(reset=False)), ValueError,
     'reset=False'),
    ('morris', dict(condition={'t_end': 5}), TypeError, 'Condition'),
    ('morris', dict(model='km'), TypeError, 'KineticModel'),
    # 'P' (an amount) cannot join the base condition's '[P]'.
    ('sobol', dict(bounds={'P': (1, 2)}, metrics={'titer': Final('[P]')}),
     ValueError, r"^bounds: .*'P'.*same form"),
    # At(t) outside the base condition's span [0, 5] would fail every run
    ('morris', dict(metrics={'late': At('[P]', 100)}), ValueError,
     r'^metrics: .*outside the simulated range \[0.0, 5.0\]'),
    ('sobol', dict(metrics={'late': At('[P]', 100)}), ValueError,
     '^metrics: .*outside'),
])
def test_invalid_arguments_raise_before_any_run(inert_model, monkeypatch,
                                                call, kwargs, error, match):
    runs = []
    monkeypatch.setattr(nsk.KineticModel, 'run',
                        lambda self, *a, **k: runs.append(1))
    args = dict(model=inert_model, condition=BASE, bounds=BOUNDS,
                metrics={'P': Final('[P]')})
    args.update(kwargs)
    fn = {'morris': morris, 'sobol': sobol}[call]
    extra = {k: v for k, v in args.items()
             if k not in ('model', 'condition', 'bounds', 'metrics')}
    with pytest.raises(error, match=match):
        fn(args['model'], args['condition'], args['bounds'], args['metrics'],
           **extra)
    assert runs == []


# --- (j) CSV round trips -----------------------------------------------------
def _assert_same(a, b):
    assert type(a) is type(b)
    for name in ('parameters', 'metric_names', 'ranges', 'status', 'errors',
                 'seed', 'n_boot', 'units', 'condition'):
        assert getattr(a, name) == getattr(b, name), name
    for name in ('unit_design', 'metric_values'):
        assert np.array_equal(getattr(a, name), getattr(b, name),
                              equal_nan=True), name
    for name in a._TABLES:
        pd.testing.assert_frame_equal(getattr(a, name), getattr(b, name),
                                      check_exact=True)
    pd.testing.assert_frame_equal(a.samples, b.samples, check_exact=True)


def test_csv_round_trips_are_exact(inert_model, tmp_path):
    metrics = {'P': Final('[P]'), 'S': Final('[S]')}
    g = morris(inert_model, BASE, BOUNDS, metrics, r=3, seed=1, n_boot=10)
    g.to_csv(tmp_path / 'morris.csv')
    back = MorrisResult.from_csv(tmp_path / 'morris.csv')
    _assert_same(g, back)
    assert (back.r, back.levels) == (g.r, g.levels)
    assert np.array_equal(back.elementary_effects, g.elementary_effects)
    back.to_csv(tmp_path / 'again.csv')
    assert (tmp_path / 'again.csv').read_bytes() == \
        (tmp_path / 'morris.csv').read_bytes()
    # quantities named after the fields they hold
    assert set(pd.read_csv(tmp_path / 'morris.csv')['quantity']) == {
        'meta', 'range', 'unit_design', 'metric_value', 'status', 'error',
        'table'}
    s = sobol(inert_model, BASE, {'k': BOUNDS['k'], 'z': BOUNDS['z']},
              metrics, n=4, seed=2, n_boot=0)
    s.to_csv(tmp_path / 'sobol.csv')
    back = SobolResult.from_csv(tmp_path / 'sobol.csv')
    _assert_same(s, back)
    assert back.n == 4 and np.isnan(back.S1_conf.to_numpy()).all()
    with pytest.raises(ValueError, match='SobolResult'):
        SobolResult.from_csv(tmp_path / 'morris.csv')
    with pytest.raises(ValueError, match='MorrisResult'):
        MorrisResult.from_csv(tmp_path / 'sobol.csv')


def test_from_csv_rejects_incomplete_files(inert_model, tmp_path):
    # A truncated or edited file is refused, never filled with
    # uninitialized memory (or a KeyError for a missing meta key).
    g = morris(inert_model, BASE, BOUNDS, {'P': Final('[P]')}, r=3, seed=1,
               n_boot=0)
    path = tmp_path / 'morris.csv'
    g.to_csv(path)
    lines = path.read_text().splitlines(keepends=True)
    for prefix, match in (('unit_design,', 'unit_design: expected'),
                          ('metric_value,', 'metric_value: expected'),
                          ('status,', 'status: expected'),
                          ('range,', 'range: expected'),
                          ('table,mu_star,', 'mu_star: expected'),
                          ('meta,levels,', r"meta keys \['levels'\]")):
        first = next(line for line in lines if line.startswith(prefix))
        cut = tmp_path / 'cut.csv'
        cut.write_text(''.join(line for line in lines if line != first))
        with pytest.raises(ValueError, match=f'incomplete.*{match}'):
            MorrisResult.from_csv(cut)


def test_csv_round_trip_with_failures(tmp_path):
    km = nsk.KineticModel.from_antimony(BLOW_TOY, units=UNITS)
    with pytest.warns(UserWarning):
        s = sobol(km, nsk.Condition(t_end=2, n_steps=5), {'k': (0.1, 0.7)},
                  {'S': Final('[S]')}, n=4, seed=0, n_boot=5)
    s.to_csv(tmp_path / 'fail.csv')
    _assert_same(s, SobolResult.from_csv(tmp_path / 'fail.csv'))


# --- (k) plots ---------------------------------------------------------------
def test_plot_smoke(inert_model):
    import matplotlib.pyplot as plt
    metrics = {'P': Final('[P]'), 'S': Final('[S]')}
    g = morris(inert_model, BASE, BOUNDS, metrics, r=3, seed=0, n_boot=5)
    fig, ax = g.plot('S')
    assert ax.get_title() == 'S' and ax.get_xlabel() == 'mu*'
    assert {t.get_text() for t in ax.texts} <= set(g.parameters)
    plt.close(fig)
    s = sobol(inert_model, BASE, {'k': BOUNDS['k'], 'z': BOUNDS['z']},
              metrics, n=4, seed=0, n_boot=5)
    fig, ax = s.plot()
    assert ax.get_title() == 'P' and len(ax.patches) == 4
    assert [t.get_text() for t in ax.get_xticklabels()] == ['k', 'z']
    plt.close(fig)
    for result in (g, s):
        fig0, ax0 = plt.subplots()
        assert result.plot('P', ax=ax0) == (fig0, ax0)
    plt.close('all')
    with pytest.raises(KeyError, match='not a metric'):
        s.plot('Q')
