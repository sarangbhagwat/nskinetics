# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests of nskinetics.analysis.calibration (datasets, problems, fitting,
FitResult statistics, CSV round trip, apply, Laplace samples, prediction
bands and plots)."""
import dataclasses
import math
import re
import warnings

import matplotlib
import numpy as np
import pandas as pd
import pytest
import tellurium as te
from scipy.optimize import least_squares
from scipy.stats import chi2 as chi2_distribution

import nskinetics as nsk
from nskinetics.analysis import (CalibrationProblem, Dataset, Experiment,
                                 FitParameter, FitResult, Observable,
                                 fit_parameters, profile_likelihood)
from nskinetics.analysis import calibration as cal
from nskinetics.exceptions import FittingError, KineticSimulationError

matplotlib.use('Agg')

UNITS = {'time': 'h', 'conc': 'g/L'}
Z = 1.959963984540054

# S => P: S = S0 e^(-k t), so dS/dk = -t S and, with sigma known,
# SE_k = 1 / sqrt(sum (t_i S_i / sigma)^2).
DECAY_TOY = """
model decay()
  S = 10; P = 0; k = 0.5;
  J: S => P; k*S;
end
"""

# Only the product k1*k2 is identifiable.
PRODUCT_TOY = """
model prod()
  S = 10; k1 = 0.5; k2 = 1;
  J: S => ; k1*k2*S;
end
"""

# dx/dt = k x^2 blows up at t = 1/k: for k > 1/0.9 a run to 0.9 fails.
BLOW_TOY = """
model blow()
  x = 1; k = 0.5;
  J: => x; k*x^2;
end
"""

FEEDSPIKE_TOY = """
model fb()
  compartment env = 1; species S in env, P in env;
  S = 50; P = 0; k = 0.5; threshold = 5; target = 50; feed_conc = 500;
  n_max = 3; n_spk = 0; last_vol = 0; tot_vol = 0;
  J: S => P; k*S*env;
end
"""

T = np.linspace(0.5, 6.0, 12)


def _model(text):
    return nsk.KineticModel(te.loadAntimonyModel(text), units=UNITS)


def _feedspike_model():
    km = _model(FEEDSPIKE_TOY)
    spike = nsk.FeedSpike('S', 'S <= threshold', 'target', 'feed_conc',
                          volume_var='env', max_count='n_max',
                          count_var='n_spk', last_vol_var='last_vol',
                          tot_vol_var='tot_vol')
    for event in spike.expand():
        km.add_event(event)
    km.compile_events()
    return km


def _model_state(km):
    """Everything A7 says an analysis must leave unchanged."""
    integrator = km._te.getIntegrator()
    return (km._user_parameter_values(),
            {v: km._te[v] for v in km.baseline_variables},
            list(km.events), km._events_compiled,
            {k: integrator.getValue(k) for k in integrator.getSettings()},
            km._te.getSBML())


def _decay_data(k=0.5, S0=10.0, sigma=None, seed=None, times=T):
    S = S0 * np.exp(-k * times)
    if seed is not None:
        S = S + sigma * np.random.default_rng(seed).standard_normal(S.size)
    return pd.DataFrame({'observable': 'S', 'time': times, 'value': S})


def _decay_problem(km=None, data=None, sigma=0.1, params=None, **kwargs):
    km = km if km is not None else _model(DECAY_TOY)
    data = data if data is not None else _decay_data()
    exp = Experiment('E', nsk.Condition(t_end=6), data)
    params = params or [FitParameter('k', (0.01, 10), initial=0.2)]
    return CalibrationProblem(km, [exp], [Observable('S', '[S]', sigma=sigma)],
                              params, **kwargs)


@pytest.fixture(scope='module')
def decay_fit():
    return fit_parameters(_decay_problem(), n_starts=3, seed=0)


@pytest.fixture(scope='module')
def blow_problem():
    t = np.linspace(0.1, 0.9, 9)
    data = pd.DataFrame({'observable': 'x', 'time': t,
                         'value': 1 / (1 - 0.5 * t)})
    exp = Experiment('E', nsk.Condition(t_end=0.9), data)
    return CalibrationProblem(_model(BLOW_TOY), [exp],
                              [Observable('x', 'x', sigma=0.05)],
                              [FitParameter('k', (0.1, 3), initial=0.2)])


@pytest.fixture(scope='module')
def blow_fit(blow_problem):
    return fit_parameters(blow_problem, n_starts=4, seed=0)


# --- (a) recovery and determinism --------------------------------------------
def test_noise_free_decay_recovers_k_deterministically(decay_fit):
    fit = decay_fit
    assert fit.names == ['k']
    assert abs(fit.estimates['k'] / 0.5 - 1) < 1e-6
    assert fit.success and list(fit.flags['flag']) == []
    assert fit.n_data == 12 and fit.n_free == 1 and fit.dof == 11
    # the FitParameter bounds exactly (exp(ln 0.01) is not 0.01)
    assert fit.lower[0] == 0.01 and fit.upper[0] == 10.0
    assert len(fit.starts) == 3
    assert fit.starts['best'].sum() == 1
    # start 0 is the user's initial value
    assert math.isclose(fit.starts['initial:k'].iloc[0], 0.2, rel_tol=1e-12)
    again = fit_parameters(_decay_problem(), n_starts=3, seed=0)
    assert np.array_equal(again.x, fit.x)
    pd.testing.assert_frame_equal(again.starts, fit.starts)


def test_start_zero_defaults_to_the_value_the_run_starts_from(run_counter):
    # the model's value (k = 0.5), read on first use of x0, not at
    # construction (validation runs nothing)
    prob = _decay_problem(params=[FitParameter('k', (0.01, 1))])
    assert run_counter == []
    assert np.allclose(prob.x0, [math.log(0.5)])
    # a model value outside the bounds: the midpoint in the fitting space
    prob = _decay_problem(params=[FitParameter('k', (1, 100))])
    assert np.allclose(prob.x0, [math.log(10)])
    # a start on a bound is kept (no margin: the optimizer leaves it)
    prob = _decay_problem(params=[FitParameter('k', (0.5, 2), initial=0.5)])
    lo, hi = prob.bounds
    assert prob.x0[0] == lo[0] == math.log(0.5)
    # a local '[S]' starts at each condition's recipe (else the model's
    # value); an observable scale at the midpoint of its bounds
    exps = [Experiment('E1', nsk.Condition(t_end=6, initial={'[S]': 20}),
                       _decay_data()),
            Experiment('E2', nsk.Condition(t_end=6), _decay_data())]
    prob = CalibrationProblem(
        _model(DECAY_TOY), exps, [Observable('S', '[S]', scale='a',
                                             sigma=0.1)],
        [FitParameter('k', (0.01, 10)),
         FitParameter('[S]', (1, 100), scope='local'),
         FitParameter('a', (0.1, 10), log=False)])
    assert prob.names == ['k', '[S]@E1', '[S]@E2', 'a']
    assert np.allclose(prob.x0, [math.log(0.5), math.log(20), math.log(10),
                                 5.05])
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert np.allclose(fit.starts.filter(like='initial:').iloc[0],
                       [0.5, 20, 10, 5.05])
    # the start-value run fails (k = 2 blows up before t_end): the midpoint
    data = pd.DataFrame({'observable': 'x', 'time': [0.1, 0.6],
                         'value': [1.1, 1.6]})
    exp = Experiment('E', nsk.Condition(t_end=0.9, parameters={'k': 2}),
                     data)
    prob = CalibrationProblem(_model(BLOW_TOY), [exp],
                              [Observable('x', 'x', sigma=0.05)],
                              [FitParameter('x', (0.5, 8), scope='local')])
    assert np.allclose(prob.x0, [math.log(2)])     # not ln(1), the model's


def test_a_start_near_one_in_log_space_does_not_stall():
    # the log-midpoint of (0.1, 10) is ln 1 up to round-off (about 2e-16);
    # SciPy's trf scales its first trust radius with the start's norm, so
    # it used to stop there after 2 evaluations, reporting convergence
    data = _decay_data()
    data['value'] *= 2.5
    exp = Experiment('E', nsk.Condition(t_end=6), data)
    prob = CalibrationProblem(_model(DECAY_TOY), [exp],
                              [Observable('S', '[S]', scale='a', sigma=0.1)],
                              [FitParameter('a', (0.1, 10))])
    assert prob.x0[0] != 0 and abs(prob.x0[0]) < 1e-12
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert abs(fit.estimates['a'] / 2.5 - 1) < 1e-6
    assert fit.starts['nfev'].iloc[0] > 2
    assert fit.starts['initial:a'].iloc[0] == pytest.approx(1, rel=1e-12)


@pytest.mark.parametrize('b0', [0.0, -5e7])
def test_wide_linear_bounds_do_not_loosen_the_stopping_test(b0):
    # trf stops when |step| < xtol (xtol + |u|): a shift of the optimizer's
    # coordinates that grew with the bounds (+-1e8 here) made it stop after
    # one step (chi2 251) and report convergence
    data = _decay_data()
    data['value'] += 0.37
    exp = Experiment('E', nsk.Condition(t_end=6), data)
    prob = CalibrationProblem(
        _model(DECAY_TOY), [exp],
        [Observable('S', '[S]', sigma=0.1, offset='b')],
        [FitParameter('k', (0.01, 10), initial=0.2),
         FitParameter('b', (-1e8, 1e8), log=False, initial=b0)])
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert fit.estimates['k'] == pytest.approx(0.5, rel=1e-6)
    assert fit.estimates['b'] == pytest.approx(0.37, abs=1e-6)
    assert fit.starts['nfev'].iloc[0] > 2
    assert fit.starts['initial:b'].iloc[0] == b0
    # the rank test does not depend on the width of the bounds either
    assert fit.numerical_rank == 2
    assert not {'rank_deficient', 'ill_conditioned'} & set(fit.flags.flag)


def test_rank_test_does_not_depend_on_the_width_of_linear_bounds():
    # the column scale of a linear parameter was max(|p|, 1e-3 (hi - lo)):
    # with an offset near 0 on +-1e9 it dominated, and a well-determined
    # k was reported rank deficient (se inf, se of b 4x too small)
    data = _decay_data(sigma=0.1, seed=0)
    data['value'] += 0.37
    fits = []
    for width in (1.0, 1e9):
        prob = CalibrationProblem(
            _model(DECAY_TOY),
            [Experiment('E', nsk.Condition(t_end=6), data)],
            [Observable('S', '[S]', sigma=0.1, offset='b')],
            [FitParameter('k', (0.01, 10), initial=0.2),
             FitParameter('b', (-width, width), log=False, initial=0.0)])
        fits.append(fit_parameters(prob, n_starts=1, seed=0))
    a, b = fits
    for fit in fits:
        assert fit.numerical_rank == 2
        assert np.all(np.isfinite(fit.se))
    # the two optimizations stop a few 1e-6 apart (xtol vs ftol)
    assert np.allclose(a.se, b.se, rtol=1e-4, atol=0)
    assert a.condition_number == pytest.approx(b.condition_number, rel=1e-4)
    assert list(a.flags.flag) == list(b.flags.flag) == ['high_correlation']
    # closed form: Gauss-Newton SEs of (ln k, b) from the analytic Jacobian
    # of 10 e^(-k t) + b at the least-squares optimum
    t, y = data['time'].to_numpy(), data['value'].to_numpy()
    x = least_squares(lambda x: (10 * np.exp(-np.exp(x[0]) * t) + x[1] - y)
                      / 0.1, [np.log(0.5), 0.0], xtol=1e-15, ftol=1e-15,
                      gtol=1e-15).x
    kh = np.exp(x[0])
    jac = np.column_stack([-kh * t * 10 * np.exp(-kh * t),
                           np.ones_like(t)]) / 0.1
    expected = np.sqrt(np.diag(np.linalg.inv(jac.T @ jac)))
    assert np.allclose([a.se['k'], a.se['b']], expected, rtol=1e-4, atol=0)
    # the same Jacobian with wider bounds: identical results
    wide = dataclasses.replace(a, lower=np.array([0.01, -1e9]),
                               upper=np.array([10.0, 1e9]))
    assert wide.numerical_rank == 2
    assert np.array_equal(wide.se, a.se)
    assert wide.condition_number == a.condition_number
    assert list(wide.flags.flag) == list(a.flags.flag)


def test_an_interior_initial_is_used_verbatim_on_a_wide_linear_range():
    # no margin: 0.2 is not moved to 0.1% of the range (100) inside 0
    prob = _decay_problem(params=[FitParameter('k', (0, 1e5), log=False,
                                               initial=0.2)])
    assert prob.x0[0] == 0.2
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert fit.starts['initial:k'].iloc[0] == 0.2
    assert fit.estimates['k'] == pytest.approx(0.5, rel=1e-6)


def test_local_initial_mapping_accepts_none_and_foreign_experiments():
    # None and an experiment not in the problem both leave the default
    # (here E2's own condition value)
    p = FitParameter('[S]', (1, 100), scope='local',
                     initial={'E1': 20, 'E2': None, 'elsewhere': 7})
    assert dict(p.initial) == {'E1': 20.0, 'E2': None, 'elsewhere': 7.0}
    assert p.initial_for('E2') is None and p.initial_for('E3') is None
    exps = [Experiment(name, nsk.Condition(t_end=6, initial={'[S]': 30}),
                       _decay_data()) for name in ('E1', 'E2')]
    prob = CalibrationProblem(_model(DECAY_TOY), exps,
                              [Observable('S', '[S]', sigma=0.1)],
                              [FitParameter('k', (0.01, 10), initial=0.2), p])
    assert prob.names == ['k', '[S]@E1', '[S]@E2']
    assert np.allclose(prob.x0, [math.log(0.2), math.log(20),
                                 math.log(30)])


def test_fit_parameters_hash_pickle_and_copy():
    import copy
    import pickle
    kwargs = dict(scope='local', initial={'E1': 20, 'E2': None})
    p = FitParameter('[S]', (1, 100), **kwargs)
    q = FitParameter('[S]', (1, 100), **kwargs)
    assert p == q and hash(p) == hash(q)
    assert hash(p) != hash(FitParameter('[S]', (1, 100), scope='local'))
    assert len({p, q, FitParameter('k', (0.01, 10), initial=0.2)}) == 2
    for back in (pickle.loads(pickle.dumps(p)), copy.deepcopy(p),
                 copy.copy(p)):
        assert back == p and hash(back) == hash(p)
        with pytest.raises(TypeError):
            back.initial['E1'] = 5          # still read-only
    g = FitParameter('k', (0.01, 10), log=False, initial=0.2)
    assert pickle.loads(pickle.dumps(g)) == g and copy.deepcopy(g) == g


# --- (b) Wald SE, CI, scaling and information criteria -----------------------
def test_wald_se_matches_closed_form_with_seeded_noise():
    sigma = 0.1
    data = _decay_data(sigma=sigma, seed=0)
    fit = fit_parameters(_decay_problem(data=data, sigma=sigma), n_starts=2,
                         seed=0)
    k = fit.estimates['k']
    S = 10 * np.exp(-k * T)
    se_k = 1 / math.sqrt(np.sum((T * S / sigma) ** 2))
    assert math.isclose(fit.se_natural['k'], se_k, rel_tol=1e-3)
    # log-space SE (the fitting space) is SE_k / k; the CI is asymmetric
    assert math.isclose(fit.se['k'], se_k / k, rel_tol=1e-3)
    lo, hi = fit.ci.loc['k']
    assert math.isclose(lo, k * math.exp(-Z * fit.se['k']), rel_tol=1e-12)
    assert math.isclose(hi, k * math.exp(Z * fit.se['k']), rel_tol=1e-12)
    assert hi - k > k - lo
    # chi2 is the plain weighted SSR; statistics follow their formulas
    r = fit.residuals['weighted_residual'].to_numpy()
    assert math.isclose(fit.chi2, np.sum(r ** 2), rel_tol=1e-12)
    model = fit.residuals['model'].to_numpy()
    assert np.allclose(model, S, rtol=1e-6)
    assert np.allclose(r, (model - data['value']) / sigma, rtol=1e-9)
    n, kf = 12, 1
    assert math.isclose(fit.reduced_chi2, fit.chi2 / (n - kf))
    assert math.isclose(fit.aic, fit.chi2 + 2 * kf)
    assert math.isclose(fit.bic, fit.chi2 + kf * math.log(n))
    assert math.isclose(fit.aicc, fit.aic + 2 * kf * (kf + 1) / (n - kf - 1))
    # absolute_sigma=False: same estimate, covariance scaled by reduced chi2
    rel = fit_parameters(_decay_problem(data=data, sigma=sigma,
                                        absolute_sigma=False),
                         n_starts=2, seed=0)
    assert np.allclose(rel.x, fit.x, rtol=1e-10)
    assert math.isclose(rel.se['k'], fit.se['k'] * math.sqrt(
        fit.reduced_chi2), rel_tol=1e-6)
    # the estimated noise level counts: K = n_free + 1 (Hurvich-Tsai AICc)
    kr = kf + 1
    base = n * math.log(rel.chi2 / n)
    assert math.isclose(rel.aic, base + 2 * kr)
    assert math.isclose(rel.bic, base + kr * math.log(n))
    assert math.isclose(rel.aicc, rel.aic + 2 * kr * (kr + 1) / (n - kr - 1))
    # = n ln(sigma^2) + n (n + k) / (n - k - 2), up to a constant
    assert math.isclose(rel.aicc - (base + n * (n + kf) / (n - kf - 2)),
                        -n, rel_tol=1e-12)


def test_robust_loss_reports_plain_chi2_and_resists_an_outlier():
    data = _decay_data()
    data.loc[5, 'value'] += 3.0                 # a 30-sigma outlier
    lin = fit_parameters(_decay_problem(data=data), n_starts=1, seed=0)
    rob = fit_parameters(_decay_problem(data=data, loss='soft_l1'),
                         n_starts=1, seed=0)
    assert abs(rob.estimates['k'] - 0.5) < abs(lin.estimates['k'] - 0.5)
    r = rob.residuals['weighted_residual'].to_numpy()
    assert math.isclose(rob.chi2, np.sum(r ** 2), rel_tol=1e-12)
    assert rob.cost < 0.5 * rob.chi2         # soft_l1 rho(z) < z
    rho = 2 * (np.sqrt(1 + r ** 2) - 1)
    assert math.isclose(rob.cost, 0.5 * np.sum(rho), rel_tol=1e-9)
    assert rob.loss == 'soft_l1'


def test_relative_sigmas_need_the_linear_loss_and_ignore_the_sigma_scale():
    # relative weights: rescaling every sigma changes nothing (linear loss)
    data = _decay_data(sigma=0.1, seed=0)
    a, b = (fit_parameters(_decay_problem(data=data, sigma=sigma,
                                          absolute_sigma=False),
                           n_starts=1, seed=0) for sigma in (1.0, 0.01))
    assert np.allclose(a.x, b.x, rtol=1e-8)
    assert np.allclose(a.se, b.se, rtol=1e-6)
    assert np.allclose(a.ci, b.ci, rtol=1e-6)
    assert math.isclose(a.chi2 * 1e4, b.chi2, rel_tol=1e-8)
    # a robust loss needs absolute sigmas; its no-sigma hint says so
    for loss in ('soft_l1', 'huber', 'cauchy'):
        with pytest.raises(ValueError, match='needs absolute sigmas'):
            _decay_problem(loss=loss, absolute_sigma=False)
    with pytest.raises(ValueError, match='no usable sigma') as err:
        CalibrationProblem(_model(DECAY_TOY), [_exp()],
                           [Observable('S', '[S]')], _K, loss='soft_l1')
    assert 'absolute_sigma=False' not in str(err.value)
    with pytest.raises(ValueError, match='absolute_sigma=False'):
        CalibrationProblem(_model(DECAY_TOY), [_exp()],
                           [Observable('S', '[S]')], _K)


def test_array_theta_outside_the_bounds_warns():
    prob = _decay_problem()                 # k in (0.01, 10): ln k in bounds
    with pytest.warns(UserWarning, match=r"theta lies outside the bounds for "
                      r"\['k'\]: an array is read in the fitting space") as w:
        r = prob.residuals(np.array([12.0]))    # natural values: exp(12)
    assert w[0].filename == __file__
    for call in (prob.cost, prob.to_natural):
        with pytest.warns(UserWarning, match='fitting space') as w:
            call([12.0])
        assert w[0].filename == __file__
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        inside = prob.residuals(np.array([math.log(0.5)]))
        assert np.array_equal(inside, prob.residuals({'k': 0.5}))
        # natural values by name: no warning, and k = 12 rather than e^12
        assert not np.allclose(prob.residuals({'k': 12.0}), r)


def test_initial_sets_start_zero_from_a_previous_fit(run_counter,
                                                     two_parameter_fit):
    data = _decay_data(sigma=0.1, seed=0)
    data.loc[5, 'value'] += 3.0                 # a 30-sigma outlier
    robust = fit_parameters(_decay_problem(data=data, loss='soft_l1'),
                            n_starts=2, seed=0)
    linear = _decay_problem(data=data)
    refit = fit_parameters(linear, initial=robust, n_starts=1, seed=0)
    assert math.isclose(refit.starts['initial:k'].iloc[0],
                        robust.estimates['k'], rel_tol=1e-12)
    full = fit_parameters(linear, n_starts=3, seed=0)
    assert math.isclose(refit.estimates['k'], full.estimates['k'],
                        rel_tol=1e-6)
    # natural values by name, or an array in the fitting space
    for initial in ({'k': 0.3}, np.array([math.log(0.3)])):
        fit = fit_parameters(linear, initial=initial, n_starts=2, seed=0,
                             max_nfev=1)
        assert math.isclose(fit.starts['initial:k'].iloc[0], 0.3,
                            rel_tol=1e-12)
        assert fit.starts['initial:k'].iloc[1] != 0.3   # sampled starts
    # clipped into the bounds like every start (a start on a bound is kept)
    fit = fit_parameters(linear, initial={'k': 0.01}, n_starts=1, seed=0,
                         max_nfev=1)
    lo, hi = linear.bounds
    assert math.isclose(math.log(fit.starts['initial:k'].iloc[0]), lo[0],
                        rel_tol=1e-12)
    with pytest.warns(UserWarning, match='initial lies outside the bounds'):
        fit_parameters(linear, initial=[12.0], n_starts=1, max_nfev=1)
    # a fit of a problem with other free parameters: before any run
    n_runs = len(run_counter)
    with pytest.raises(ValueError, match=r"initial must give every free "
                       r"parameter \['k'\]; unknown \['\[S\]'\]"):
        fit_parameters(linear, initial=two_parameter_fit)
    assert len(run_counter) == n_runs


@pytest.mark.parametrize('loss', cal.LOSSES)
def test_cost_matches_least_squares_for_every_loss(loss):
    prob = _decay_problem(loss=loss, f_scale=2.5)
    theta = {'k': 0.45}
    r = prob.residuals(theta)
    z = np.abs(r) / 2.5
    assert z.min() < 1 < z.max()              # both branches of huber
    ref = least_squares(lambda x: r, [0.0], loss=loss, f_scale=2.5,
                        max_nfev=1).cost
    assert math.isclose(prob.cost(theta), ref, rel_tol=1e-12)


def test_cached_residuals_cannot_be_modified_by_the_optimizer():
    """least_squares may rescale the array fun returns in place (SciPy 1.15
    does for robust losses); the evaluator's cache must not change."""
    prob = _decay_problem()
    prob._ready()
    ev = cal._Evaluator(prob._private, prob._spec)
    x = np.array([math.log(0.4)])
    r0 = ev.evaluate(x)[0].copy()
    r = ev.residuals(x)
    r *= 1e3
    assert np.array_equal(ev.evaluate(x)[0], r0)
    assert np.array_equal(ev.fd_jacobian(x).f0, r0)
    assert np.array_equal(ev.residuals(x), r0)
    with pytest.raises(ValueError, match='read-only'):
        ev.evaluate(x)[0][0] = 0.0
    assert np.array_equal(prob.residuals(x), r0)


def test_observable_scale_is_estimated_as_a_virtual_parameter():
    data = _decay_data()
    data['value'] *= 2.5
    km = _model(DECAY_TOY)
    exp = Experiment('E', nsk.Condition(t_end=6), data)
    prob = CalibrationProblem(
        km, [exp], [Observable('S', '[S]', scale='a', sigma=0.1)],
        [FitParameter('k', (0.01, 10), initial=0.3),
         FitParameter('a', (0.5, 10), initial=1.0)])
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert fit.kinds == ['model', 'observable']
    assert np.allclose(fit.estimates.to_numpy(), [0.5, 2.5], rtol=1e-6)
    assert fit.fd_flags.loc['a', 'scheme'] == 'central'
    # the scale is returned, not applied
    km2 = _model(DECAY_TOY)
    rest = fit.apply(km2)
    assert list(rest) == ['a'] and math.isclose(rest['a'], 2.5, rel_tol=1e-6)
    assert math.isclose(km2._te['k'], fit.estimates['k'])


def test_local_observable_scale_names_and_residuals():
    km = _model(DECAY_TOY)
    scaled = _decay_data()
    scaled['value'] *= 3.0
    exps = [Experiment('E1', nsk.Condition(t_end=6), scaled),
            Experiment('E2', nsk.Condition(t_end=6), _decay_data()),
            Experiment('E3', nsk.Condition(t_end=6), _decay_data(),
                       observables=[Observable('S', '[S]', sigma=0.1)])]
    prob = CalibrationProblem(
        km, exps, [Observable('S', '[S]', scale='a', sigma=0.1)],
        [FitParameter('k', (0.01, 10)),
         FitParameter('a', (0.1, 10), scope='local')])
    # E3 replaces the shared observable by one without a scale
    assert prob.names == ['k', 'a@E1', 'a@E2']
    r = prob.residuals({'k': 0.5, 'a@E1': 3.0, 'a@E2': 1.0})
    assert np.allclose(r, 0.0, atol=1e-6)
    r = prob.residuals({'k': 0.5, 'a@E1': 3.0, 'a@E2': 2.0})
    assert np.allclose(r[:12], 0.0, atol=1e-6)
    assert np.allclose(r[12:24], _decay_data()['value'] / 0.1, rtol=1e-6)
    assert np.allclose(r[24:], 0.0, atol=1e-6)


def test_noise_models_and_sigma_column_override():
    y = np.array([0.0, 2.0, 5.0, 10.0])
    t = np.array([1.0, 2.0, 3.0, 4.0])
    rows = pd.DataFrame({'experiment': 'E', 'observable': 'P', 'time': t,
                         'value': y, 'sigma': [np.nan, np.nan, 0.7, np.nan]})
    km = _model(DECAY_TOY)

    def sigmas(obs, floor=cal.RELATIVE_FLOOR):
        exp = Experiment('E', nsk.Condition(t_end=6), Dataset(rows))
        prob = CalibrationProblem(km, [exp], [obs],
                                  [FitParameter('k', (0.01, 10))],
                                  relative_floor=floor)
        return prob.data['sigma'].to_numpy()

    prop = 0.05 * np.maximum(y, cal.RELATIVE_FLOOR * 10.0)
    expected = prop.copy()
    expected[2] = 0.7
    got = sigmas(Observable('P', '[P]', noise='proportional', rel_sigma=0.05))
    assert np.allclose(got, expected, rtol=1e-15)
    prop = 0.05 * np.maximum(y, 0.1 * 10.0)
    expected = np.sqrt(0.2 ** 2 + prop ** 2)
    expected[2] = 0.7
    got = sigmas(Observable('P', '[P]', noise='combined', sigma=0.2,
                            rel_sigma=0.05), floor=0.1)
    assert np.allclose(got, expected, rtol=1e-15)
    got = sigmas(Observable('P', '[P]', sigma=0.3))
    assert np.allclose(got, [0.3, 0.3, 0.7, 0.3])
    # the floor is relative to the largest |value| of that observable in
    # that experiment only
    t3 = [1.0, 2.0, 3.0]
    data = Dataset(pd.DataFrame({
        'experiment': ['A'] * 3 + ['B'] * 6,
        'observable': ['P'] * 6 + ['S'] * 3, 'time': t3 * 3,
        'value': [0, 1, 10, 0, 100, 1000, 5, 50, 500]}))
    exps = [Experiment(e, nsk.Condition(t_end=6), data[e]) for e in 'AB']
    obs = [Observable(n, f'[{n}]', noise='proportional', rel_sigma=0.1)
           for n in 'PS']
    frame = CalibrationProblem(km, exps, obs, [FitParameter('k', (0.01, 10))],
                               relative_floor=0.1).data

    def block(e, o):
        rows = frame[(frame['experiment'] == e) & (frame['observable'] == o)]
        return rows.sort_values('time')['sigma'].to_numpy()

    assert np.allclose(block('A', 'P'), [0.1, 0.1, 1.0], rtol=1e-15)
    assert np.allclose(block('B', 'P'), [10, 10, 100], rtol=1e-15)
    assert np.allclose(block('B', 'S'), [5, 5, 50], rtol=1e-15)


def test_observable_values_on_a_simulation_result():
    km = _model(DECAY_TOY)
    res = km.copy().run(nsk.Condition(t_end=6, n_steps=6,
                                      selections=('[S]',)))
    obs = Observable('S', '[S]', scale='a', offset=1.0, sigma=0.1)
    assert obs.parameters == ('a',)
    assert np.allclose(obs.values(res, {'a': 2}), 2 * res['[S]'] + 1)
    with pytest.raises(ValueError, match="'a'"):
        obs.values(res)                       # a named scale needs its value
    with pytest.raises(KeyError):
        Observable('P', '[P]', sigma=0.1).values(res)     # not recorded


def test_transform_is_applied_before_scale_and_offset():
    data = _decay_data()
    data['value'] = 2 * np.log(data['value']) + 1
    km = _model(DECAY_TOY)
    exp = Experiment('E', nsk.Condition(t_end=6), data)
    prob = CalibrationProblem(
        km, [exp], [Observable('S', '[S]', scale=2.0, offset=1.0,
                               transform=np.log, sigma=0.01)],
        [FitParameter('k', (0.01, 10), initial=0.3)])
    assert np.allclose(prob.residuals({'k': 0.5}), 0.0, atol=1e-6)
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert math.isclose(fit.estimates['k'], 0.5, rel_tol=1e-6)


def _log_above_detection(x):
    # a transform with a domain: NaN (not an exception) outside it
    with np.errstate(all='ignore'):
        return np.where(x > 1e-3, np.log(x), np.nan)


def test_a_transform_returning_nan_is_penalized_and_counted():
    data = _decay_data()
    data['value'] = np.log(data['value'])
    prob = CalibrationProblem(
        _model(DECAY_TOY), [Experiment('E', nsk.Condition(t_end=6), data)],
        [Observable('S', '[S]', transform=_log_above_detection,
                    sigma=0.01)],
        [FitParameter('k', (0.01, 10), initial=0.3)])
    assert np.allclose(prob.residuals({'k': 0.5}), 0.0, atol=1e-6)
    n = prob.n_failed_evals
    # S(6) = 10 exp(-30) is below the domain at k = 5
    r = prob.residuals({'k': 5.0})
    assert np.all(r == cal.PENALTY_RESIDUAL)
    assert prob.n_failed_evals == n + 1


# --- (e) fed-batch experiments with local initial values ---------------------
def test_two_feedspike_experiments_with_local_initial_substrate():
    km = _feedspike_model()
    conds = {
        'FB1': nsk.Condition(initial={'[S]': 50}, parameters={'threshold': 5},
                             t_end=12),
        'FB2': nsk.Condition(initial={'[S]': 30},
                             parameters={'threshold': 10}, t_end=12)}
    t = np.arange(0.5, 12.01, 0.5)
    gen = km.copy()
    rows = []
    for name, c in conds.items():
        res = gen.run(c.with_parameters(k=0.4).replace(
            times=np.concatenate([[0.0], t]), selections=('[S]', '[P]',
                                                          'n_spk')))
        assert res.final['n_spk'] >= 2          # the feed spikes fire
        for sel, ob in (('[S]', 'S'), ('[P]', 'P')):
            rows += [(name, ob, ti, v) for ti, v in zip(t, res[sel][1:])]
    data = Dataset(pd.DataFrame(rows, columns=['experiment', 'observable',
                                               'time', 'value']))
    obs = [Observable('S', '[S]', sigma=0.5),
           Observable('P', '[P]', noise='proportional', rel_sigma=0.05)]
    exps = [Experiment(n, c, data[n]) for n, c in conds.items()]
    params = [FitParameter('k', (0.05, 2), initial=0.3),
              FitParameter('[S]', (10, 80), scope='local',
                           initial={'FB1': 45, 'FB2': 35})]
    before = _model_state(km)
    prob = CalibrationProblem(km, exps, obs, params)
    assert prob.names == ['k', '[S]@FB1', '[S]@FB2']
    assert prob.event_counters == ('n_spk',)
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert np.allclose(fit.estimates.to_numpy(), [0.4, 50.0, 30.0],
                       rtol=1e-4)
    assert fit.experiments == [None, 'FB1', 'FB2']
    assert _model_state(km) == before
    # No override leaks: FB2 alone gives the FB2 block of a joint problem
    # whose FB1 also sets parameters and an initial value FB2 does not.
    alone = CalibrationProblem(km, [exps[1]], obs, params)
    fb1x = Experiment('FB1', conds['FB1'].with_parameters(
        target=80, feed_conc=300).with_initial({'[P]': 3.0}), data['FB1'])
    joint = CalibrationProblem(km, [fb1x, exps[1]], obs, params)
    theta = {'k': 0.37, '[S]@FB1': 44.0, '[S]@FB2': 33.0}
    r_joint, r_prob = joint.residuals(theta), prob.residuals(theta)
    n1 = len(fb1x.data)
    assert np.array_equal(r_joint[n1:], alone.residuals(
        {'k': 0.37, '[S]@FB2': 33.0}))
    assert np.array_equal(r_joint[n1:], r_prob[n1:])
    assert not np.allclose(r_joint[:n1], r_prob[:n1])  # FB1's extras apply
    assert list(prob.data['experiment'][n1:]) == ['FB2'] * len(exps[1].data)
    # apply(): writes k, returns the local initial values
    target = _feedspike_model()
    rest = fit.apply(target)
    assert math.isclose(target._te['k'], fit.estimates['k'])
    assert set(rest) == {'[S]@FB1', '[S]@FB2'}
    assert target._te['init([S])'] == 50.0      # locals are not applied
    assert _model_state(km) == before


# --- (f) integration failures ------------------------------------------------
def test_integration_failures_are_penalized_not_raised(blow_fit,
                                                        blow_problem):
    fit = blow_fit
    assert fit.n_failed_evals > 0
    assert 'failed' in set(fit.starts['status'])
    assert math.isclose(fit.estimates['k'], 0.5, rel_tol=1e-6)
    assert 'failed_evaluations' in set(fit.flags['flag'])
    failed = fit.starts[fit.starts['status'] == 'failed']
    assert failed['cost'].isna().all()
    assert failed['message'].str.contains('start point failed').all()
    # residuals() never raises either: a constant penalty, counted
    n = blow_problem.n_failed_evals
    r = blow_problem.residuals({'k': 2.5})
    assert np.all(r == cal.PENALTY_RESIDUAL)
    assert blow_problem.n_failed_evals == n + 1


def test_every_start_failing_raises_fitting_error():
    t = np.linspace(0.1, 0.9, 9)
    data = pd.DataFrame({'observable': 'x', 'time': t, 'value': 1 + t})
    exp = Experiment('E', nsk.Condition(t_end=0.9), data)
    prob = CalibrationProblem(_model(BLOW_TOY), [exp],
                              [Observable('x', 'x', sigma=0.05)],
                              [FitParameter('k', (2, 5))])
    with pytest.raises(FittingError, match='All 2 start'):
        fit_parameters(prob, n_starts=2, seed=0)


# --- (g) validation before any simulation ------------------------------------
@pytest.fixture
def run_counter(monkeypatch):
    calls = []
    run, simulate = nsk.KineticModel.run, nsk.KineticModel.simulate

    def counted_run(self, *args, **kwargs):
        calls.append('run')
        return run(self, *args, **kwargs)

    def counted_simulate(self, *args, **kwargs):
        calls.append('simulate')
        return simulate(self, *args, **kwargs)

    monkeypatch.setattr(nsk.KineticModel, 'run', counted_run)
    monkeypatch.setattr(nsk.KineticModel, 'simulate', counted_simulate)
    return calls


def _exp(data=None, **cond):
    cond = cond or {'t_end': 6}
    return Experiment('E', nsk.Condition(**cond),
                      data if data is not None else _decay_data())


_OBS = [Observable('S', '[S]', sigma=0.1)]
_K = [FitParameter('k', (0.01, 10))]


@pytest.mark.parametrize('build, match', [
    (lambda km: CalibrationProblem(km, [_exp()], _OBS,
                                   [FitParameter('kk', (0.1, 1))]),
     r"unknown model id\(s\) \['kk'\]"),
    (lambda km: CalibrationProblem(km, [_exp()], _OBS, []),
     'parameters must not be empty'),
    (lambda km: CalibrationProblem(
        km, [_exp(_decay_data(times=np.array([1.0, 2.0])))], _OBS,
        [FitParameter('k', (0.01, 10)), FitParameter('[S]', (1, 20)),
         FitParameter('[P]', (0.01, 1))]),
     'exceed the 2 data points'),
    (lambda km: CalibrationProblem(km, [_exp()], _OBS,
                                   [FitParameter('k', (0.1, 1)),
                                    FitParameter('k', (0.2, 2))]),
     'Duplicate target'),
    (lambda km: CalibrationProblem(km, [_exp()],
                                   [Observable('S', '[Q]', sigma=0.1)], _K),
     r"observable 'S' selection: unknown model id\(s\) \['\[Q\]'\]"),
    (lambda km: CalibrationProblem(km, [_exp()],
                                   [Observable('P', '[P]', sigma=0.1)], _K),
     r"data for observable\(s\) \['S'\]"),
    (lambda km: CalibrationProblem(km, [_exp()],
                                   [Observable('S', '[S]', scale='a',
                                               sigma=0.1)], _K),
     r"\['a'\] have no FitParameter"),
    (lambda km: CalibrationProblem(km, [_exp()], [Observable('S', '[S]')],
                                   _K),
     'no usable sigma'),
    (lambda km: CalibrationProblem(km, [_exp(), _exp()], _OBS, _K),
     'Duplicate experiment'),
    (lambda km: CalibrationProblem(km, [_exp(t_end=3)], _OBS, _K),
     'outside the simulated span'),
    (lambda km: CalibrationProblem(km, [_exp(t_end=6, reset=False)], _OBS,
                                   _K),
     'reset=False'),
    (lambda km: CalibrationProblem(km, [_exp()], _OBS, _K, rel_step=0),
     'rel_step'),
    (lambda km: CalibrationProblem(km, [_exp()], _OBS, _K,
                                   tolerances=(-1, 1e-6)),
     'tolerances'),
    (lambda km: CalibrationProblem(km, [_exp()], _OBS, _K,
                                   max_refinements=-1),
     'max_refinements'),
    (lambda km: CalibrationProblem(km, [_exp()], _OBS, _K, loss='soft_l1',
                                   absolute_sigma=False),
     r"loss='soft_l1' needs absolute sigmas"),
    # a scale name that is also a model id
    (lambda km: CalibrationProblem(km, [_exp()],
                                   [Observable('S', '[S]', scale='k',
                                               sigma=0.1)],
                                   [FitParameter('k', (0.1, 10))]),
     'also model ids'),
    # the condition sets the amount 'S', the fit the concentration '[S]'
    (lambda km: CalibrationProblem(km, [_exp(t_end=6, initial={'S': 20})],
                                   _OBS, [FitParameter('[S]', (1, 50))]),
     'cannot be merged'),
    # the scale of an observable without data is in no residual (it used
    # to be reported as an unknown model id)
    (lambda km: CalibrationProblem(
        km, [_exp()], _OBS + [Observable('Q', '[P]', scale='a', sigma=0.1)],
        _K + [FitParameter('a', (0.1, 10))]),
     r"\['a'\] are the scale/offset of observable\(s\) \['Q'\] that have "
     'no data'),
    # the same for an experiment's own observable
    (lambda km: CalibrationProblem(
        km, [Experiment('E', nsk.Condition(t_end=6), _decay_data(),
                        observables=['S', Observable('Q', '[P]', offset='b',
                                                     sigma=0.1)])],
        _OBS, _K + [FitParameter('b', (-1, 1), log=False)]),
     r"\['b'\] are the scale/offset of observable\(s\) \['Q'\] that have "
     'no data'),
])
def test_invalid_problems_raise_before_any_simulation(run_counter, build,
                                                      match):
    km = _model(DECAY_TOY)
    with pytest.raises(ValueError, match=match):
        build(km)
    assert run_counter == []


@pytest.mark.parametrize('kwargs, match', [
    ({'target': 'k', 'bounds': (0.0, 1.0)}, 'lo > 0'),
    ({'target': 'k', 'bounds': (-1.0, 1.0)}, r'lo > 0.*Pass log=False'),
    ({'target': 'k', 'bounds': (2.0, 1.0)}, 'lo < hi'),
    ({'target': 'k', 'bounds': (0.1, 1.0), 'initial': 2.0}, 'outside'),
    ({'target': 'k', 'bounds': (0.1, 1.0), 'scope': 'both'}, 'scope'),
    ({'target': '[S]', 'bounds': (1, 5), 'scope': 'local',
      'initial': {'E': 9.0}}, 'outside'),
])
def test_invalid_fit_parameters(kwargs, match):
    with pytest.raises(ValueError, match=match):
        FitParameter(**kwargs)
    FitParameter('k', (-1.0, 1.0), log=False)      # linear may cross 0
    with pytest.raises(TypeError, match='scope=.local.'):
        FitParameter('k', (0.1, 1), initial={'E': 0.5})


@pytest.mark.parametrize('kwargs, match', [
    ({'noise': 'poisson'}, 'noise must be one of'),
    ({'noise': 'proportional'}, 'needs rel_sigma'),
    ({'noise': 'combined', 'rel_sigma': 0.1}, 'needs sigma'),
    ({'noise': 'additive', 'rel_sigma': 0.1}, 'rel_sigma is not used'),
    ({'noise': 'proportional', 'rel_sigma': 0.1, 'sigma': 1.0},
     'sigma is not used'),
    ({'sigma': -1.0}, 'must be > 0'),
    ({'scale': 'a', 'offset': 'a', 'sigma': 1.0}, 'same fit parameter'),
    ({'scale': 0.0, 'sigma': 1.0}, 'nonzero'),
])
def test_invalid_observables(kwargs, match):
    with pytest.raises(ValueError, match=match):
        Observable('S', '[S]', **kwargs)
    with pytest.raises(TypeError, match='callable'):
        Observable('S', '[S]', sigma=1.0, transform='log')


def test_invalid_fit_arguments(run_counter):
    prob = _decay_problem()
    with pytest.raises(ValueError, match='method'):
        fit_parameters(prob, method='mcmc')
    with pytest.raises(ValueError, match='sampler'):
        fit_parameters(prob, sampler='grid')
    with pytest.raises(ValueError, match='n_starts'):
        fit_parameters(prob, n_starts=0)
    with pytest.raises(ValueError, match="method='multistart'"):
        fit_parameters(prob, method='hybrid', workers=2)
    with pytest.raises(ValueError, match="only used by method='hybrid'"):
        fit_parameters(prob, global_options={'de_maxiter': 3})
    with pytest.raises(ValueError, match='cannot set'):
        fit_parameters(prob, method='hybrid',
                       global_options={'parallel': 'threads'})
    with pytest.raises(ValueError, match='cannot set'):
        fit_parameters(prob, method='hybrid', global_options={'bounds': []})
    # unknown keys (typos, scipy.optimize.minimize arguments) and counts
    # below 1 fail before the global stage runs
    for options in ({'de_max_iter': 100}, {'tol': 1e-8}, {'nope': 1}):
        with pytest.raises(ValueError,
                           match=r'Unknown global_options key.*allowed'):
            fit_parameters(prob, method='hybrid', global_options=options)
    with pytest.raises(ValueError, match=r"global_options\['de_popsize'\]"):
        fit_parameters(prob, method='hybrid', global_options={'de_popsize': -3})
    with pytest.raises(TypeError, match='local_maxiter'):
        fit_parameters(prob, method='hybrid',
                       global_options={'local_maxiter': 2.5})
    assert set(cal._GLOBAL_OPTIONS) >= {'de_maxiter', 'de_popsize',
                                        'do_basinhop', 'local_method'}
    assert not set(cal._GLOBAL_OPTIONS) & set(cal._RESERVED_OPTIONS)
    # initial: every free parameter by name, the right shape, not hybrid
    with pytest.raises(ValueError, match=r"initial must give every free "
                       r"parameter \['k'\]; unknown \['kk'\]"):
        fit_parameters(prob, initial={'kk': 0.4})
    with pytest.raises(ValueError, match=r'initial must have shape \(1,\)'):
        fit_parameters(prob, initial=[0.1, 0.2])
    with pytest.raises(ValueError, match="method='hybrid' starts from"):
        fit_parameters(prob, method='hybrid', initial={'k': 0.4})
    with pytest.raises(TypeError, match='CalibrationProblem'):
        fit_parameters(object())
    with pytest.raises(ValueError, match='every free parameter'):
        prob.residuals({'kk': 1.0})
    with pytest.raises(ValueError, match='> 0'):
        prob.residuals({'k': -1.0})
    assert run_counter == []


def test_global_target_set_by_a_condition_warns():
    km = _model(DECAY_TOY)
    exps = [Experiment('E1', nsk.Condition(t_end=6, parameters={'k': 0.3}),
                       _decay_data()),
            Experiment('E2', nsk.Condition(t_end=6, parameters={'k': 0.7}),
                       _decay_data())]
    with pytest.warns(UserWarning, match=r"replaces the value set by the "
                      r"condition\(s\) of experiment\(s\) \['E1', 'E2'\]"):
        CalibrationProblem(km, exps, _OBS, _K)
    exps = [Experiment('E1', nsk.Condition(t_end=6, initial={'[S]': 12}),
                       _decay_data()),
            Experiment('E2', nsk.Condition(t_end=6), _decay_data())]
    with pytest.warns(UserWarning, match=r"'\[S\]'.*replaces the value"):
        CalibrationProblem(km, exps, _OBS, [FitParameter('[S]', (1, 20))])
    # a local target set by the conditions is what scope='local' is for
    with warnings.catch_warnings():
        warnings.simplefilter('error', UserWarning)
        CalibrationProblem(km, exps, _OBS,
                           [FitParameter('[S]', (1, 20), scope='local')])


def test_many_free_parameters_warn():
    t = np.linspace(0.5, 6, 40)
    km = _model(DECAY_TOY)
    exps = [Experiment(f'E{i}', nsk.Condition(t_end=6), _decay_data(times=t))
            for i in range(16)]
    with pytest.warns(UserWarning, match='16 free parameters'):
        prob = CalibrationProblem(km, exps, _OBS,
                                  [FitParameter('[S]', (1, 20),
                                                scope='local')])
    assert prob.n_free == 16 and prob.names[0] == '[S]@E0'


# --- flags -------------------------------------------------------------------
def test_flags_rank_deficiency_correlation_and_bounds():
    km = _model(PRODUCT_TOY)
    exp = Experiment('E', nsk.Condition(t_end=6), _decay_data())
    prob = CalibrationProblem(km, [exp], _OBS,
                              [FitParameter('k1', (0.05, 5), initial=0.3),
                               FitParameter('k2', (0.1, 10), initial=1.2)])
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert math.isclose(fit.estimates['k1'] * fit.estimates['k2'], 0.5,
                        rel_tol=1e-6)
    flags = fit.flags.set_index('flag')
    assert flags.loc['rank_deficient', 'parameters'] == ('k1', 'k2')
    # only k1*k2 is determined: the estimates are anti-correlated
    assert flags.loc['high_correlation', 'value'] < -0.99
    assert np.all(np.isinf(fit.se.to_numpy()))
    assert not fit.table['identifiable'].any()
    text = fit.summary()
    assert text.index('Warnings (') < text.index('Parameters (')
    assert text.index('Parameters (') < text.index('Fit statistics')
    norm = ' '.join(text.split())          # warnings are wrapped
    assert 'rank 1 < 2' in norm and 'time [h]' in norm
    assert 'rank / cond(J D)' in text and _long_lines(text) == []
    assert '  - [rank_deficient] The residual Jacobian' in text
    # an optimum beyond the bounds ends at the bound
    at = fit_parameters(_decay_problem(params=[FitParameter('k', (0.6, 2))]),
                        n_starts=1, seed=0)
    assert at.active_mask.tolist() == [-1]
    row = at.flags.set_index('flag').loc['at_bound']
    assert row['parameters'] == ('k',) and 'lower bound' in row['message']
    # the natural bounds are the FitParameter's exactly (exp(ln 0.35) is
    # not 0.35): estimates and draws never pass them, and a draw beyond a
    # bound is moved onto it exactly
    at = fit_parameters(_decay_problem(params=[FitParameter('k',
                                                            (0.01, 0.35))]),
                        n_starts=1, seed=0)
    assert at.active_mask.tolist() == [1]
    assert at.lower[0] == 0.01 and at.upper[0] == 0.35
    assert math.isclose(at.estimates['k'], 0.35, rel_tol=1e-9)
    assert at.estimates['k'] <= 0.35
    draws = at.sample(200, seed=0)['k']
    assert draws.between(0.01, 0.35).all() and (draws == 0.35).any()


def test_lack_of_fit_is_flagged_for_absolute_sigmas(robust_hybrid_fit):
    data = _decay_data(sigma=0.1, seed=0)
    good = fit_parameters(_decay_problem(data=data, sigma=0.1), n_starts=1,
                          seed=0)
    assert 'lack_of_fit' not in set(good.flags['flag'])
    # sigma 10x too small: chi2 ~ 100 dof, far above the 99.9% quantile
    bad = fit_parameters(_decay_problem(data=data, sigma=0.01), n_starts=1,
                         seed=0)
    assert np.allclose(bad.x, good.x, rtol=1e-6)
    row = bad.flags.iloc[0]
    assert row['flag'] == 'lack_of_fit' and row['parameters'] == ()
    assert row['value'] == bad.reduced_chi2 > 50
    assert bad.chi2 > chi2_distribution.ppf(cal.LACK_OF_FIT_LEVEL, bad.dof)
    assert "method='hybrid'" in row['message']
    text = bad.summary()
    assert _long_lines(text) == []
    norm = ' '.join(text.split())               # warnings are wrapped
    assert norm.index(' '.join(row['message'].split())) < norm.index(
        'Parameters (')
    # absolute_sigma=False: no test
    fit = fit_parameters(_decay_problem(data=data, sigma=0.01,
                                        absolute_sigma=False),
                         n_starts=1, seed=0)
    assert 'lack_of_fit' not in set(fit.flags['flag'])
    # a robust loss: 2 * cost (<= chi2) is tested instead
    rob = fit_parameters(_decay_problem(data=data, sigma=0.01,
                                        loss='soft_l1'), n_starts=1, seed=0)
    row = rob.flags.set_index('flag').loc['lack_of_fit']
    assert row['value'] == pytest.approx(2 * rob.cost / rob.dof, rel=1e-12)
    assert 2 * rob.cost < rob.chi2
    assert '2*cost' in row['message'] and 'outliers' in row['message']
    # with the right sigmas, a robust fit passes
    assert 2 * robust_hybrid_fit.cost < chi2_distribution.ppf(
        cal.LACK_OF_FIT_LEVEL, robust_hybrid_fit.dof)
    assert 'lack_of_fit' not in set(robust_hybrid_fit.flags['flag'])


def test_relative_sigma_without_degrees_of_freedom_has_no_covariance():
    data = _decay_data(sigma=0.1, seed=0, times=np.array([2.0]))
    fit = fit_parameters(_decay_problem(data=data, absolute_sigma=False),
                         n_starts=1, seed=0)
    assert fit.dof == 0 and math.isnan(fit.se['k'])
    row = fit.flags.set_index('flag').loc['no_covariance']
    assert 'dof = 0' in row['message']
    assert 'Warnings: none' not in fit.summary()
    with pytest.raises(ValueError, match='not finite'):
        fit.sample(5, seed=0)
    # absolute sigmas need no degrees of freedom
    fit = fit_parameters(_decay_problem(data=data), n_starts=1, seed=0)
    assert np.isfinite(fit.se['k']) and len(fit.flags) == 0


@pytest.fixture(scope='module')
def two_parameter_fit():
    data = _decay_data(S0=8.0, k=0.3, sigma=0.1, seed=0)
    return fit_parameters(_decay_problem(data=data, params=[
        FitParameter('k', (0.01, 10), initial=0.5),
        FitParameter('[S]', (1, 20), initial=10.0)]), n_starts=1, seed=0)


@pytest.mark.parametrize('c, present, absent', [
    (1e3, (), ('ill_conditioned', 'rank_deficient')),
    (1e5, ('ill_conditioned', 'high_rel_se'), ('rank_deficient',)),
    (2e6, ('rank_deficient',), ('ill_conditioned',)),
])
def test_conditioning_flags(two_parameter_fit, c, present, absent):
    # two unit columns at the angle th = 2 arctan(1/c): the singular values
    # of [a, b] are sqrt(1 +- cos th), so cond = cot(th/2) = c. The columns
    # are then given very different norms, which the equilibration undoes:
    # the flags measure collinearity, not scale
    fit = two_parameter_fit
    u, _ = np.linalg.qr(np.random.default_rng(0).standard_normal(
        (fit.n_data, 2)))
    th = 2 * math.atan(1 / c)
    a, b = u[:, 0], math.cos(th) * u[:, 0] + math.sin(th) * u[:, 1]
    g = dataclasses.replace(fit, jacobian=np.column_stack([100 * a,
                                                           1e-3 * b]))
    assert math.isclose(g.condition_number, c, rel_tol=1e-9)
    flags = set(g.flags['flag'])
    assert set(present) <= flags and not set(absent) & flags
    if 'ill_conditioned' in flags:
        row = g.flags.set_index('flag').loc['ill_conditioned']
        assert math.isclose(row['value'], c ** 2, rel_tol=1e-9)


def test_precision_disparity_is_not_collinearity(two_parameter_fit):
    # orthogonal columns, one determined 1e4 times less precisely: full
    # rank and well conditioned once equilibrated; the weak parameter is
    # reported by its relative SE instead
    fit = two_parameter_fit
    u, _ = np.linalg.qr(np.random.default_rng(0).standard_normal(
        (fit.n_data, 2)))
    g = dataclasses.replace(fit, jacobian=u * [100.0, 100.0 / 1e4])
    assert g.numerical_rank == 2
    assert g.condition_number == pytest.approx(1.0, rel=1e-12)
    assert np.allclose(g.se, [1 / 100.0, 1e4 / 100.0], rtol=1e-12)
    assert not {'rank_deficient', 'ill_conditioned'} & set(g.flags.flag)
    rows = g.flags[g.flags.flag == 'high_rel_se']
    assert list(rows.parameters) == [('[S]',)]
    assert not g.table.loc['[S]', 'identifiable']
    assert g.table.loc['k', 'identifiable']


def test_a_column_within_the_noise_of_the_largest_is_null(two_parameter_fit):
    # the same steps (1e-3 in ln p) and a 2e6 disparity: the [S] column
    # changes the residuals by 5e-7 (< rank_rtol) of the k column over its
    # step, i.e. by no more than integration noise, so it is a null
    # direction and k is as if [S] were not fitted
    fit = two_parameter_fit
    assert np.allclose(fit.fd_flags.loc[fit.names, 'step'], 1e-3)
    u, _ = np.linalg.qr(np.random.default_rng(0).standard_normal(
        (fit.n_data, 2)))
    g = dataclasses.replace(fit, jacobian=u * [100.0, 100.0 / 2e6])
    assert g.numerical_rank == 1 and g.condition_number == math.inf
    assert g.se['k'] == pytest.approx(1 / 100.0, rel=1e-12)
    assert math.isinf(g.se['[S]'])
    assert g.correlation.loc['k', '[S]'] == 0.0
    row = g.flags.set_index('flag').loc['rank_deficient']
    assert row['parameters'] == ('[S]',)
    assert "do not measurably depend on ['[S]']" in row['message']
    assert set(g.flags.flag) == {'rank_deficient'}
    # an observable scale or offset involves no integration: never null
    h = dataclasses.replace(g, kinds=['model', 'observable'])
    assert h.numerical_rank == 2 and np.isfinite(h.se['[S]'])
    # an unreliable column (e.g. across an event jump, whose difference is
    # large) does not set the reference while a reliable one exists
    assert fit.fd_flags['reliable'].all()
    flags = fit.fd_flags.copy()
    flags.loc['k', 'reliable'] = False
    assert dataclasses.replace(g, fd_flags=flags).numerical_rank == 2


@pytest.fixture(scope='module')
def offset_fit():
    data = _decay_data(sigma=0.1, seed=0)
    data['value'] += 0.37
    prob = CalibrationProblem(
        _model(DECAY_TOY), [Experiment('E', nsk.Condition(t_end=6), data)],
        [Observable('S', '[S]', sigma=0.1, offset='b')],
        [FitParameter('k', (0.01, 10), initial=0.2),
         FitParameter('b', (-1, 1), log=False, initial=0.0)])
    return fit_parameters(prob, n_starts=1, seed=0)


def test_rank_test_does_not_depend_on_the_origin_or_the_precision(
        offset_fit):
    # an offset b far from 0 (a detector baseline) or data 1e4 times more
    # precise must not change the rank or the flags: with columns scaled
    # by max(|b|, 1/||J_b||), b = 1e7 was reported rank deficient (se of k
    # inf, se of b 4x too small) and b = 0 with precise data likewise
    fit = dataclasses.replace(offset_fit, lower=np.array([0.01, -1e8]),
                              upper=np.array([10.0, 1e8]))
    reference = None
    for b in (0.0, 0.37, 1e7):
        for factor in (1.0, 1e4):
            jac = fit.jacobian * factor
            g = dataclasses.replace(fit, x=np.array([fit.x[0], b]),
                                    jacobian=jac)
            assert g.numerical_rank == 2
            expected = np.sqrt(np.diag(np.linalg.inv(jac.T @ jac)))
            assert np.allclose(g.se, expected, rtol=1e-9, atol=0)
            assert g.table['identifiable'].all()
            state = (list(g.flags.flag), g.condition_number)
            reference = reference or state
            assert state[0] == reference[0] == ['high_correlation']
            assert state[1] == pytest.approx(reference[1], rel=1e-12)


def _unit_fit(u):
    """Decay toy with S and P measured, concentrations in units of 1/u;
    [P]0 (lower bound 0) is fitted in linear space."""
    km = _model(f"""
    model decay()
      S = {10 * u}; P = {0.5 * u}; k = 0.5;
      J: S => P; k*S;
    end
    """)
    rng = np.random.default_rng(0)
    S = 10 * u * np.exp(-0.5 * T)
    P = 10.5 * u - S
    rows = ([('S', t, v + 0.1 * u * rng.standard_normal())
             for t, v in zip(T, S)]
            + [('P', t, v + 0.1 * u * rng.standard_normal())
               for t, v in zip(T, P)])
    data = pd.DataFrame(rows, columns=['observable', 'time', 'value'])
    prob = CalibrationProblem(
        km, [Experiment('E', nsk.Condition(t_end=6), data)],
        [Observable('S', '[S]', sigma=0.1 * u),
         Observable('P', '[P]', sigma=0.1 * u)],
        [FitParameter('k', (0.01, 10)), FitParameter('[S]', (u, 50 * u)),
         FitParameter('[P]', (0, 5 * u), log=False)])
    return fit_parameters(prob, n_starts=2, seed=0)


def test_rank_test_does_not_depend_on_the_units_of_linear_parameters():
    # dr/dp of a linear parameter scales with 1/unit of p: on the raw
    # Jacobian, units of 1e-7 made k and [S] rank deficient (infinite SE)
    # and units of 1e-5 ill-conditioned
    a, b = _unit_fit(1.0), _unit_fit(1e-7)
    for fit in (a, b):
        assert fit.numerical_rank == 3
        assert not {'rank_deficient', 'ill_conditioned'} & set(fit.flags.flag)
    assert list(a.flags.flag) == list(b.flags.flag)
    assert np.allclose(a.rel_se, b.rel_se, rtol=1e-3)
    assert np.allclose(a.ci.div(a.estimates, axis=0),
                       b.ci.div(b.estimates, axis=0), rtol=1e-3)
    assert a.condition_number == pytest.approx(b.condition_number, rel=1e-3)
    assert np.allclose(a.correlation, b.correlation, atol=1e-3)
    # closed form: S = S0 e^(-k t), P = P0 + S0 (1 - e^(-k t)), Gauss-Newton
    # SEs of (ln k, ln S0, P0) from the analytic Jacobian at the optimum
    res = a.residuals
    t, y = res['time'].to_numpy(), res['data'].to_numpy()
    is_s = (res['observable'] == 'S').to_numpy()

    def model(x):
        k, s0 = np.exp(x[:2])
        decay = np.exp(-k * t)
        return np.where(is_s, s0 * decay, x[2] + s0 * (1 - decay)), k, s0

    x = least_squares(lambda x: (model(x)[0] - y) / 0.1,
                      [np.log(0.5), np.log(10), 0.5], xtol=1e-15,
                      ftol=1e-15, gtol=1e-15).x
    _, k, s0 = model(x)
    decay = np.exp(-k * t)
    jac = np.column_stack([
        np.where(is_s, -k * t * s0 * decay, k * t * s0 * decay),
        np.where(is_s, s0 * decay, s0 * (1 - decay)),
        np.where(is_s, 0.0, 1.0)]) / 0.1
    se_p0 = math.sqrt(np.linalg.inv(jac.T @ jac)[2, 2])
    assert a.rel_se['[P]'] == pytest.approx(se_p0 / x[2], rel=1e-3)


def test_relative_se_is_not_judged_for_a_zero_capable_linear_parameter():
    # an offset b (truth 0, bounds (-1, 1)) estimated near 0: its rel_se
    # se/|b| depends on the origin, not on how well the data determine it
    data = _decay_data(sigma=0.2, seed=0)
    prob = CalibrationProblem(
        _model(DECAY_TOY), [Experiment('E', nsk.Condition(t_end=6), data)],
        [Observable('S', '[S]', offset='b', sigma=0.2)],
        [FitParameter('k', (0.01, 10)),
         FitParameter('b', (-1, 1), log=False)])
    fit = fit_parameters(prob, n_starts=3, seed=0)
    assert fit.rel_se['b'] > cal.REL_SE_MAX       # still reported
    assert fit.ci.loc['b', 'lower'] < 0 < fit.ci.loc['b', 'upper']
    assert not {'high_rel_se', 'wide_interval'} & set(fit.flags.flag)
    # judged against its bounds instead: the 95% interval fits in them
    assert 2 * Z * fit.se['b'] <= 2
    assert fit.table.loc['b', 'identifiable']
    assert fit.table.loc['k', 'identifiable']


@pytest.fixture(scope='module')
def k_only_noisy_fit():
    return fit_parameters(_decay_problem(data=_decay_data(sigma=0.1, seed=0)),
                          n_starts=1, seed=0)


@pytest.mark.parametrize('case', ['unobserved', 'unobserved_at_bound',
                                  'weak'])
def test_an_undetermined_zero_capable_linear_parameter_is_flagged(
        case, k_only_noisy_fit):
    # a log=False parameter with a bound at 0 has no relative scale: one
    # the data do not determine within its bounds must still be reported
    # (it was identifiable with no warning, while its profile does not
    # close). [P] does not reach the observed S at all: its column is
    # integration noise, a null direction (it was a weak one, with se of k
    # 2.1x too large); [X] changes S by a small but real kx * X (near its
    # estimate X = 0 the step is 5e-6, over which its column changes the
    # residuals by 2.7e-5 of k's, well above rank_rtol)
    if case == 'weak':
        km = _model('''
        model weak()
          species X; S = 10; P = 0; k = 0.5; X = 0; kx = 1e-2;
          J: S => P; k*S + kx*X;
        end
        ''')
        params = [FitParameter('k', (0.01, 10)),
                  FitParameter('[X]', (0, 5), log=False)]
        prob = _decay_problem(km=km, params=params)
    else:
        initial = 5.0 if case == 'unobserved' else None
        params = [FitParameter('k', (0.01, 10)),
                  FitParameter('[P]', (0, 10), log=False, initial=initial)]
        prob = _decay_problem(data=_decay_data(sigma=0.1, seed=0),
                              params=params)
    fit = fit_parameters(prob, n_starts=1, seed=0)
    name = params[1].target
    width = params[1].bounds[1] - params[1].bounds[0]
    assert fit.estimates['k'] == pytest.approx(0.5, rel=0.05)
    assert not fit.table.loc[name, 'identifiable']
    assert fit.table.loc['k', 'identifiable']
    assert 'Warnings: none' not in fit.summary()
    assert 'high_correlation' not in set(fit.flags.flag)
    if case == 'weak':
        assert np.isfinite(fit.se[name])
        row = fit.flags.set_index('flag').loc['wide_interval']
        assert row['parameters'] == (name,)
        assert row['value'] == pytest.approx(2 * Z * fit.se[name] / width,
                                             rel=1e-12)
        assert row['value'] > 1
        assert '[wide_interval]' in fit.summary()
        return
    assert math.isinf(fit.se[name])
    assert fit.se['k'] == pytest.approx(k_only_noisy_fit.se['k'], rel=1e-3)
    row = fit.flags.set_index('flag').loc['rank_deficient']
    assert row['parameters'] == (name,)
    assert 'wide_interval' not in set(fit.flags.flag)
    assert '[rank_deficient]' in fit.summary()
    if case == 'unobserved_at_bound':
        assert 'at_bound' in set(fit.flags.flag)


UNOBSERVED_TOY = """
model unobserved()
  S = 10; P = 0; k = 0.5; Q = 5; kq = 0.3;
  J: S => P; k*S;
  J2: Q => ; kq*Q;
end
"""


@pytest.mark.parametrize('target', ['kq', '[P]'])
def test_a_parameter_the_outputs_do_not_depend_on_is_a_null_direction(
        target, k_only_noisy_fit):
    # kq (an unobserved side branch) and [P] (not fed back on S), fitted
    # in log space: their columns are integration noise, which follows the
    # dynamics. Equilibrated to unit norm, it looked like a genuine
    # direction nearly collinear with k: se of k 6x too large (with a
    # tolerance-dependent factor), a false high_correlation flag and too
    # wide samples. It is a null direction: k is as if it were not fitted
    # (as in identifiability and in the profile of k), whatever the
    # tolerances
    ref = k_only_noisy_fit
    for tolerances in (None, (1e-14, 1e-10)):
        prob = _decay_problem(
            km=_model(UNOBSERVED_TOY), data=_decay_data(sigma=0.1, seed=0),
            params=[FitParameter('k', (0.01, 10), initial=0.2),
                    FitParameter(target, (0.01, 10))],
            tolerances=tolerances)
        fit = fit_parameters(prob, n_starts=1, seed=0)
        assert fit.estimates['k'] == pytest.approx(ref.estimates['k'],
                                                   rel=1e-4)
        assert fit.se['k'] == pytest.approx(ref.se['k'], rel=1e-3)
        assert np.allclose(fit.ci.loc['k'], ref.ci.loc['k'], rtol=1e-3)
        assert math.isinf(fit.se[target]) and fit.numerical_rank == 1
        assert fit.correlation.loc['k', target] == 0.0
        flags = fit.flags.set_index('flag')
        assert flags.loc['rank_deficient', 'parameters'] == (target,)
        assert not {'high_correlation', 'high_rel_se'} & set(flags.index)
        assert fit.table.loc['k', 'identifiable']
        assert not fit.table.loc[target, 'identifiable']
        with pytest.warns(UserWarning, match='infinite standard errors'):
            draws = fit.sample(400, seed=0)
        assert np.std(np.log(draws['k'])) == pytest.approx(ref.se['k'],
                                                           rel=0.1)


def test_correlation_on_a_null_direction_is_the_vanishing_ridge_limit(
        two_parameter_fit):
    fit = two_parameter_fit
    a, b, c = np.random.default_rng(1).standard_normal((3, fit.n_data))
    names = ['p1', 'p2', 'p3']

    def with_jacobian(jac):
        return dataclasses.replace(
            fit, names=names, targets=names, experiments=[None] * 3,
            kinds=['model'] * 3, log=np.ones(3, dtype=bool),
            lower=np.full(3, 0.01), upper=np.full(3, 100.0), x=np.zeros(3),
            active_mask=np.zeros(3, dtype=int), jacobian=jac)

    # p1 and p2 enter only as p1*p2 (equal columns in log space)
    jac = np.column_stack([a, a, b])
    g = with_jacobian(jac)
    assert g.numerical_rank == 2
    assert np.isinf(g.se[['p1', 'p2']]).all() and np.isfinite(g.se['p3'])
    corr = g.correlation
    assert corr.loc['p1', 'p2'] == pytest.approx(-1.0, abs=1e-9)
    assert corr.loc['p1', 'p3'] == 0.0 and corr.loc['p2', 'p3'] == 0.0
    assert np.array_equal(np.diag(corr), np.ones(3))
    assert np.allclose(corr, corr.T, rtol=0, atol=0)
    flags = g.flags.set_index('flag')
    assert flags.loc['high_correlation', 'parameters'] == ('p1', 'p2')
    # the limit of a vanishing ridge, (J^T J + eps I)^-1
    w, v = np.linalg.eigh(jac.T @ jac)
    ridge = (v / (w + 1e-8)) @ v.T
    d = np.sqrt(np.diag(ridge))
    assert np.allclose(corr, ridge / np.outer(d, d), rtol=0, atol=1e-4)
    # continuous with a full-rank, nearly unidentifiable Jacobian
    near = with_jacobian(np.column_stack([a, a + 1e-3 * c, b]))
    assert near.numerical_rank == 3 and np.isfinite(near.se).all()
    assert near.correlation.loc['p1', 'p2'] < -0.99


def test_no_covariance_and_not_converged_flags(two_parameter_fit):
    fit = two_parameter_fit
    broken = dataclasses.replace(fit,
                                 jacobian=np.full_like(fit.jacobian, np.nan))
    assert 'no_covariance' in set(broken.flags['flag'])
    assert 'no_covariance' not in set(fit.flags['flag'])
    short = fit_parameters(fit.problem, n_starts=1, seed=0, max_nfev=1)
    assert short.flags.iloc[0]['flag'] == 'not_converged'


def test_fd_flags_propagate_into_flags(decay_fit):
    fd = decay_fit.fd_flags.copy()
    fd.loc['k', 'event_count_changed'] = True
    fd.loc['k', 'reliable'] = False
    fit = dataclasses.replace(decay_fit, fd_flags=fd)
    assert list(fit.flags['flag']) == ['event_count_changed']
    assert fit.table.loc['k', 'event_count_changed']
    fd.loc['k', 'event_count_changed'] = False
    fd.loc['k', 'error'] = 'every run failed'
    fit = dataclasses.replace(decay_fit, fd_flags=fd)
    row = fit.flags.iloc[0]
    assert row['flag'] == 'unreliable_derivative'
    assert 'every run failed' in row['message']


# --- methods and workers -----------------------------------------------------
def test_event_count_change_from_a_real_fit_is_flagged():
    # a sample 1e-5 h after the first feed spike (at ln 3 / k for [S] 30,
    # threshold 10): perturbing k moves the spike across it, which the FD
    # refinements cannot resolve
    km = _feedspike_model()
    cond = nsk.Condition(initial={'[S]': 30}, parameters={'threshold': 10},
                         t_end=12)
    t = np.sort(np.append(np.arange(0.5, 12.01, 0.5),
                          math.log(3) / 0.4 + 1e-5))
    res = km.copy().run(cond.with_parameters(k=0.4).replace(
        times=np.concatenate([[0.0], t]), selections=('[S]', '[P]')))
    data = pd.DataFrame(
        [(ob, ti, v) for sel, ob in (('[S]', 'S'), ('[P]', 'P'))
         for ti, v in zip(t, res[sel][1:])],
        columns=['observable', 'time', 'value'])
    obs = [Observable('S', '[S]', sigma=0.5),
           Observable('P', '[P]', noise='proportional', rel_sigma=0.05)]
    prob = CalibrationProblem(km, [Experiment('FB', cond, data)], obs,
                              [FitParameter('k', (0.05, 2), initial=0.41)])
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert math.isclose(fit.estimates['k'], 0.4, rel_tol=1e-4)
    assert fit.fd_flags.loc['k', 'event_count_changed']
    assert not fit.fd_flags.loc['k', 'reliable']
    assert 'event_count_changed' in set(fit.flags['flag'])
    assert fit.table.loc['k', 'event_count_changed']
    text = fit.summary()
    assert 0 < text.index('a perturbation changes how many times') < \
        text.index('Parameters (')


def _long_lines(text):
    """Summary lines wider than SUMMARY_WIDTH, the parameter table aside."""
    lines = text.splitlines()
    start = next(i for i, ln in enumerate(lines)
                 if ln.startswith('Parameters ('))
    end = lines.index('Fit statistics:')
    return [ln for i, ln in enumerate(lines)
            if not start < i < end and len(ln) > cal.SUMMARY_WIDTH]


def test_summary_wraps_lines_and_counts_starts_at_the_best_cost(
        decay_fit, two_parameter_fit):
    starts = decay_fit.starts.copy()
    best = int(starts.index[starts['best']][0])
    others = [i for i in starts.index if i != best]
    starts.loc[others[0], 'cost'] = decay_fit.cost + 0.4  # delta chi2 0.8
    starts.loc[others[1], 'cost'] = decay_fit.cost + 0.6  # delta chi2 1.2
    fit = dataclasses.replace(decay_fit, starts=starts)
    assert re.search(r'\n  starts at the best cost +2 of 3\n', fit.summary())
    # absolute_sigma=False: delta chi2 in units of the reduced chi-square (4)
    rel = dataclasses.replace(fit, absolute_sigma=False,
                              chi2=4.0 * fit.dof)
    assert re.search(r'starts at the best cost +3 of 3', rel.summary())
    # a failed start is not counted
    starts.loc[others[0], 'status'] = 'failed'
    failed = dataclasses.replace(decay_fit, starts=starts)
    assert re.search(r'starts at the best cost +1 of 2', failed.summary())
    # one start: no count
    assert 'starts at the best cost' not in two_parameter_fit.summary()
    # summary() leaves the result's frames as they were
    for frame in (two_parameter_fit.estimates, two_parameter_fit.covariance,
                  two_parameter_fit.correlation, two_parameter_fit.ci,
                  two_parameter_fit.table):
        assert frame.index.name == 'parameter'
    # long warnings and messages are wrapped, with hanging indents
    long = 'x' * 30 + ' word' * 40
    flagged = dataclasses.replace(decay_fit, success=False, message=long)
    text = flagged.summary()
    assert _long_lines(text) == []
    lines = text.splitlines()
    i = next(i for i, ln in enumerate(lines)
             if ln.startswith('  - [not_converged] The optimizer did not '
                              'report'))
    assert lines[i + 1].startswith('    ') and lines[i + 1][4] != ' '
    j = next(i for i, ln in enumerate(lines) if ln.startswith('  message'))
    indent = lines[j].index('x' * 30)
    assert lines[j + 1].startswith(' ' * indent + 'word')
    # the flags keep their messages unwrapped
    assert flagged.warnings[0].endswith(long)
    # the header: the method and loss, then the counts, never split mid-line
    assert two_parameter_fit.summary().splitlines()[:2] == [
        "Calibration fit (multistart, loss='linear')",
        '2 free parameter(s), 12 data points, 1/1 start(s) evaluated']


def test_repr_shows_the_estimates(two_parameter_fit):
    text = repr(two_parameter_fit)
    k, s0 = two_parameter_fit.estimates
    assert f"estimates={{k: {k:.4g}, [S]: {s0:.4g}}}" in text
    assert text.startswith("FitResult(method='multistart', names=['k', '[S]']")


def test_hybrid_method_recovers_k():
    prob = _decay_problem()
    fit = fit_parameters(prob, method='hybrid', seed=0,
                         global_options={'de_maxiter': 5, 'de_popsize': 5})
    assert fit.method == 'hybrid' and len(fit.starts) == 1
    assert math.isclose(fit.estimates['k'], 0.5, rel_tol=1e-6)
    # the count is the whole global stage's, not the L-BFGS-B polish's
    assert re.match(r'Global stage \(\d+ objective evaluations\): DE ',
                    fit.message)
    # nfev includes the global stage, not only the polish
    assert fit.starts['nfev'].iloc[0] > 10
    text = ' '.join(fit.summary().split())
    assert 'global stage + least-squares polish' in text
    assert 'start(s) evaluated' not in text
    assert 'starts at the best cost' not in text
    again = fit_parameters(prob, method='hybrid', seed=0,
                           global_options={'de_maxiter': 5, 'de_popsize': 5})
    assert np.array_equal(again.x, fit.x)


def test_lhs_starts_are_deterministic():
    fit = fit_parameters(_decay_problem(), n_starts=3, sampler='lhs', seed=0)
    assert math.isclose(fit.estimates['k'], 0.5, rel_tol=1e-6)
    again = fit_parameters(_decay_problem(), n_starts=3, sampler='lhs',
                           seed=0)
    pd.testing.assert_frame_equal(again.starts, fit.starts)
    sobol = fit_parameters(_decay_problem(), n_starts=3, seed=0)
    assert not np.allclose(fit.starts['initial:k'].iloc[1:],
                           sobol.starts['initial:k'].iloc[1:])


def test_explicit_evaluation_options_are_used():
    prob = _decay_problem(tolerances=(1e-10, 1e-8), rel_step=1e-4,
                          max_refinements=0, event_counters=())
    assert prob.tolerances == (1e-10, 1e-8) and prob.rel_step == 1e-4
    assert prob.max_refinements == 0 and prob.event_counters == ()
    integrator = prob._private._te.getIntegrator()
    assert integrator.getValue('relative_tolerance') == 1e-8
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert math.isclose(fit.estimates['k'], 0.5, rel_tol=1e-5)
    assert fit.fd_flags.loc['k', 'step'] == 1e-4


def test_two_workers_match_serial_starts(blow_problem, blow_fit):
    fit = fit_parameters(blow_problem, n_starts=4, seed=0, workers=2)
    assert np.allclose(fit.x, blow_fit.x, rtol=1e-10)
    assert list(fit.starts['status']) == list(blow_fit.starts['status'])
    assert fit.n_failed_evals == blow_fit.n_failed_evals
    assert np.allclose(fit.starts['cost'], blow_fit.starts['cost'],
                       rtol=1e-6, equal_nan=True)


def test_spawn_workers_name_an_unpicklable_transform(monkeypatch):
    import multiprocessing

    from nskinetics.engine import execution
    exp = Experiment('E', nsk.Condition(t_end=6), _decay_data())
    params = [FitParameter('k', (0.01, 10), initial=0.2)]
    obs = [Observable('S', '[S]', sigma=0.1, transform=lambda x: 1.0 * x)]
    prob = CalibrationProblem(_model(DECAY_TOY), [exp], obs, params)
    fit = fit_parameters(prob, n_starts=1, seed=0)   # serial: not pickled
    # an experiment's own observable is checked too
    own = Experiment('E', nsk.Condition(t_end=6), _decay_data(),
                     observables=obs)
    shared = [Observable('S', '[S]', sigma=0.1, transform=np.sqrt)]
    own_prob = CalibrationProblem(_model(DECAY_TOY), [own], shared, params)
    monkeypatch.setattr(execution, '_mp_context',
                        lambda: multiprocessing.get_context('spawn'))
    match = "Observable 'S': transform .* is not picklable.*module-level"
    for call in (lambda: fit_parameters(prob, n_starts=2, workers=2),
                 lambda: fit_parameters(own_prob, n_starts=2, workers=2),
                 lambda: profile_likelihood(fit, ['k'], workers=2)):
        with pytest.raises(ValueError, match=match):
            call()
    # a module-level function passes the check
    cal._check_transforms_picklable(CalibrationProblem(
        _model(DECAY_TOY), [exp], shared, params))


# --- (h) CSV round trip and apply --------------------------------------------
def _assert_same_fit(a, b):
    for name in ('names', 'targets', 'experiments', 'kinds', 'cost', 'chi2',
                 'n_data', 'n_failed_evals', 'success', 'message', 'method',
                 'loss', 'f_scale', 'absolute_sigma', 'seed', 'units',
                 'rank_rtol', 'dof', 'reduced_chi2', 'aic', 'bic',
                 'numerical_rank', 'condition_number'):
        assert getattr(a, name) == getattr(b, name), name
    assert (a.aicc == b.aicc) or (math.isnan(a.aicc) and math.isnan(b.aicc))
    for name in ('log', 'lower', 'upper', 'x', 'jacobian', 'active_mask'):
        assert np.array_equal(getattr(a, name), getattr(b, name)), name
    for name in ('residuals', 'starts', 'fd_flags', 'flags', 'covariance',
                 'correlation', 'ci', 'table'):
        pd.testing.assert_frame_equal(getattr(a, name), getattr(b, name))
    for name in ('estimates', 'se', 'se_natural', 'rel_se'):
        pd.testing.assert_series_equal(getattr(a, name), getattr(b, name))
    assert a.summary() == b.summary()


@pytest.fixture(scope='module')
def local_scale_fit():
    """Local linear-space observable scales, absolute_sigma=False."""
    scaled = _decay_data()
    scaled['value'] *= 3.0
    exps = [Experiment('E1', nsk.Condition(t_end=6), scaled),
            Experiment('E2', nsk.Condition(t_end=6),
                       _decay_data(sigma=0.1, seed=1))]
    prob = CalibrationProblem(
        _model(DECAY_TOY), exps,
        [Observable('S', '[S]', scale='a', sigma=0.1)],
        [FitParameter('k', (0.01, 10)),
         FitParameter('a', (0.1, 10), scope='local', log=False)],
        absolute_sigma=False)
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert fit.experiments == [None, 'E1', 'E2']
    assert fit.kinds == ['model', 'observable', 'observable']
    assert fit.log.tolist() == [True, False, False]
    return fit


@pytest.fixture(scope='module')
def robust_hybrid_fit():
    data = _decay_data(sigma=0.1, seed=0)
    fit = fit_parameters(_decay_problem(data=data, loss='soft_l1',
                                        f_scale=2.0),
                         method='hybrid', seed=0)
    assert (fit.method, fit.loss, fit.f_scale) == ('hybrid', 'soft_l1', 2.0)
    return fit


@pytest.mark.parametrize('which', ['decay_fit', 'blow_fit',
                                   'two_parameter_fit', 'local_scale_fit',
                                   'robust_hybrid_fit'])
def test_csv_round_trip_is_exact(which, request, tmp_path):
    fit = request.getfixturevalue(which)
    path = tmp_path / 'fit.csv'
    fit.to_csv(path)
    back = FitResult.from_csv(path)
    assert back.problem is None
    _assert_same_fit(fit, back)
    back.to_csv(tmp_path / 'again.csv')
    assert path.read_text() == (tmp_path / 'again.csv').read_text()


def test_from_csv_rejects_foreign_and_truncated_files(decay_fit, tmp_path):
    path = tmp_path / 'fit.csv'
    decay_fit.to_csv(path)
    lines = path.read_text().splitlines()
    truncated = tmp_path / 'cut.csv'
    truncated.write_text('\n'.join(lines[:-3]) + '\n')
    with pytest.raises(ValueError, match='incomplete'):
        FitResult.from_csv(truncated)
    other = tmp_path / 'other.csv'
    other.write_text('a,b\n1,2\n')
    with pytest.raises(ValueError, match='not written by FitResult'):
        FitResult.from_csv(other)


def test_apply_writes_global_estimates_through_the_initial_route():
    data = _decay_data(S0=8.0, k=0.3)
    prob = _decay_problem(data=data, params=[
        FitParameter('k', (0.01, 10), initial=0.5),
        FitParameter('[S]', (1, 20), initial=10.0)])
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert np.allclose(fit.estimates.to_numpy(), [0.3, 8.0], rtol=1e-6)
    km = _model(DECAY_TOY)
    assert fit.apply(km) == {}
    assert math.isclose(km._te['k'], fit.estimates['k'])
    assert math.isclose(km._te['init([S])'], fit.estimates['[S]'])
    km.reset()                              # persists through baseline reset
    res = km.run(nsk.Condition(times=np.concatenate([[0.0], T])))
    assert math.isclose(res['[S]'][0], fit.estimates['[S]'])
    assert np.allclose(res['[S]'][1:], data['value'], rtol=1e-5)
    # a model without the target: ValueError before any write
    other = _model(BLOW_TOY)
    before = _model_state(other)
    with pytest.raises(ValueError, match='unknown model id'):
        fit.apply(other)
    assert _model_state(other) == before


# --- (h) Laplace samples and prediction bands --------------------------------
@pytest.fixture(scope='module')
def noisy_fit():
    """k of S = 10 e^(-0.5 t) with seeded noise (sigma 0.3)."""
    data = _decay_data(sigma=0.3, seed=0)
    return fit_parameters(_decay_problem(data=data, sigma=0.3), n_starts=2,
                          seed=0)


def test_sample_is_reproducible_log_normal_and_within_bounds(noisy_fit):
    fit = noisy_fit
    draws = fit.sample(2000, seed=1)
    assert draws.equals(fit.sample(2000, seed=1))
    assert not draws.equals(fit.sample(2000, seed=2))
    assert list(draws.columns) == ['k'] and draws.index.name == 'sample'
    assert len(draws) == 2000
    # normal in ln k around the estimate, with the Wald SE
    z = np.log(draws['k'].to_numpy())
    assert abs(z.mean() - fit.x[0]) < 4 * fit.se['k'] / math.sqrt(2000)
    assert math.isclose(z.std(), fit.se['k'], rel_tol=0.05)
    # draws beyond a bound are moved onto it (in natural units, exactly)
    k = fit.estimates['k']
    narrow = dataclasses.replace(fit, lower=np.array([0.98 * k]),
                                 upper=np.array([1.02 * k]))
    clipped = narrow.sample(500, seed=1)['k']
    assert clipped.between(narrow.lower[0], narrow.upper[0]).all()
    assert (clipped == narrow.lower[0]).sum() > 50
    assert (clipped == narrow.upper[0]).sum() > 50


def test_sample_draws_student_t_when_the_noise_level_is_estimated():
    # 4 points and 1 parameter: dof 3, so t(0.975, 3) = 3.18 (not 1.96)
    data = _decay_data(sigma=0.3, seed=0, times=np.array([0.5, 1.5, 3, 5]))
    fit = fit_parameters(_decay_problem(data=data, sigma=0.3,
                                        absolute_sigma=False),
                         n_starts=1, seed=0)
    assert fit.dof == 3
    draws = fit.sample(20000, seed=3)
    assert draws.equals(fit.sample(20000, seed=3))
    z = np.log(draws['k'].to_numpy())
    lo, hi = np.log(fit.ci.loc['k'].to_numpy())
    assert lo > np.log(fit.lower[0]) + 50 * fit.se['k']     # no clipping
    # the draws' 95% range is the Student-t Wald interval
    width = hi - lo
    assert abs(np.percentile(z, 2.5) - lo) < 0.03 * width
    assert abs(np.percentile(z, 97.5) - hi) < 0.03 * width
    assert np.percentile(np.abs(z - fit.x[0]) / fit.se['k'], 95) > 3.0


def test_sample_needs_a_finite_covariance_and_warns_on_null_directions(
        noisy_fit):
    broken = dataclasses.replace(
        noisy_fit, jacobian=np.full_like(noisy_fit.jacobian, np.nan))
    with pytest.raises(ValueError, match='not finite'):
        broken.sample(10, seed=0)
    for n, error in ((0, ValueError), (1.5, TypeError)):
        with pytest.raises(error, match='n must be'):
            noisy_fit.sample(n)
    exp = Experiment('E', nsk.Condition(t_end=6), _decay_data())
    fit = fit_parameters(CalibrationProblem(
        _model(PRODUCT_TOY), [exp], _OBS,
        [FitParameter('k1', (0.05, 5), initial=0.3),
         FitParameter('k2', (0.1, 10), initial=1.2)]), n_starts=1, seed=0)
    with pytest.warns(UserWarning, match='infinite standard errors'):
        draws = fit.sample(20, seed=0)
    assert list(draws.columns) == ['k1', 'k2']


def test_predict_bands_are_ordered_reproducible_and_hold_the_nominal(
        noisy_fit):
    fit = noisy_fit
    cond = nsk.Condition(t_end=8, n_steps=40)
    bands = fit.predict(cond, n_samples=100, seed=1)
    assert list(bands.columns) == ['time', 'output', 'nominal', 'q5', 'q50',
                                   'q95']
    # the model's default outputs (floating species), time-ordered blocks
    assert list(dict.fromkeys(bands['output'])) == ['[S]', '[P]']
    assert bands.attrs == {'n_samples': 100, 'n_failed': 0,
                           'quantiles': (5.0, 50.0, 95.0),
                           'parameters': ['k'], 'seed': 1,
                           'units': UNITS}
    assert (bands['q5'] <= bands['q50']).all()
    assert (bands['q50'] <= bands['q95']).all()
    tol = 1e-9 * (1 + bands['nominal'].abs())
    assert (bands['q5'] - tol <= bands['nominal']).all()
    assert (bands['nominal'] <= bands['q95'] + tol).all()
    assert (bands['q95'] - bands['q5']).max() > 0.1      # a real band
    pd.testing.assert_frame_equal(bands, fit.predict(cond, n_samples=100,
                                                     seed=1))
    # the nominal is the model at the estimate; the bands are percentiles
    # over the draws of sample(seed=...)
    S = bands[bands['output'] == '[S]']
    k_hat = fit.estimates['k']
    assert np.allclose(S['nominal'], 10 * np.exp(-k_hat * S['time']),
                       rtol=1e-6)
    k = fit.sample(100, seed=1)['k'].to_numpy()
    t_end = S['time'].iloc[-1]
    for q in (5, 50, 95):
        assert math.isclose(S[f'q{q}'].iloc[-1], np.percentile(
            10 * np.exp(-k * t_end), q), rel_tol=1e-6)
    custom = fit.predict(cond, n_samples=20, q=(97.5, 2.5), seed=3,
                         outputs='[P]')
    assert list(custom.columns) == ['time', 'output', 'nominal', 'q2.5',
                                    'q97.5']
    assert set(custom['output']) == {'[P]'}


def test_predict_applies_local_estimates_only_for_their_experiment():
    rows = []
    for name, S0 in (('E1', 8.0), ('E2', 12.0)):
        S = S0 * np.exp(-0.5 * T)
        rows += [(name, 'S', t, v) for t, v in zip(T, S)]
        if name == 'E1':
            rows += [(name, 'P', t, v) for t, v in zip(T, S0 - S)]
    data = Dataset(pd.DataFrame(rows, columns=['experiment', 'observable',
                                               'time', 'value']))
    exps = [Experiment(n, nsk.Condition(t_end=6), data[n])
            for n in ('E1', 'E2')]
    fit = fit_parameters(CalibrationProblem(
        _model(DECAY_TOY), exps, _OBS + [Observable('P', '[P]', sigma=0.1)],
        [FitParameter('k', (0.01, 10), initial=0.3),
         FitParameter('[S]', (1, 20), scope='local')]), n_starts=1, seed=0)
    assert np.allclose(fit.estimates, [0.5, 8.0, 12.0], rtol=1e-6)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        own = fit.predict('E2', n_samples=10, seed=0, outputs=['[S]'])
        other = fit.predict(nsk.Condition(t_end=5), n_samples=10, seed=0,
                            outputs=['[S]'])
    assert not [w for w in caught if 'local estimates' in str(w.message)]
    assert other.attrs['parameters'] == ['k']
    assert own.attrs['parameters'] == ['k', '[S]@E2']
    assert math.isclose(own['nominal'].iloc[0], 12.0, rel_tol=1e-9)
    assert own['q5'].iloc[0] < 12.0 < own['q95'].iloc[0]
    # an experiment's own condition: its local estimates are left out,
    # with a warning (both experiments share Condition(t_end=6))
    with pytest.warns(UserWarning, match=r"experiment\(s\) \['E1', 'E2'\]; "
                      r"their local estimates \['\[S\]@E1', '\[S\]@E2'\] "
                      r"are not applied .* Pass an experiment's name"):
        shared = fit.predict(exps[1].condition, n_samples=10, seed=0,
                             outputs=['[S]'])
    assert shared.attrs['parameters'] == ['k']
    assert shared['nominal'].iloc[0] == 10.0   # the model's own S0
    # the plot grid: experiments x observables, a panel without data is off
    fig, axes = fit.plot()
    assert axes.shape == (2, 2)
    assert [ax.get_title() for ax in axes.ravel()[:3]] == [
        'E1: S', 'E1: P', 'E2: S']
    assert not axes[1, 1].axison
    import matplotlib.pyplot as plt
    plt.close('all')


def test_predict_counts_and_warns_failed_samples(blow_fit):
    # A 10x wider Gaussian: draws with k > 1/t_end blow up before t_end.
    wide = dataclasses.replace(blow_fit, jacobian=blow_fit.jacobian / 10)
    cond = nsk.Condition(t_end=1.8, n_steps=20)
    with pytest.warns(UserWarning, match='sampled runs failed'):
        bands = wide.predict(cond, n_samples=40, seed=0)
    k = wide.sample(40, seed=0)['k'].to_numpy()
    n_failed = bands.attrs['n_failed']
    assert n_failed == (k > 1 / 1.8).sum() and 0 < n_failed < 40
    end = bands.iloc[-1]
    ok = k[k <= 1 / 1.8]
    assert math.isclose(end['q50'], np.percentile(1 / (1 - ok * 1.8), 50),
                        rel_tol=1e-5)
    # every run failing raises
    with pytest.warns(UserWarning), pytest.raises(KineticSimulationError,
                                                  match='Every run'):
        wide.predict(nsk.Condition(t_end=5, n_steps=20), n_samples=5, seed=0)


def test_predict_validation_and_from_csv(noisy_fit, run_counter, tmp_path):
    cond = nsk.Condition(t_end=6, n_steps=10)
    cases = [((object(),), {}, TypeError, 'Condition or an experiment'),
             (('nope',), {}, ValueError, 'Unknown experiment'),
             ((cond,), {'n_samples': 0}, ValueError, 'n_samples'),
             ((cond,), {'q': (5, 101)}, ValueError, 'percentiles'),
             ((cond,), {'q': (5, 5)}, ValueError, 'Duplicate'),
             ((cond,), {'outputs': ['[S]', '[S]']}, ValueError, 'Duplicate'),
             ((cond,), {'outputs': []}, ValueError, 'outputs'),
             ((cond,), {'outputs': ['[Q]']}, ValueError, 'Q'),
             ((cond,), {'workers': 0}, ValueError, 'workers'),
             ((cond.replace(reset=False),), {}, ValueError, 'reset=False'),
             ((cond,), {'model': 'km'}, TypeError, 'KineticModel')]
    for args, kwargs, error, match in cases:
        with pytest.raises(error, match=match):
            noisy_fit.predict(*args, **kwargs)
    assert run_counter == []
    path = tmp_path / 'fit.csv'
    noisy_fit.to_csv(path)
    back = FitResult.from_csv(path)
    with pytest.raises(ValueError, match='pass model='):
        back.predict(cond)
    with pytest.raises(ValueError, match='experiment name needs'):
        back.predict('E', model=_model(DECAY_TOY))
    assert run_counter == []
    # with a model, a result read from a file predicts the same bands
    pd.testing.assert_frame_equal(
        back.predict(cond, n_samples=20, seed=4, model=_model(DECAY_TOY)),
        noisy_fit.predict(cond, n_samples=20, seed=4), rtol=1e-4)
    # sampling needs no model or problem
    assert back.sample(5, seed=0).equals(noisy_fit.sample(5, seed=0))


def test_predict_needs_a_fitted_model_parameter(run_counter):
    # only an observable scale was fitted: nothing to apply to a run
    data = _decay_data()
    data['value'] *= 2.0
    exp = Experiment('E', nsk.Condition(t_end=6), data)
    prob = CalibrationProblem(_model(DECAY_TOY), [exp],
                              [Observable('S', '[S]', scale='a', sigma=0.1)],
                              [FitParameter('a', (0.1, 10), initial=1.5)])
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert abs(fit.estimates['a'] / 2 - 1) < 1e-6
    run_counter.clear()
    with pytest.raises(ValueError, match='No fitted model parameter applies'):
        fit.predict(nsk.Condition(t_end=6), n_samples=5, seed=0)
    assert run_counter == []


def test_predict_warns_that_an_observable_name_is_a_model_selection(
        noisy_fit):
    # the observable 'S' records '[S]'; the output 'S' is the amount of S
    cond = nsk.Condition(t_end=6, n_steps=6)
    match = (r"output 'S' is read as the model selection 'S', not as the "
             r"observable 'S', which records '\[S\]'; pass '\[S\]'")
    with pytest.warns(UserWarning, match=match):
        bands = noisy_fit.predict(cond, n_samples=5, seed=0,
                                  outputs=['S', '[P]'])
    assert list(bands['output'].unique()) == ['S', '[P]']
    with pytest.warns(UserWarning, match=match):
        noisy_fit.predict('E', n_samples=5, seed=0, outputs='S')
    # the observable's own selection, or another output: no warning
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        noisy_fit.predict(cond, n_samples=5, seed=0, outputs=['[S]', 'P'])


# --- (k) plot ----------------------------------------------------------------
def test_plot_marks_feed_spikes_on_a_dense_model_curve(tmp_path):
    import matplotlib.pyplot as plt
    km = _feedspike_model()
    cond = nsk.Condition(initial={'[S]': 30}, parameters={'threshold': 10},
                         t_end=12)
    t = np.arange(0.5, 12.01, 0.5)
    res = km.copy().run(cond.with_parameters(k=0.4).replace(
        times=np.concatenate([[0.0], t]), selections=('[S]', '[P]')))
    data = pd.DataFrame({'observable': ['S'] * t.size + ['P'] * t.size,
                         'time': np.concatenate([t, t]),
                         'value': np.concatenate([res['[S]'][1:],
                                                  res['[P]'][1:]])})
    obs = [Observable('S', '[S]', sigma=0.5),
           Observable('P', '[P]', noise='proportional', rel_sigma=0.05)]
    before = _model_state(km)
    fit = fit_parameters(CalibrationProblem(
        km, [Experiment('FB', cond, data)], obs,
        [FitParameter('k', (0.05, 2), initial=0.3)]), n_starts=1, seed=0)
    assert math.isclose(fit.estimates['k'], 0.4, rel_tol=1e-4)
    fig, axes = fit.plot(n_points=241)
    # A7: plotting, predicting and profiling leave the user's model alone
    fit.predict('FB', n_samples=3, seed=0, outputs=['n_spk'])
    profile_likelihood(fit, 'k', n_points=3)
    assert _model_state(km) == before
    assert axes.shape == (1, 2)
    ax = axes[0, 0]
    assert ax.get_title() == 'FB: S'
    lines = {line.get_label(): line for line in ax.get_lines()}
    grid = lines['model'].get_xdata()
    assert grid.size >= 241 and grid[0] == 0.0 and grid[-1] == 12.0
    # the curve goes through the fitted model values at the data times
    at_data = np.isin(grid, t)
    assert np.allclose(lines['model'].get_ydata()[at_data],
                       fit.residuals['model'][:t.size], rtol=1e-6)
    # dotted lines at the first output time after each feed spike: S falls
    # from 30 to the threshold 10, then from the target 50 to 10
    events = np.array([line.get_xdata()[0] for line in ax.get_lines()
                       if line.get_linestyle() == ':'])
    spikes = (math.log(3) + math.log(5) * np.arange(3)) / 0.4
    assert events.size == 3
    assert np.all((events >= spikes) & (events - spikes <= 12 / 240 + 1e-9))
    # data with sigma error bars
    assert len(ax.containers) == 1
    # without the problem: the model at the data times, no event lines
    path = tmp_path / 'fit.csv'
    fit.to_csv(path)
    fig2, axes2 = FitResult.from_csv(path).plot()
    lines2 = axes2[0, 0].get_lines()
    model2 = [line for line in lines2 if line.get_label() == 'model'][0]
    assert np.array_equal(model2.get_xdata(), t)
    assert not [line for line in lines2 if line.get_linestyle() == ':']
    with pytest.raises(ValueError, match='n_points'):
        fit.plot(n_points=1)
    plt.close('all')


def test_plot_falls_back_to_the_data_times_when_the_dense_run_fails():
    import matplotlib.pyplot as plt
    # the data end at 0.9, the condition at 2.5: past the blow-up at 1/k = 2
    t = np.linspace(0.1, 0.9, 9)
    data = pd.DataFrame({'observable': 'x', 'time': t,
                         'value': 1 / (1 - 0.5 * t)})
    exp = Experiment('E', nsk.Condition(t_end=2.5), data)
    prob = CalibrationProblem(_model(BLOW_TOY), [exp],
                              [Observable('x', 'x', sigma=0.05)],
                              [FitParameter('k', (0.1, 1), initial=0.4)])
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert math.isclose(fit.estimates['k'], 0.5, rel_tol=1e-6)
    assert fit.n_failed_evals == 0      # the fit's runs end at 0.9
    with pytest.warns(UserWarning, match='could not be computed'):
        fig, axes = fit.plot()
    lines = axes[0, 0].get_lines()
    model = [line for line in lines if line.get_label() == 'model'][0]
    assert np.array_equal(model.get_xdata(), t)
    assert np.allclose(model.get_ydata(), fit.residuals['model'])
    plt.close('all')


# --- (j) datasets and experiments --------------------------------------------
def test_dataset_from_long_csv(tmp_path):
    path = tmp_path / 'hplc.csv'
    path.write_text('run,t_h,analyte,g_L,sd\n'
                    'FB1,0.5,S,9.1,0.2\n'
                    'FB1,1.0,S,8.0,\n'
                    '2,0.5,P,1.0,0.1\n')
    d = Dataset.from_long_csv(path, experiment='run', time='t_h',
                              observable='analyte', value='g_L', sigma='sd')
    assert len(d) == 3
    assert d.experiments == ['FB1', '2'] and d.observables == ['S', 'P']
    f = d.frame
    assert list(f.columns) == list(cal._DATA_COLUMNS)
    assert np.allclose(f['time'], [0.5, 1.0, 0.5])
    assert np.isnan(f['sigma'][1]) and f['sigma'][0] == 0.2
    sub = d['FB1']
    assert len(sub) == 2 and sub.experiments == ['FB1']
    assert '2' in d and 'FB9' not in d
    with pytest.raises(KeyError, match='FB9'):
        d['FB9']
    with pytest.raises(ValueError, match=r"\['nope'\] not found"):
        Dataset.from_long_csv(path, experiment='run', time='t_h',
                              observable='analyte', value='g_L',
                              sigma='nope')
    plain = Dataset.from_long_csv(path, experiment='run', time='t_h',
                                  observable='analyte', value='g_L')
    assert plain.frame['sigma'].isna().all()
    # labels are stored as str; numeric keys are converted
    assert 2 in d and d[2] == d['2']
    # a file with the canonical column names needs no keywords
    canonical = tmp_path / 'canonical.csv'
    d.frame.to_csv(canonical, index=False)
    assert Dataset.from_long_csv(canonical) == d


def test_dataset_from_frame_uses_a_sigma_column_like_the_constructor():
    df = pd.DataFrame({'experiment': [1, 1], 'observable': 'S',
                       'time': [1.0, 2.0], 'value': [5.0, 4.0],
                       'sigma': [0.1, np.nan], 'sd': [0.3, 0.4]})
    assert Dataset.from_frame(df) == Dataset(df)
    assert Dataset.from_frame(df).frame['sigma'][0] == 0.1
    assert np.allclose(Dataset.from_frame(df, sigma='sd').frame['sigma'],
                       [0.3, 0.4])
    assert Dataset.from_frame(df.drop(columns='sigma')).frame[
        'sigma'].isna().all()
    renamed = df.rename(columns={'experiment': 'run', 'time': 't_h'})
    assert Dataset.from_frame(renamed, experiment='run', time='t_h') == \
        Dataset(df)
    ds = Dataset(df)
    assert 1 in ds and len(ds[1]) == 2 and ds[1] == ds['1']
    with pytest.raises(KeyError, match="No experiment '2'"):
        ds[2]


def test_dataset_from_wide_drops_empty_cells():
    w = pd.DataFrame({'exp': ['A', 'A', 'B'], 't': [1, 2, 1],
                      'S': [1.0, np.nan, 3.0], 'P': [0.1, 0.2, 0.3],
                      'sdP': [0.01, 0.02, 0.03]})
    d = Dataset.from_wide(w, experiment='exp', time='t',
                          columns={'S': 'substrate', 'P': 'product'},
                          sigma={'P': 'sdP'})
    f = d.frame
    assert list(f['observable']) == ['substrate'] * 2 + ['product'] * 3
    assert list(f['experiment']) == ['A', 'B', 'A', 'A', 'B']
    assert np.allclose(f['value'], [1.0, 3.0, 0.1, 0.2, 0.3])
    assert np.allclose(f['sigma'][2:], [0.01, 0.02, 0.03])
    assert f['sigma'][:2].isna().all()
    d2 = Dataset.from_wide(w, experiment='exp', time='t', columns=['P'])
    assert d2.observables == ['P'] and len(d2) == 3
    with pytest.raises(ValueError, match='not among the measured'):
        Dataset.from_wide(w, experiment='exp', time='t', columns=['S'],
                          sigma={'P': 'sdP'})


@pytest.mark.parametrize('change, match', [
    ({'value': [1.0, np.nan]}, r'value is not finite in rows \[1\]'),
    ({'time': [0.0, np.inf]}, 'time is not finite'),
    ({'sigma': [0.1, 0.0]}, 'sigma must be finite and > 0'),
    ({'sigma': [-0.1, np.nan]}, 'sigma must be finite and > 0'),
    ({'experiment': ['E', None]}, 'experiment labels are missing'),
    ({'observable': ['S', '']}, 'observable labels are empty'),
    ({'time': ['a', 'b']}, 'time must be numeric'),
])
def test_dataset_validation(change, match):
    base = {'experiment': ['E', 'E'], 'observable': ['S', 'S'],
            'time': [0.0, 1.0], 'value': [1.0, 2.0], 'sigma': [0.1, np.nan]}
    base.update(change)
    with pytest.raises(ValueError, match=match):
        Dataset(pd.DataFrame(base))
    with pytest.raises(ValueError, match='missing column'):
        Dataset(pd.DataFrame({'experiment': ['E'], 'time': [1.0]}))
    with pytest.raises(ValueError, match='no rows'):
        Dataset(pd.DataFrame(columns=['experiment', 'observable', 'time',
                                      'value']))
    with pytest.raises(TypeError, match='DataFrame'):
        Dataset({'experiment': ['E']})


def test_experiment_output_grid_and_data_selection():
    rows = pd.DataFrame({'experiment': ['A', 'A', 'A', 'B'],
                         'observable': ['S', 'S', 'P', 'S'],
                         'time': [2.0, 1.0, 2.0, 3.0],
                         'value': [1.0, 2.0, 3.0, 4.0]})
    data = Dataset(rows)
    exp = Experiment('A', nsk.Condition(t0=0.5, t_end=5), data)
    assert exp.times.tolist() == [0.5, 1.0, 2.0]       # {t0} U data times
    assert len(exp.data) == 3 and exp.data.experiments == ['A']
    # a one-experiment dataset is relabelled
    relabelled = Experiment('X', nsk.Condition(t_end=5), data['B'])
    assert relabelled.data.experiments == ['X']
    with pytest.raises(ValueError, match="none named 'C'"):
        Experiment('C', nsk.Condition(t_end=5), data)
    with pytest.raises(ValueError, match=r'outside the simulated span'):
        Experiment('A', nsk.Condition(t0=1.5, t_end=5), data)
    with pytest.raises(ValueError, match='at least one later time'):
        Experiment('B', nsk.Condition(t0=3.0, t_end=5), data)
    with pytest.raises(ValueError, match='listed more than once'):
        Experiment('A', nsk.Condition(t_end=5), data, observables=['S', 'S'])
    # times= conditions span times[0]..times[-1]
    with pytest.raises(ValueError, match=r'\[0, 2\.5\]'):
        Experiment('B', nsk.Condition(times=[0.0, 2.5]), data)


def test_experiment_observables_override_the_shared_ones():
    data = _decay_data()
    data['value'] = 10 - data['value']          # P = S0 - S
    data['observable'] = 'X'
    km = _model(DECAY_TOY)
    exp = Experiment('E', nsk.Condition(t_end=6), data,
                     observables=[Observable('X', '[P]', sigma=0.1)])
    prob = CalibrationProblem(km, [exp], [Observable('X', '[S]', sigma=0.1)],
                              [FitParameter('k', (0.01, 10))])
    assert np.allclose(prob.residuals({'k': 0.5}), 0.0, atol=1e-6)
    with pytest.raises(ValueError, match="names observable 'Y'"):
        CalibrationProblem(km, [Experiment('E', nsk.Condition(t_end=6), data,
                                           observables=['Y'])],
                           [Observable('X', '[S]', sigma=0.1)],
                           [FitParameter('k', (0.01, 10))])
    # an experiment's list is its complete set: overriding only S's sigma
    # must also name the shared P (the documented example)
    both = pd.concat([_decay_data(),
                      _decay_data().assign(observable='P')])
    shared = [Observable('S', '[S]', sigma=0.1),
              Observable('P', '[P]', sigma=0.1)]
    with pytest.raises(ValueError, match=r"\['P'\] that it does not define"
                       r".*complete set; name the shared observables"):
        CalibrationProblem(km, [Experiment(
            'E', nsk.Condition(t_end=6), both,
            observables=[Observable('S', '[S]', sigma=0.5)])], shared,
            [FitParameter('k', (0.01, 10))])
    prob = CalibrationProblem(km, [Experiment(
        'E', nsk.Condition(t_end=6), both,
        observables=['P', Observable('S', '[S]', sigma=0.5)])], shared,
        [FitParameter('k', (0.01, 10))])
    sigma = prob.data.groupby('observable')['sigma'].unique()
    assert list(sigma['S']) == [0.5] and list(sigma['P']) == [0.1]
